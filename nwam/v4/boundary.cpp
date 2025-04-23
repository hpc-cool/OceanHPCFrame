#include <stdlib.h>
#include <stdio.h>
#include "svars.h"
#include "ksvml.h"
#define VDOUBLE svfloat64_t
#define DOUBLE float64_t
#define ADD(A,B) svadd_f64_m(pm,A,B)
#define MLA(A,B,C) svmla_f64_m(pm,A,B,C)
#define MUL(A,B) svmul_f64_m(pm,A,B)
#define VVD(P) *(svfloat64_t*)((DOUBLE*)(P))
#define VSUM(v) svaddv_f64(pm,v)
#define LD(P) svld1_f64(pm,P)
#define DUP(A) svdup_f64(A)
#define LDO(P,l) svld1_f64(pm,P+kpos[l*8])
#define LDWK(O) DUP(lpp->wkal[O])
#define FABS(v) svabs_f64_x(pm,v)
#define NKP 8
typedef struct _FBoundary{
  double cosths[CJNTHET],sinths[CJNTHET];
  double wk[CKL];
  double windfield;
  int nwpc;
}FBoundary;
float*bwxy;
__BEGIN_DECLS
void initc_boundary_(FBoundary *bdy_,double*pvws,int *nsp,float*wxy);
void setspec_(double *eet,int *iacb,int *iace,int *type) ;
void setspec(threadGroup *gi,double *e,int nwpb,int nwpe,int n) ;
__END_DECLS

void initc_boundary_(FBoundary *bdy_,double*pvws,int *nsp,float*wxy){
  int j;
  for(j=0;j<CJNTHET;j++)md.bdy.cosths[j]=bdy_->cosths[j];
  for(j=0;j<CJNTHET;j++)md.bdy.sinths[j]=bdy_->sinths[j];
  for(j=0;j<CKL    ;j++)md.bdy.rwk[j]=1./bdy_->wk[j];
  md.bdy.windfield=bdy_->windfield;
  md.bdy.nwpc=bdy_->nwpc;
  md.bdy.g=9.81;
  md.bdy.rg=1./9.81;
  md.bdy.gama=3.3;
  md.bdy.zpi=2*3.1415926535897932384626433832795;
  md.bdy.xj0=md.bdy.windfield*1000.;
  md.bdy.xj_=md.bdy.g*md.bdy.xj0;
  md.bdy.arlfa_=0.076*2/md.bdy.zpi;
  md.bdy.wsj_=22.*md.bdy.g;
  md.bdy.F_1=1.;
  md.bdy.F_m1d25=-1.25;
  md.bdy.F_m0d5=-0.5;
  md.bdy.rsigma1=1/0.07;
  md.bdy.rsigma2=1./0.09;
  md.bdy.F_m0d4=-0.4;
  md.bdy.F_m0d33=-0.33;
  md.bdy.wv0=0.9;
  md.bdy.pvws=pvws;
}
void g_initSetspec(threadGroup *gi){
  memcpy(&gi->bdy,&md.bdy,sizeof(gi->bdy));
  gi->bdy.pvws=(double*)HNEWN(double,kl*(gi->nwpc+1));
  gi->bdy.nwpc=gi->nwpc;
  double *d,*s;
  for(int i=1;i<=gi->nwpc;i++){
    int l2g=gi->l2gind[i];
    s=md.bdy.pvws+kl*l2g;
    d=gi->bdy.pvws+kl*i;
    VVD(d+ 0)=VVD(s+ 0);
    VVD(d+ 8)=VVD(s+ 8);
    VVD(d+16)=VVD(s+16);
  }				
}
__attribute__ ((optnone)) void calwl(threadGroup *gi,DOUBLE vx,DOUBLE vy,DOUBLE rwv,double *wl){
  VDOUBLE vvx,vvy,vvt,vc,vs,vrwv;
  vvx=DUP(vx);            vvy=DUP(vy); 
  vrwv=DUP(rwv);          vrwv=vrwv*vrwv;
 
  vc=VVD(gi->bdy.cosths);     vs=VVD(gi->bdy.sinths);
  vvt=(vc*vvx+vs*vvy);    
  VVD(wl)=vvt*vvt*vrwv;
 
  vc=VVD(gi->bdy.cosths+NKP); vs=VVD(gi->bdy.sinths+NKP);
  vvt=(vc*vvx+vs*vvy);    
  VVD(wl+NKP)=vvt*vvt*vrwv;
}
//风浪谱，n=1 全场初始化，n=2 边界
void setspec(threadGroup *gi,double *e,int nwpb,int nwpe,int n) {
  int iac;
  Boundary *pbdy=&gi->bdy;

  for(iac=nwpb;iac<nwpe;iac++){
    if(gi->nsp[iac]<n)continue;
    windvs *wvs=&gi->pisp.pip[iac].wvs;
    double vx=wvs->wx;
    double vy=wvs->wy;
    double wv=wvs->wv;
    double*pvws=pbdy->pvws+iac*kl;

    if (wv<=0) wv=pbdy->wv0;
    double rwv=pbdy->F_1/wv;
    double xj=pbdy->xj_*rwv*rwv;
    double arlfa=pbdy->arlfa_*pow(xj,pbdy->F_m0d4);
    double wsj=pbdy->wsj_*pow(xj,pbdy->F_m0d33)*rwv;
    double wkj=wsj*wsj*pbdy->rg;
    double rwsj=pbdy->F_1/wsj;
    int j,k,kj;
    double alpha,rwk0,wsk,rsigma,vt,vt1,vt2;
    double wl[jnthet+4];
    svbool_t pm=svptrue_b64();
    VDOUBLE vrwsj=DUP(rwsj);
#if 1
    calwl(gi,vx,vy,rwv,wl);
#else
    for (j=0;j<jnthet;j++){
      double vt;
        vt=(vx*pbdy->cosths[j]+vy*pbdy->sinths[j])*rwv;
        wl[j]=vt*vt;
    }
#endif

    double *eei=e+mkj*iac;
    for (k=0;k<kl;k+=8){
      VDOUBLE vrwk0,vwsk,va1,vvt,vvt1,vvt2,vrsigma;
      vrwk0=VVD(pbdy->rwk+k);
      vwsk=VVD(pvws+k);
      svbool_t pq=svcmple(pm,vwsk,DUP(wsj));
      vrsigma=svsel_f64(pq,DUP(pbdy->rsigma1),DUP(pbdy->rsigma2));
      //alpha=arlfa/wk0**4 * exp(-1.25*(wkj/wk0)**2) * gama**(exp(-0.5*((1.-wsk/wsj)/sigma)**2))*(wl/wv)**2
      vvt= DUP(wkj)*vrwk0;
      vvt=svml_exp_f64(pbdy->F_m1d25*vvt*vvt);
      vrwk0=vrwk0*vrwk0*vrwk0*vrwk0;
      vvt1=vrsigma-vrsigma*vwsk*vrwsj;
      vvt1=svml_exp_f64(pbdy->F_m0d5*vvt1*vvt1);
      va1=DUP(arlfa)*(vrwk0)*vvt*svml_pow_f64(DUP(pbdy->gama),vvt1);
      for (j=0,kj=k;j<jnthet;j++,kj+=kl){
        VVD(eei+kj)=DUP(wl[j])*va1;
      } 
    }
  }
}
void setspec_(double *e,int *nwpb,int*nwpe,int*n) {
  setspec(NULL,e,*nwpb,*nwpe,*n) ;
}
