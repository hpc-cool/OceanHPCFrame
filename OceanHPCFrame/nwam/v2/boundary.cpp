#include <stdlib.h>
#include <stdio.h>
#include "ksvml.h"
#include "svars.h"
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
Boundary bdy;
__BEGIN_DECLS
void initc_boundary_(FBoundary *bdy_,double*pvws,int *nsp,float*wxy);
void setspec_(double *eet,int *iacb,int *iace,int *type) ;
__END_DECLS

void initc_boundary_(FBoundary *bdy_,double*pvws,int *nsp,float*wxy){
  int j;
  for(j=0;j<CJNTHET;j++)bdy.cosths[j]=bdy_->cosths[j];
  for(j=0;j<CJNTHET;j++)bdy.sinths[j]=bdy_->sinths[j];
  for(j=0;j<CKL    ;j++)bdy.rwk[j]=1./bdy_->wk[j];
  bdy.windfield=bdy_->windfield;
  bdy.nwpc=bdy_->nwpc;
  bdy.g=9.81;
  bdy.rg=1./9.81;
  bdy.gama=3.3;
  bdy.zpi=2*3.1415926535897932384626433832795;
  bdy.xj0=bdy.windfield*1000.;
  bdy.xj_=bdy.g*bdy.xj0;
  bdy.arlfa_=0.076*2/bdy.zpi;
  bdy.wsj_=22.*bdy.g;
  bdy.F_1=1.;
  bdy.F_m1d25=-1.25;
  bdy.F_m0d5=-0.5;
  bdy.rsigma1=1/0.07;
  bdy.rsigma2=1./0.09;
  bdy.F_m0d4=-0.4;
  bdy.F_m0d33=-0.33;
  bdy.wv0=0.9;
  bdy.pvws=pvws;
  bwxy=wxy;
}
__attribute__ ((optnone)) void calwl(DOUBLE vx,DOUBLE vy,DOUBLE rwv,double *wl){
  VDOUBLE vvx,vvy,vvt,vc,vs,vrwv;
  vvx=DUP(vx);            vvy=DUP(vy); 
  vrwv=DUP(rwv);          vrwv=vrwv*vrwv;
 
  vc=VVD(bdy.cosths);     vs=VVD(bdy.sinths);
  vvt=(vc*vvx+vs*vvy);    
  VVD(wl)=vvt*vvt*vrwv;
 
  vc=VVD(bdy.cosths+NKP); vs=VVD(bdy.sinths+NKP);
  vvt=(vc*vvx+vs*vvy);    
  VVD(wl+NKP)=vvt*vvt*vrwv;
}
//风浪谱，n=1 全场初始化，n=2 边界
void setspec(double *e,int nwpb,int nwpe,int n) {
  int iac;
  Boundary *pbdy=&bdy;

  for(iac=nwpb;iac<nwpe;iac++){
    if(gppar.nsp[iac]<n)continue;
    float*wxy=bwxy+iac*4;
    double vx=wxy[0];
    double vy=wxy[1];
    double wv=wxy[2];
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
    calwl(vx,vy,rwv,wl);
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
      vvt=EXP(pbdy->F_m1d25*vvt*vvt);
      vrwk0=vrwk0*vrwk0*vrwk0*vrwk0;
      vvt1=vrsigma-vrsigma*vwsk*vrwsj;
      vvt1=EXP(pbdy->F_m0d5*vvt1*vvt1);
      va1=DUP(arlfa)*(vrwk0)*vvt*POW(DUP(pbdy->gama),vvt1);
      for (j=0,kj=k;j<jnthet;j++,kj+=kl){
        VVD(eei+kj)=DUP(wl[j])*va1;
      } 
    }
  }
}
void setspec_(double *e,int *nwpb,int*nwpe,int*n) {
  setspec(e,*nwpb,*nwpe,*n) ;
}
