#include <stdlib.h>
#include "svars.h"
#ifndef CONSTKGRID
#error Must use CONSTKGRID
#endif
//!-----------------------------------------------
//!计算波高周期波向等
//!-----------------------------------------------
#define  IWFDK   0
#define  IIWF2DK 1
#define  IWS2DK  2
#define  IWF     3
#define  IMDPS   4
#define  BVIWS2  0
#define  BVIMDPS 1
#define LDO(P,l) svld1_f64(pm,P+kpos[l*8])
#define LDWK(O) DUP(lpp->wkal[O])
#define FABS(v) svabs_f64_x(pm,v)
//#define DBGINF //printf("%s %4d %04d ",__FILE__,__LINE__,pisp.mpi_id)
#define DBGINFE //printf("%s %4d %04d\n",__FILE__,__LINE__,pisp.mpi_id)
extern "C" {
  void initc_output_cal_(int *nwpc_,int *ndep_,double*pvkd_,double*pebdep_);
  void mean1_(HTHREADINFO ti,double *ea,int *iacb,int *iace, float *ape ,float *tpf ,float *aet ,float *h1_3);
  void mixture_(HTHREADINFO ti,double *ea,int *iacb,int *iace,float*bv);
  void mean1t_(threadGroup *gi,double *ea,int *iac_,float*apet,float*tpft,float*aett,float*h1_3t);
  void initc_output_cal_(int *nwpc_,int *ndep_,double*pvkd_,double*pebdep_){
    md.nwpc=*nwpc_;
    md.pvkd=pvkd_;
    md.pebdep=pebdep_;
    md.ndep =*ndep_;
    md.npdep=BVIMDPS+md.ndep;
  }
  void g_initOutput_cal(threadGroup *gi){
    gi->ndep=md.ndep;
    gi->npdep=md.npdep;
    gi->pebdep=HNEWN(double,(md.nwpc+1)*gi->npdep*kld);
    gi->pvkd  =HNEWN(double,(md.nwpc+1)*IMDPS*kl);
    //write through ape,tpf,aet,h1_3,bv; 
    for(int iac=1;iac<=gi->nwpc;iac++){
      double *d,*s;
      d=gi->pvkd+iac*IMDPS*kl;
      s=md.pvkd+gi->l2gind[iac]*IMDPS*kl;
      for(int j=0;j<IMDPS*kl;j+=8){
        VVD(d+j) =VVD(s+j);
      }
      d=gi->pebdep+iac*gi->npdep*kl;
      s=md.pebdep+gi->l2gind[iac]*gi->npdep*kl;
      for(int j=0;j<gi->npdep*kl;j+=8){
        VVD(d+j) =VVD(s+j);
      }
    }
  }

  void mean1_1(threadGroup *gi,double *eei,double*pvkdt, float *ape ,float *tpf ,float *aet ,float *h1_3){
    int iac;
    double wfk ,wfk1,wsk ,wsk1,wkk ,wkk1;
    ImplschPar *lpp=&gi->pisp;
    double aess=0,awfss=0,asiss=0,apess=0,aets=0,aetc=0;
    VDOUBLE vaess=DUP(0.),vawfss=DUP(0.),vasiss=DUP(0.),vapess=DUP(0.),vaets=DUP(0.),vaetc=DUP(0.);
    int i,k;
    double tztp=1.2,tztz=1.099314;
    //double pi=3.1415926535897932384626433832795d0;
    double zpi=2*3.1415926535897932384626433832795;
    double pi2d=360/zpi;
    svbool_t pm=svptrue_b64();
    for(  k=0;k<kl ;k+=NKP){
      double ekjs=0, aetst=0.,aetct=0.;
      VDOUBLE vekjs, vaetst, vaetct, vwkdkk;
      vekjs=vaetst=vaetct=vwkdkk=DUP(0.);
      double wkdkk;
      int kj,j;
      for(j=0,kj=k;j<jnthet;kj+=kl,j++){
        VDOUBLE vekj=VVD(eei+kj);
        vekjs+=vekj;
        vaetst+=vekj*VVD(lpp->sinths+j);
        vaetct+=vekj*VVD(lpp->cosths+j);
        //ekjs [i]+=ekj;
        //aetst[i]+=ekj*lpp->sinths[j];
        //aetct[i]+=ekj*lpp->cosths[j];
      }
      vwkdkk=VVD(lpp->wkdk+k);
      vaess +=vekjs  *VVD(lpp->dwk+k);
      vapess+=vekjs  *VVD((pvkdt+kl*IWS2DK )+k);
      vasiss+=vekjs  *VVD((pvkdt+kl*IIWF2DK)+k);
      vawfss+=vekjs  *VVD((pvkdt+kl*IWFDK  )+k);
      vaets +=vaetst *vwkdkk;
      vaetc +=vaetct *vwkdkk;
    }
    aets=VSUM(vaets);
    aetc=VSUM(vaetc);
    awfss=VSUM(vawfss);
    asiss=VSUM(vasiss);
    apess=VSUM(vapess);
    aess=VSUM(vaess);
    //!^^^^^^^^^^^^^^^^^^
    //!ape: tz; tpf: tp; aet: th; h1_3:hs
    *aet=atan2(aets,aetc)*pi2d;
    if (*aet<0.)*aet=360.+*aet;
    *h1_3=4.*sqrt(aess);
    if(aess>1e-30){
      *ape=tztz*zpi*sqrt(aess/apess);
      *tpf=asiss*awfss/(aess*aess);
    }else{
      *ape=0;
      *tpf=0;
    }
  }
  void mean1_(HTHREADINFO ti,double *ea,int *iacb,int *iace, float *ape ,float *tpf ,float *aet ,float *h1_3){
    int iac;
    threadGroup *gi=ti->pg;
    for(iac=*iacb;iac<*iace;iac++){
      int iacg=gi->l2gind[iac];
      mean1_1(gi,ea+iac*mkj,gi->pvkd+iac*IMDPS*kld,
              ape +iacg,
              tpf +iacg,
              aet +iacg,
              h1_3+iacg);
    }
  }
#if 1
  void mean1t_(threadGroup *gi,double *ea,int *iac_,float*apet,float*tpft,float*aett,float*h1_3t){
    int iac=*iac_;
    mean1_1(gi,ea+iac*mkj,md.pvkd+iac*IMDPS*kld,apet,tpft,aett,h1_3t);
  }
#endif
  void mixture_1(threadGroup *gi,double *eei,double *pebdept,float *bvo,int nskip){
    double wkk,wsk2;
    int kh,k,j;
    ImplschPar *lpp=&gi->pisp;
    svbool_t pm=svptrue_b64();
    //printf("mixtrue_1 aaaaaaaaaaaaaa\n");
    for(kh=0;kh<gi->ndep;kh++){
      double bv1,bv2,bv3,bvt,ebdep;
      bv1=0;bv2=0;bv3=0;
      VDOUBLE vbvt=DUP(0.),vbv1=DUP(0.),vbv2=DUP(0.),vbv3=DUP(0.);
      VDOUBLE vwkk=DUP(0.),vwsk2=DUP(0.),vebdep=DUP(0.);
      for(  k=0;k<kl ;k+=8){
        double ekjs=0;
        VDOUBLE vekjs=DUP(0.0);
        int kj;
        for(kj=k;kj<mkj;kj+=kl){
          vekjs +=VVD(eei+kj);
        }
        vwkk=VVD(lpp->wk+k);
        vwsk2=VVD(pebdept+k+BVIWS2*kl);
        vebdep=VVD(pebdept+k+(BVIMDPS+kh)*kl);
        vbvt=vekjs*vebdep  ;vbv1+=vbvt;
        vbvt*=vwsk2    ;vbv2+=vbvt;
        vbvt*=vwkk     ;vbv3+=vbvt;
      }
      bv1=VSUM(vbv1);
      bv2=VSUM(vbv2);
      bv3=VSUM(vbv3);
      if(bv2>1e-30){
        bvo[kh*nskip]=bv1/sqrt(bv2)*bv3;
      }else{
        bvo[kh*nskip]=0;
      }
      //if(kh*nskip<5 && bvo[kh*nskip]!=0) printf("bvo[kh*nskip] %f \n",bvo[kh*nskip]);
    }
    //exit(3);
  }
  void mixture_(HTHREADINFO ti,double *ea,int *iacb,int *iace,float*bv){
    int iac;
    threadGroup *gi=ti->pg;
    for(iac=*iacb;iac<*iace;iac++){
      int iacg=gi->l2gind[iac];
      mixture_1(gi,ea+iac*mkj,
                gi->pebdep+iac*gi->npdep*kld,
                md.bv+iacg,md.nwpc+1);
    }
  }
  void addea_(double*ea,double*ee,int*iacb,int*iace){
    int iac,kj;
    for(iac=*iacb;iac<*iace;iac++){
      double *pea=ea+mkj*iac;
      double *pee=ea+mkj*iac;
      for(kj=0;kj<mkj;kj+=8){
        VVD(pea+kj)+=VVD(pee+kj);
      }
    }
  }
  void averageea_(double*ea,int*nea,int*iacb,int*iace){
    int iac,kj;
    double rv=1./ (*nea);
    VDOUBLE vrv=DUP(rv);
    for(iac=*iacb;iac<*iace;iac++){
      double *pe=ea+mkj*iac;
      for(kj=0;kj<mkj;kj+=8)VVD(pe+kj)*=vrv;
    }
  }
  void zeroea_(double*ea,int*iacb,int*iace){
    memset(ea,0,kl*jnthet*(*iace-*iacb+1)*sizeof(double));
  }

}
