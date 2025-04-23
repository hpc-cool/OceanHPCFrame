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
struct OutPutP{
  double *pvkd,*pebdep;
  int nwpc,ndep,npdep;
}vop={0};
struct OutPutP*op=&vop;
extern "C" {
  void initc_output_cal_(int *nwpc_,int *ndep_,double*pvkd_,double*pebdep_);
  void mean1_(double *ea,int *iacb,int *iace, float *ape ,float *tpf ,float *aet ,float *h1_3);
  void mixture_(double *ea,int *iacb,int *iace,float*bv);
  void mean1t_(double *ea,int *iac_,float*apet,float*tpft,float*aett,float*h1_3t);
  void initc_output_cal_(int *nwpc_,int *ndep_,double*pvkd_,double*pebdep_){
    op->nwpc=*nwpc_;
    //op->pvkd=HNEWN(double,(op->nwpc+1)*IMDPS*kl);
    //memcpy(op->pvkd,pvkd_,sizeof(double)*(op->nwpc+1)*IMDPS*kl);
    op->pvkd=pvkd_;
    op->ndep =*ndep_;
    op->npdep=BVIMDPS+op->ndep;
    op->pebdep=pebdep_;
  }
  void mean1_1(double *eei,double*pvkdt, float *ape ,float *tpf ,float *aet ,float *h1_3){
    double wfk ,wfk1,wsk ,wsk1,wkk ,wkk1;
    ImplschPar *lpp=&pisp;
    double aess=0,awfss=0,asiss=0,apess=0,aets=0,aetc=0;
    VDOUBLE vaess=DUP(0.),vawfss=DUP(0.),vasiss=DUP(0.),vapess=DUP(0.),vaets=DUP(0.),vaetc=DUP(0.);
    int k;
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
  void mean1_(double *ea,int *iacb,int *iace, float *ape ,float *tpf ,float *aet ,float *h1_3){
    int iac;
    for(iac=*iacb;iac<*iace;iac++){
      mean1_1(ea+iac*mkj,op->pvkd+iac*IMDPS*kld,ape+iac,tpf+iac,aet+iac,h1_3+iac);
    }
  }
  void mean1t_(double *ea,int *iac_,float*apet,float*tpft,float*aett,float*h1_3t){
    int iac=*iac_;
    mean1_1(ea+iac*mkj,op->pvkd+iac*IMDPS*kld,apet,tpft,aett,h1_3t);
  }

  void mixture_1(double *eei,double *pebdept,float *bvo,int nskip){
    double wkk,wsk2;
    int kh,k,j;
    ImplschPar *lpp=&pisp;
    svbool_t pm=svptrue_b64();
    for(kh=0;kh<op->ndep;kh++){
      double bv1,bv2,bv3,bvt,ebdep;
      bv1=0;bv2=0;bv3=0;
      VDOUBLE vbvt,vbv1,vbv2,vbv3,vwkk,vwsk2,vebdep;
      vbvt=vbv1=vbv2=vbv3=vwkk=vwsk2=vebdep=DUP(0.);
      for(  k=0;k<kl ;k+=8){
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
    }
  }
  void mixture_(double *ea,int *iacb,int *iace,float*bv){
    int iac;
    for(iac=*iacb;iac<*iace;iac++){
      mixture_1(ea+iac*mkj,op->pebdep+iac*op->npdep*kld,bv+iac,op->nwpc+1);
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
  void zeroea(double*ea,int iacb,int iace){
    memset(ea,0,kl*jnthet*(iace-iacb+1)*sizeof(double));
  }

}
