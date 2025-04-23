#include "svars.h"
#define NODELTEE
#define LOGID 0
#define DBGINFN printf("%s %d %d %d:",__FILE__,__LINE__,itt,iat);
void Checkee(double*ee,int iac,const char*fn,int ln);
#define CHKEE(ee,iat) //Checkee(ee,iat,__FILE__,__LINE__)
#define CNAN(v,i1,i2,i3,i4) if(isnan(v)){printf("NAN %s %d %d %d %d:%s %d\n",#v ,i1,i2,i3,i4,__FILE__,__LINE__);}
void mean2(ImplschPar *lpp,double*ee,ImplschP*pi,double *sds,double*awk) ;
void mean3(ImplschPar *lpp,double*ee,ImplschP*pi) ;
void InitCImplsch(){

}
#define register 
#define NKP 8
#define KLP (CKL+NKP)
#if 0
#define  PWK(I0)   \
    { wka00=lpp->wkal[I0*4  ];wka01=lpp->wkal[I0*4+1];wka02=lpp->wkal[I0*4+2];wka03=lpp->wkal[I0*4+3];\
      wka10=lpp->wkal[I0*4+4];wka11=lpp->wkal[I0*4+5];wka12=lpp->wkal[I0*4+6];wka13=lpp->wkal[I0*4+7];\
      register int lt0,lt1,ii,kt,jt; ZBYTE *kjs=lpp->kjs+j*48+I0*4;\
      kt=kjs[0];jt=kjs[1]; lt0=jt*KLP+kt;\
      kt=kjs[2];jt=kjs[3]; lt1=jt*KLP+kt;\
      pe00=tee+NKP+lt0;pe01=tee+NKP+lt0+1;pe02=tee+NKP+lt1;pe03=tee+NKP+lt1+1;\
      kjs+=4;\
      kt=kjs[0];jt=kjs[1]; lt0=jt*KLP+kt;\
      kt=kjs[2];jt=kjs[3]; lt1=jt*KLP+kt;\
      pe10=tee+NKP+lt0;pe11=tee+NKP+lt0+1;pe12=tee+NKP+lt1;pe13=tee+NKP+lt1+1;\
    }
#define INTERE12(coff1,coff2) \
      e1=coff1*( wka00*pe00[k]+wka01*pe01[k]+wka02*pe02[k]+wka03*pe03[k]); \
      e2=coff2*( wka10*pe10[k]+wka11*pe11[k]+wka12*pe12[k]+wka13*pe13[k]);
#define CALC12(I0)                                               \
    PWK(I0);                                                     \
    for(k=0;k<CKL;k++){                                         \
      register double e0,e1,e2;e0=pe0[k]; INTERE12(lmdpd,lmdmd);  \
      pse[k]=pse[k]+LIV_m2*e0*(e0*e1+e0*e2+e1*e2);               \
    }
#define CALC34(I0)                                               \
    PWK(I0);                                                     \
    for(k=0;k<CKL;k++){                                         \
      register double e0,e1,e2;e0=pe0[k]; INTERE12(LIV_1,lmdmd);  \
      pse[k]=pse[k]+lmdpd20*e1*(e1*(lmdpd*e0+e2)+lmdpd*e0*e2);   \
    }
#define CALC56(I0)                                               \
    PWK(I0);                                                     \
    for(k=0;k<CKL;k++){                                         \
      register double e0,e1,e2;e0=pe0[k]; INTERE12(lmdpd,LIV_1);  \
      pse[k]=pse[k]+lmdmd20*e2*(e2*(lmdmd*e0+e1)+lmdmd*e0*e1);   \
    }

void winter(ImplschPar *lpp,double *ee,double *se,double enh){
  register double lmdpd,lmdmd,lmdpd20,lmdmd20,LIV_m2,LIV_1;
  int j;
  double tee[KLP*CJNTHET+NKP]; 
  if(1){
    double *ped=tee,*pes=ee;
    int k,j;
    for(k=0;k<NKP;k++) *ped++=0;
    for(j=0;j<CJNTHET;j++){
      for(k=0;k<CKL;k++) *ped++=*pes++;
      for(k=0;k<NKP;k++) *ped++=0;
    }
  }
  LIV_m2=-2.;LIV_1=1.;
  lmdpd=lpp->lmdpd;lmdmd=lpp->lmdmd;lmdpd20=lpp->lmdpd20;lmdmd20=lpp->lmdmd20;
  memset(se,0,sizeof(double)*mkj);
  for(j=jnthet-1;j>=0;j--){
    register double wka00,wka01,wka02,wka03;
    register double wka10,wka11,wka12,wka13;
    register double*pe00,*pe01,*pe02,*pe03;
    register double*pe10,*pe11,*pe12,*pe13;
    register double*pse;
    double *pe0=tee+NKP+j*KLP;
    register int k;
    pse=se+j*(CKL);
    CALC12( 0);//* resonate 1 // kjoff : -2,-4   (16+8)*kl
    CALC12( 2);//* resonate 2 // kjoff : -2,-4   (16+8)*kl
    CALC34( 4);//* resonate 3 // kjoff : -3,-6   (16+9)*kl
    CALC34( 6);//* resonate 4 // kjoff : -3,-6   (16+9)*kl
    CALC56( 8);//* resonate 5 // kjoff :  3, 5   (16+9)*kl
    CALC56(10);//* resonate 6 // kjoff :  3, 5   (16+9)*kl
    if(1){ //2*kl
      double *cwks17=lpp->cwks17;
      for(k=0;k<CKL;k++){// 9 click 3*4 op
        pse[k]=pse[k]*enh*cwks17[k];
      }
    }
  }
}
#else
#define LDO(P,l) svld1_f64(pm,P+kpos[l*8])
#define LDWK(O) DUP(lpp->wkal[O])
#define FABS(v) svabs_f64_x(pm,v)
#define  PWK(I0){\
    wka00=LDWK(I0*4+0);wka01=LDWK(I0*4+1);wka02=LDWK(I0*4+2);wka03=LDWK(I0*4+3);\
    wka10=LDWK(I0*4+4);wka11=LDWK(I0*4+5);wka12=LDWK(I0*4+6);wka13=LDWK(I0*4+7);\
    register int lt0,lt1,ii,kt,jt; \
    ZBYTE *kjs=lpp->kjs+j*48+I0*4;\
    kt=kjs[0];jt=kjs[1]; lt0=jt*KLP+kt;\
    kt=kjs[2];jt=kjs[3]; lt1=jt*KLP+kt;\
    pe00=tee+NKP+lt0;pe01=tee+NKP+lt0+1;\
    pe02=tee+NKP+lt1;pe03=tee+NKP+lt1+1;\
    kjs+=4;\
    kt=kjs[0];jt=kjs[1]; lt0=jt*KLP+kt;\
    kt=kjs[2];jt=kjs[3]; lt1=jt*KLP+kt;\
    pe10=tee+NKP+lt0;pe11=tee+NKP+lt0+1;\
    pe12=tee+NKP+lt1;pe13=tee+NKP+lt1+1;\
  }
// e1=MUL(coff1,MLA(MLA(MLA(MUL( wka00,LD(pe00+k)),wka01,LD(pe01+k)),wka02,LD(pe02+k)),wka03,LD(pe03+k))); \
// e2=MUL(coff2,MLA(MLA(MLA(MUL( wka10,LD(pe10+k)),wka11,LD(pe11+k)),wka12,LD(pe12+k)),wka13,LD(pe13+k)));
#define INTERE12(coff1,coff2) \
      e1=coff1*( wka00*VVD(pe00+k)+wka01*VVD(pe01+k)+wka02*VVD(pe02+k)+wka03*VVD(pe03+k)); \
      e2=coff2*( wka10*VVD(pe10+k)+wka11*VVD(pe11+k)+wka12*VVD(pe12+k)+wka13*VVD(pe13+k));
#define CALC12(I0)                                                 \
    PWK(I0);                                                       \
    for(k=0;k<CKL;k+=NKP){                                          \
      register VDOUBLE e0,e1,e2;e0=VVD(pe0+k);INTERE12(lmdpd,lmdmd);\
      VVD(pse+k)+=LIV_m2*e0*(e1*e2+e0*(e1+e2));                     \
    }
#define CALC34(I0)                                                  \
    PWK(I0);                                                        \
    for(k=0;k<CKL;k+=NKP){                                          \
      register VDOUBLE e0,e1,e2;e0=VVD(pe0+k);INTERE12(LIV_1,lmdmd);\
      VVD(pse+k)+=lmdpd20*e1*(e1*e2+lmdpd*e0*(e1+e2));              \
    }
#define CALC56(I0)                                                  \
    PWK(I0);                                                        \
    for(k=0;k<CKL;k+=NKP){                                          \
      register VDOUBLE e0,e1,e2;e0=VVD(pe0+k);INTERE12(lmdpd,LIV_1);\
      VVD(pse+k)+=lmdmd20*e2*(e2*e1+lmdmd*e0*(e1+e2));              \
    }

void winter(ImplschPar *lpp,double *ee,double *se,double enh){
  register VDOUBLE lmdpd,lmdmd,lmdpd20,lmdmd20,LIV_m2,LIV_1;
  int j;
  double tee[KLP*CJNTHET+NKP]; 
  svbool_t pm=svptrue_b64();
  if(1){
    double *ped=tee,*pes=ee;
    int k,j;
#if 1
    VVD(tee)=DUP(0.);
    for(j=0;j<CJNTHET;j++){
      VVD(tee+j*(CKL+NKP)+NKP)=VVD(ee+j*CKL);
      VVD(tee+j*(CKL+NKP)+2*NKP)=VVD(ee+j*CKL+NKP);
      VVD(tee+j*(CKL+NKP)+3*NKP)=VVD(ee+j*CKL+2*NKP);
      VVD(tee+j*(CKL+NKP)+4*NKP)=VVD(ee+j*CKL+3*NKP);
      VVD(tee+j*(CKL+NKP)+5*NKP)=DUP(0.);
    }
#else
    for(k=0;k<NKP;k++) *ped++=0;
    for(j=0;j<CJNTHET;j++){
      for(k=0;k<CKL;k++) *ped++=*pes++;
      for(k=0;k<NKP;k++) *ped++=0;
    }
#endif
  }
  LIV_m2=DUP(-2.);LIV_1=DUP(1.);
  lmdpd=DUP(lpp->lmdpd);lmdmd=DUP(lpp->lmdmd);lmdpd20=DUP(lpp->lmdpd20);lmdmd20=DUP(lpp->lmdmd20);
  memset(se,0,sizeof(double)*mkj);
  for(j=jnthet-1;j>=0;j--){
    register VDOUBLE wka00,wka01,wka02,wka03;
    register VDOUBLE wka10,wka11,wka12,wka13;
    register double*pe00,*pe01,*pe02,*pe03;
    register double*pe10,*pe11,*pe12,*pe13;
    register double*pse;
    double *pe0=tee+NKP+j*KLP;
    register int k;
    pse=se+j*(CKL);
    CALC12( 0);//* resonate 1 // kjoff : -2,-4   (16+8)*kl
    CALC12( 2);//* resonate 2 // kjoff : -2,-4   (16+8)*kl
    CALC34( 4);//* resonate 3 // kjoff : -3,-6   (16+9)*kl
    CALC34( 6);//* resonate 4 // kjoff : -3,-6   (16+9)*kl
    CALC56( 8);//* resonate 5 // kjoff :  3, 5   (16+9)*kl
    CALC56(10);//* resonate 6 // kjoff :  3, 5   (16+9)*kl
    if(1){ //2*kl
      double *cwks17=lpp->cwks17;
      VDOUBLE venh=DUP(enh);
      for(k=0;k<CKL;k+=NKP){// 9 click 3*4 op
        VVD(pse+k)*=venh*VVD(cwks17+k);
      }
    }
  }
}
#endif
#undef PWKA
#undef PWKB
#undef INTERE12
#undef CALC12
#undef CALC34
#undef CALC56
void implschs_(double *eec,double *eet,int *piacb,int *piace) {
  implschs(eec,eet,*piacb,*piace) ;
}
void implschs(double *eec,double *eet,int iacb,int iace) {
  int iac;
  int k,j,kj;
  double se[mkj];
  double sds,awk;
  ImplschPar *lpp=&pisp;
  //iacb=1; iace=gppar.nwpc;
  double*ee=eet+mkj*iacb;
  double*ec=eec+mkj*iacb;
  svbool_t pm=svptrue_b64();
  svbool_t pq;
  for(iac=iacb;iac<=iace;iac++,ee+=mkj,ec+=mkj){
    double enh;
    if(gppar.nsp[iac]!=1)continue;
    ImplschP*pi=&lpp->pip[iac];

    mean2(lpp,ee,pi,&sds,&awk);
    //!*************************************
    //! snonlin(e)
    //!*************************************
    {
      double xx=0.75*pi->depth*awk;
      double enh;
      if (xx<0.5)xx=0.5;
      enh=1.+(5.5/xx)*(1.-0.833*xx)*exp(-1.25*xx);
      winter(lpp,ee,se,enh);
    }
    {   
      register double deltts;
      double scd=sqrt((0.80+0.065*pi->wvs.wv)*0.001);
      double bett =lpp->beta10*(1.+lpp->beta11*pi->wvs.wi);
#ifndef NODELTEE
      double wstarm=pi->wvs.wv*scd;
#endif  
      //Wind Input

      deltts=lpp->deltts;
      VDOUBLE vzero=DUP(0.);
      for(j=0;j<jnthet;j++) {
        //!sinput(e)
        double *pe =ee+j*kl;
        double *pc =ec+j*kl;
        double *pse=se+j*kl;
        VDOUBLE betta;
        {
          double wl=pi->wvs.wx*lpp->cosths[j]+pi->wvs.wy*lpp->sinths[j];
          betta=DUP(bett*28.*wl*scd);
        }
        for(k=0;k<kl;k+=NKP){ // vector
          //! sinput(e)
          VDOUBLE wk=VVD(lpp->wk+k);
          VDOUBLE beta=betta*wk-bett*VVD(pi->ws+k);
          //if(beta<0)beta=0;
          //!*************************************
          //!      source terms end
          //!beta wind input
          //!ssds:sdissip
          //!ISSBO:sbottom
          //!cg*cosths[j]*Rsd_tanLat spherical coords
          //! ZZP cg=pi->CCG[k]  ! Big Circle
          VDOUBLE dset=beta-sds*wk+VVD(pi->ssbo+k); //! +pi->CCG[k]*cosths[j]*tRsd_tanLat
          {
            VDOUBLE eev=VVD(pe+k);
            VDOUBLE set=VVD(pse+k)+dset*eev;
            svbool_t pq;
#define NODELTEE
#ifndef NODELTEE
            {
              VDOUBLE deltee=wstarm*VVD(lpp->grolim+k);
              VDOUBLE delteem=-deltee;
              //if(set<-deltee) set=-deltee;
              //else if(set>deltee) set=deltee;
              pq=svcmpgt(pm,set,deltree);
              set=svsel(pq,deltee,set);
              pq=svcmplt(pm,set,deltreem);
              set=svsel(pq,delteem,set);
            }
#endif          
            eev=eev+set*deltts;
            //if(eev<0)eev=0;
            pq=svcmplt(pm,eev,vzero);
            eev=svsel(pq,vzero,eev);
            VVD(pc+k)=eev;
          }
        }
      }
    }
    mean3(lpp,ec,pi);
  }
}

void mean2(ImplschPar *lpp,double*ee,ImplschP*pi,double *psds,double*pawk) {
  double ae,asi,ark,ekspm,awk;
  int k,j,kj;
  VDOUBLE vae,vasi,vark,vawk,vpes;
  vae=vasi=vawk=vark=vpes=DUP(0.);
  for(k=0;k<kl;k+=8){
    vpes=DUP(0.);
    for(j=0,kj=k;j<jnthet;kj+=kl,j++){//vector
      VDOUBLE vpe=VVD(ee+kj);
      vpes+=vpe;
    }
      vae +=vpes*VVD(lpp->dwk+k);
      vawk+=vpes*VVD(lpp->wkdk+k);
      vark+=vpes*VVD(lpp->wkibdk+k);
      vasi+=vpes*VVD(pi->iwsdk+k);
  }
  svbool_t pm=svptrue_b64();
  ae =VSUM(vae ); awk=VSUM(vawk);
  ark=VSUM(vark); asi=VSUM(vasi);
  {
    if(ae>1.e-100){
      awk=awk/ae;
      asi=ae/asi;
      ark=(ae/ark);
      ark=ark*ark;
      //!     sds=2.36e-5*asi*ark**3*ae**2/alpm2
      //!         2.36e-5/alpm2=2.587605807
      ekspm=ae*ark*ark/0.0030162;
      *psds=lpp->ads*lpp->brkd1*asi/ark*sqrt(ekspm)*exp(-lpp->brkd2*0.64/ekspm);
    }else{
      *psds=0;
    }
    *pawk=awk;
  }
}// mean2
void mean3(ImplschPar *lpp,double*ee,ImplschP*pi) {
  int    k,j,kj;
  register double ae,awk;
  register double ae1,awk1;
  VDOUBLE vae,vawk,vev;
  vae=vawk=DUP(0.);
  for(k=0;k<kl;k+=8){
    vev=DUP(0.);
    for(kj=k,j=0;j<jnthet ;kj+=kl,j++){//vector
      vev+=VVD(ee+kj);
    }
    vae +=vev*VVD(lpp->dwk+k);
    vawk+=vev*VVD(lpp->wkdk+k);
  }
  svbool_t pm=svptrue_b64();
  ae =VSUM(vae ); awk=VSUM(vawk); 
  {
    if(ae>1e-30){
      double hs,hb,hbb;
      awk=awk/ae;
      hs=4.*sqrt(ae);
      hb=lpp->zpi/awk*0.142*tanh(pi->depth*awk); //!hb=zpi/awk*0.12*tanh(dep(ia,ic)*awk)/1.6726
      hbb =pi->depth*(0.78125/1.6726)     ; //!hbb=0.78125*dep(ia,ic)/1.5864792
      if(hb>hbb) hb=hbb;
      if(hs>hb) {
        double vt;
        double chbh;
        int kj;
        vt=(hb/hs);chbh=vt;
        chbh=chbh*chbh;
        VDOUBLE vchbh=DUP(chbh);
        for(kj=0;kj<mkj;kj+=8){//vector          
          VVD(ee+kj)*=vchbh;
        }
      }
    }
  }
}
void c_checkee_(double*ee,int *iac,const char*fn,int *ln){
  Checkee(ee,*iac,fn,*ln);
}
extern "C" void wav_abort();
void Checkee(double*ee,int iac,const char*fn,int ln){
  double ae=0,be,me=0;
  int k,j;
  for(j=0;j<jnthet;j++){
    double*pe=ee+j*kl;
    for(  k=0;k<kl ;k++){//vector
      be=fabs(pe[k]);
      if(be>me)me=be;
      ae =ae +be;
    }
  }
  if(isnan(ae)){
    printf("ee is nan %s %d %d \n",fn,ln,iac);
    wav_abort();
  } 
  if(isinf(ae)){
    printf("ee is inf %s %d %d \n",fn,ln,iac);
    wav_abort();
  } 
  if(fabs(me)>1e80){
    printf("ee is big %s %d %d %e\n",fn,ln,iac,me);
    wav_abort();
  } 
}
