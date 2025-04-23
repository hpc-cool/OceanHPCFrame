#include "svars.h"
ImplschPar pisp;
prgstate pstates;
//====================================================
void  initmtpar(){
  md.pis=gppar.pis;
  md.ipos12=gppar.ipos12;//nwpc
  md.ipos8=gppar.ipos8;//nwpc
}
void C_InitImplsch(ImplschPar *pisp_,ImplschP*pip,float*pwvs){
  int isize;
  isize=(long)&(((ImplschPar *)0)->sds);
  memcpy((void*)&pisp,(void*)pisp_,isize );//for align so copy it
  gppar.ip=pip;
  pisp.pip=pip; //gppar.ip;
  //pisp.pwvs=(windvs*)pwvs;
  pisp.pis=gppar.pis;
  pisp.ipos8=gppar.ipos8;
  pisp.ipos12=gppar.ipos12;
  pisp.Dm2=-2;pisp.D1=1.;pisp.enh=1.;
  pisp.zpi=3.1415926535897932384626433832795*2.;
  pisp.D1m30=1e-30;
  pisp.F44_24778=44.24778;
  pisp.F16_0=16.;
  pisp.F0_467=0.467;//mean3
  pisp.F18_20832=18.20832;
  pisp.F0_001930368=0.001930368;
  pisp.F0_5=0.5;
  pisp.F1_0       = 1        ;
  pisp.F0_75      = 0.75     ;
  pisp.F5_5       = 5.5      ;
  pisp.F0_833     = 0.833    ;
  pisp.Fm1_25     =-1.25     ;
  pisp.F28_       = 28       ;
  pisp.F0_80      = 0.80     ;
  pisp.F0_065     = 0.065    ;
  pisp.F0_001     = 0.001    ;
  pisp.tztp       = 1.2      ;
  pisp.tzzpi      = pisp.zpi*1.099314;//zpi*1.099314;
  InitCImplsch();
}

void C_SetPropInterg(int *iac_,int*j_,int*k_,int*iquad_,int*kpos,float*pab,float*pabp){
  p_interg *pis=gppar.pis+ ((*iac_)*mkj+(*j_)*kl)+*k_;
  printf("C_SetPropInterg N %d %d %d\n",*iac_,*j_,*k_);
  pis->iquad  =*iquad_;
  pis->kpos[0]=kpos[0];pis->kpos[1]=kpos[1];pis->kpos[2]=kpos[2];pis->kpos[3]=kpos[3];
  pis->pab [0]=pab [0];pis->pab [1]=pab [1];pis->pab [2]=pab [2];pis->pab [3]=pab [3];
  pis->pabp[0]=pabp[0];pis->pabp[1]=pabp[1];pis->pabp[2]=pabp[2];pis->pabp[3]=pabp[3];
}
void C_Setwind(){
  if(gppar.uconstwind==0){
    ImplschP*pip=pisp.pip;
    int     nwpc=gppar.nwpc;
    windvs * wxy=gppar.pwvs;
    int i;
    for(i=nwpc;i;i--)pip[i].wvs  =wxy  [i];
  }else{
    static int first=1;
    if(first){
      if(gppar.uconstwind==2){
        ImplschP*pip=pisp.pip;
        int     nwpc=gppar.nwpc;
        int i;
        for(i=nwpc;i;i--){
          pip[i].wvs  =gppar.constwind;
        }
      }else{
        ImplschP*pip=pisp.pip;
        int     nwpc=gppar.nwpc;
        windvs * wxy=gppar.pwvs;
        int i;
        for(i=nwpc;i;i--){
          pip[i].wvs  =wxy  [i];
        }
      }
      first=0;
    }
  }
}

int CheckStop(){
  return md.stop_now;
}
void caliacbe(HTHREADINFO ti){
  int nt,n1,n2,nm,nwpwt,nn1,nwps,nwpc;
  float sc,ntsc;
  int id=ti->ind-md.NSpecial;
  nt=md.Nthreads-md.NSpecial;
  sc=0.3;ntsc=nt-sc;
  nwpc=gppar.nwpc;
  nwps=gppar.nwps;
  n1=(nwpc-nwps)/ntsc;
  ti->iacb1=nwpc-(nt-id)*n1+1  ;
  ti->iace1=nwpc-(nt-id-1)*n1  ;
  if(ti->iacb1<=nwps)ti->iacb1=nwps+1;
  if(ti->iace1<ti->iacb1)ti->iace1=ti->iacb1-1;
  nwps=nwpc-n1*nt;
  if(nwps<gppar.nwps) nwps=gppar.nwps;
  n2=(nwps+ntsc-1)/ntsc;
  ti->iacb2=nwps-(nt-id)*n2+1  ;
  ti->iace2=nwps-(nt-id-1)*n2  ;
  if(ti->iacb2<=0)ti->iacb2=1;
  if(ti->iace2<ti->iacb2)ti->iace2=ti->iacb2-1;
  printf("AVE %2.2d %2.2d:%4d %4d %4d %4d:%4d %4d: %4d %4d %4d %4d %4d %4d %f\n",
         ti->igrp,ti->ind,
         ti->iacb1,ti->iace1,ti->iacb2,ti->iace2,
         ti->iace1-ti->iacb1+1,ti->iace2-ti->iacb2+1,
         nwps,nwpc,n1,n2,nt,id,ntsc
         );
  fflush(stdout);
}
double dclktime();
double dclktimeinit();
#ifdef DTIMES
#define BT(i)  tbs[i  ] =dclktime();//DBGLF;
#define ET(i)  tds[i  ]+=dclktime()-tbs[i];  //DBGLF;
#define EBT(i) tbs[i+1] =dclktime();tds[i]+=tbs[i+1]-tbs[i];
#else
#define BT(i)
#define ET(i)
#define EBT(i)
#endif
#define SSM  sSetState   (ti,RFB) /*set sub  ready */
#define SWM  sWaitState  (ti,RFB) /*wait grp ready */
#define SWMR sWaitStater (ti,RFB) /*wait grp ready */

#define MWSR mWaitSubsr  (   RFB) /*mmt */
#define MWS  mWaitSubs   (   RFB) /*mmt */
#define MSS  mSetSubs    (   RFB) /*mmt */
void _threadMain_(HTHREADINFO ti){
  int ind=ti->ind;
  double tb,td;
  int RFB=1;
  double tbs[40],tds[40];
  memset(tds,0,sizeof(tds));
  dclktimeinit();
  bindcpu((md.mpi_id%NProcPNode)*NCorePClu+ti->ind);
  //ImplschLV*plv=NULL;
  //DBGLF;
  ti->mstate=&md.state;
  caliacbe(ti);	  	
  //DBGLF;
  RFB=1;SSM;SWM;
  setspec   (gppar.wav_et,1,md.nwpc,1)  ;
  //DBGLF;
  tb=dclktime();
  ti->nstep=-1;
  RFB=10;SSM;SWM;
  for(;;){
    ti->nstep++;
    if(RFB>100000){
      RFB=10;SSM;SWMR;
    }
    RFB++; SSM; //1
    EBT( 4);    propagats (gppar.wav_et,gppar.wav_ec,ti->iacb1,ti->iace1 );
    //wait wind ok
    SWM;
    EBT( 7);    implschs (gppar.wav_ec,gppar.wav_et,ti->iacb1,ti->iace1 );
    EBT(12);    setspec  (gppar.wav_et,ti->iacb1,ti->iace1,2)  ;
    EBT(12);    accumea_ (&ind,gppar.wav_et,&ti->iacb1,&ti->iace1);
    //wait boudary ok
    RFB++;SSM;SWM; //2
    EBT(10);    propagats(gppar.wav_et,gppar.wav_ec,ti->iacb2,ti->iace2);
    EBT(11);    implschs (gppar.wav_ec,gppar.wav_et,ti->iacb2,ti->iace2);
    EBT(12);    setspec  (gppar.wav_et,ti->iacb2,ti->iace2,2);
    EBT(12);    accumea_ (&ind,gppar.wav_et,&ti->iacb2,&ti->iace2);
    if(md.hist_eot){
      //checkoutput_(ti,&ind,&RFB,gppar.wav_et,&ti->iacb1,&ti->iace1,&ti->iacb2,&ti->iace2,&md.hist_eot);
    }
    if(CheckStop())break;
  }
  RFB+=1000;SSM;usleep(100);SWM;
}
