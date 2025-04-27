#include "svars.h"
ImplschPar pisp;
prgstate pstates;
//====================================================
void C_InitImplsch(ImplschPar *pisp_,ImplschP*pip,float*pwvs){
  int isize;
  isize=(long)&(((ImplschPar *)0)->sds);
  memcpy((void*)&pisp,(void*)pisp_,isize );//for align so copy it
  gppar.ip=pip;
  pisp.pip=pip; //gppar.ip;
  pisp.pwvs=(windvs*)pwvs;
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
    windvs * wxy=pisp.pwvs;
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
        windvs * wxy=pisp.pwvs;
        int i;
        for(i=nwpc;i;i--){
          pip[i].wvs  =wxy  [i];
        }
      }
      first=0;
    }
  }
}

