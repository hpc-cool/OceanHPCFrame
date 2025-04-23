#include "svars.h"
prgstate pstates;
/*
 *iepos[ 0:nwpa ]
 *ieind[0,LXN*LYN-1],ind=ix*LXN
 *nsp [ 0:nwpa ],flags,0:land,1:water,2:open boundary
 * */
void f_initGPar(threadGroup *gi){
  //init in ouput_cal
  //double *pvkd,*pebdep;
  //float *ape , *tpf , *aet , *h1_3, *bv ;
}
Part *spart;
typedef struct _FC_PARA{
  int _kl,_jnthet,_kld,_NBVDEP,_SIGMALVL,_LOGSCURR;
  int nwps,nwpc,nwpa;
  int LXB,LYB,LXN,LYN;
  int mpi_id,mpi_comm_wav,mpi_npe;
  float constwindx,constwindy;
}FC_PARA;
struct Iepos{
  int ix,iy;
};
void  initmtpar(){
  md.pis=gppar.pis;
  md.ipos12=gppar.ipos12;//nwpc
  md.ipos8=gppar.ipos8;//nwpc
}
void InitThreadsGPar(FC_PARA *fcp,struct Iepos*iepos,int*ieind,int*nsp){
  int acwv;
  //printf("AAA %d %d\n",fcp->nwpc,fcp->nwpa);
  md.mpi_id=fcp->mpi_id;
  md.mpi_comm_wav=fcp->mpi_comm_wav;
  md.mpi_npe=fcp->mpi_npe;
  md.pisp.mpi_id=fcp->mpi_id;
  md.nwps  =fcp->nwps  ;
  md.nwpc  =fcp->nwpc  ;
  md.nwpa  =fcp->nwpa  ;
  md.nwpw  =0;
  md.LXB   =fcp->LXB  ;
  md.LYB   =fcp->LYB  ;
  md.LXN   =fcp->LXN  ;
  md.LYN   =fcp->LYN  ;
  md.nsp   =nsp;
  md.ieind =ieind;
  md.iepos =iepos;
  md.constwind.wx=fcp->constwindx;
  md.constwind.wy=fcp->constwindy;
 	md.constwind.wv=md.constwind.wx*md.constwind.wx+md.constwind.wy*md.constwind.wy;
  md.constwind.wi=0;
 	acwv=md.constwind.wv*100;
  if(acwv>=10){
  	md.uconstwind=2;
  }else if(acwv>=1){
  	md.uconstwind=1;  	
  }else {
  	md.uconstwind=0;  	
  }
  initmd();
  NEWZN(spart,4*4+1);
  GPart(spart,ieind,4,4,md.nwps,md.nwpc,md.nwpa,md.LXN,md.LYN,-(int)(1+100*1/38));
  /* partion here,alloc mems
   * ZZZ fill  group info
   * */
  md.s_segs=spart[0].send_segs;
  md.g2l=NEWPN(md.g2l,md.nwpa+1);
#if 0
  char fn[256];
  sprintf("TTA_%2.2d.txt",mpi_id);
  FILE*fo=fopen(fn,"wt");
  CPYN(md.g2l,spart[0].grp2glb,md.nwpa+1);
  for(int ii=0;ii<md.nwpa;ii++){
    fprintf(fo,"%d %d %d %d\n",ii,md.g2l[ii].ipart,md.g2l[ii].iac);
  }
  fclose(fo);
#endif
  for(int p=0;p<4*4;p++){
    threadGroup *gi=md.grps[p];
    bindcpu(p*NCorePClu );
    Part *pt=spart+p+1;
    //printf("GGGGGGEEE %d %d %5d \n",mpi_id,p,pt->nwpa);
    gi->block_id=pt->block_id;
    gi->NPART   =pt->NPART;
    gi->LXB =pt->glb_ixb;
    gi->LYB =pt->glb_jyb;
    gi->LXN     =pt->GRXN;
    gi->LYN     =pt->GRYN;
    gi->recti   =pt->recti;
    gi->recto   =pt->recto;
    gi->nwpt  =pt->nwpt;//capital
    gi->nwps  =pt->nwps;//neighbor
    gi->nwpc  =pt->nwpc;//inside
    gi->nwpo  =pt->nwpo;//outside
    gi->nwpa  =pt->nwpa;//Z
    gi->s_segs=pt->send_segs;
    size_t size=0;
    char*buf;

#define VAVARS \
    VVAR(wav_ec,double,(gi->nwpa+1)*mkj);\
    VVAR(wav_et,double,(gi->nwpa+1)*mkj);\
    VVAR(l2gind,int   ,(gi->nwpa+1)    );\
    VVAR(nsp   ,char  ,(gi->nwpc+1)    );

#define  VVAR(vn,T,N) size+=sizeof(T)*(N)+256
  VAVARS
#undef VVAR
  buf=gi->fbuff=NEWN(char,size);
#define  VVAR(vn,T,N) gi->vn=(T*)buf;buf=(char*)( (((long)buf)+sizeof(T)*(N)+255) & (-256) )
  VAVARS
#undef VVAR
#undef VAVARS


    CPYN(gi->l2gind,pt->grp2glb,gi->nwpa+1);
    for(int i=1;i<=gi->nwpc;i++){
      int l2g=gi->l2gind[i];
      gi->nsp[i]=md.nsp[l2g];
    }
    /* not used:
     * recv_segs
     * pemask
     * iepos
     * grp_ieind
     * */
  }
  bindcpu(0 );
  for(int p=0;p<4*4;p++){
    gp_freepart(spart+p);
  }
  free(spart);
}
extern int mpi_id;
void c_fillsendcomm_(){
  sendcomm*scs[16+1];
  sendcomm *sc;
  scs[0]=sc=&md.s_segs;
  md.wav_et=gppar.wav_et;
  md.wav_ec=gppar.wav_ec;
  sc->src_ee=md.wav_et;
  //printf("FILL %d  %p %p\n",mpi_id,sc->src_ee,md.wav_et); fflush(stdout);
  for(int p=0;p<4*4;p++){
    scs[p+1]=sc=&md.grps[p]->s_segs;
    sc->src_ee=md.grps[p]->wav_et;
  }
  gp_filldst(scs,16);
}
void C_InitImplsch(ImplschPar *pisp_,ImplschP*pip,float*pwvs){
	int isize;
  md.pisp.pip=pip;
  md.pwvs=(windvs*)pwvs;//nwpc
	isize=(long)&(((ImplschPar *)0)->sds);
  memcpy((void*)&md.pisp,(void*)pisp_,isize );//for align so copy it
  //md.pisp.pis=md.pis;
  //md.pisp.ipos8=md.ipos8; // nwpc
	md.pisp.Dm2=-2;md.pisp.D1=1.;md.pisp.enh=1.;
  md.pisp.zpi=3.1415926535897932384626433832795*2.;
	md.pisp.D1m30=1e-30;
	md.pisp.F44_24778=44.24778;
	md.pisp.F16_0=16.;
	md.pisp.F0_467=0.467;//mean3			
	md.pisp.F18_20832=18.20832;
	md.pisp.F0_001930368=0.001930368;
	md.pisp.F0_5=0.5;
  md.pisp.F1_0       = 1        ;
	md.pisp.F0_75      = 0.75     ;
	md.pisp.F5_5       = 5.5      ;
	md.pisp.F0_833     = 0.833    ;
	md.pisp.Fm1_25     =-1.25     ;
	md.pisp.F28_       = 28       ;
	md.pisp.F0_80      = 0.80     ;
	md.pisp.F0_065     = 0.065    ;
	md.pisp.F0_001     = 0.001    ;
  md.pisp.tztp       = 1.2      ;
  md.pisp.tzzpi      = md.pisp.zpi*1.099314;//zpi*1.099314;
}
void C_Setwind(){ // ZZZ should opt to sdma 
	if(md.uconstwind==0){
	  ImplschP*pip=md.pisp.pip;
		windvs * wxy=md.pwvs;
    for(int ig=1;ig<md.ngrp;ig++){
      threadGroup *pg=md.grps[ig];
      ImplschP*pip=pg->pisp.pip;
      for(int i=1;i<=pg->nwpc;i++){
        pip[i].wvs=wxy[pg->l2gind[i]];
      }
    }
	}else{
		static int first=1;
		if(first){
      for(int ig=1;ig<md.ngrp;ig++){
        threadGroup *pg=md.grps[ig];
        if(md.uconstwind==2){
          ImplschP*pip=pg->pisp.pip;
          int     nwpc=pg->nwpc;
          for(int i=1;i<=pg->nwpc;i++){
            pip[i].wvs  =md.constwind;
          }
        }else{
          ImplschP*pip=pg->pisp.pip;
          int     nwpc=pg->nwpc;
          windvs * wxy=md.pwvs;
          for(int i=1;i<=pg->nwpc;i++){
            pip[i].wvs=wxy[pg->l2gind[i]];
          }				
        }
        first=0;			
      }
		}
	}
}

int CheckStop(){
  return md.stop_now;
}
void caliacbe(HTHREADINFO ti){
  int nt,n1,n2,nm,nn1,nwps,nwpc;
  float sc,ntsc;
  int id=ti->ind-ti->pg->NSpecial;
  nt=ti->pg->Nthreads-ti->pg->NSpecial;
  sc=0.3;ntsc=nt-sc;
  nwpc=ti->pg->nwpc;
  nwps=ti->pg->nwps;
  n1=(nwpc-nwps)/ntsc+1;
  ti->iacb1=nwpc-(nt-id  )*n1+1  ;
  ti->iace1=nwpc-(nt-id-1)*n1  ;
  if(ti->iacb1<=nwps)ti->iacb1=nwps+1;
  if(ti->iace1<ti->iacb1)ti->iace1=ti->iacb1-1;
  nwps=nwpc-n1*nt;
  if(nwps<ti->pg->nwps) nwps=ti->pg->nwps;
  n2=nwps/ntsc+1;
  ti->iacb2=nwps-(nt-id  )*n2+1  ;
  ti->iace2=nwps-(nt-id-1)*n2  ;
  if(ti->iacb2<=0)ti->iacb2=1;
  if(ti->iace2<ti->iacb2)ti->iace2=ti->iacb2-1;
  fprintf(ti->fo,"VGT %2.2d %2.2d:%4d %4d %4d %4d:%4d %4d: %4d %4d %4d %4d %4d %4d %f\n",
         ti->igrp,ti->ind,
         ti->iacb1,ti->iace1,ti->iacb2,ti->iace2,
         ti->iace1-ti->iacb1+1,ti->iace2-ti->iacb2+1,
         nwps,nwpc,n1,n2,nt,id,ntsc
         );
}
//#define DTIMES
#ifdef DTIMES
#define BT(i)  tscb(ti,i)
#define ET(i)  tsce(ti,i)
#define EBT(i) tsceb(ti,i);
#define PH()   prth(ti);
#define PT(tag)   prtsc(ti,tag);
#else
#define BT(i)
#define ET(i)
#define EBT(i)
#define PH()
#define PT(tag)   
#endif
void c_sendboundary_(int *RFG){
  //printf("main set FLG %d\n",*RFG);
  gp_send_data(&md.s_segs,sizeof(double)*mkj,16,*RFG);
  //gp_sdma_send_data(&md.s_segs,sizeof(double)*mkj,16,*RFG);
}
void c_checkboundary_(int *RFG){
  //printf("main check FLG %d\n",*RFG);
  gp_check_recv(&md.s_segs,16,*RFG);
}
/* M:main mpi thread
 * G:group main thread
 * S:group sub thread 
 * W:wait
 * S:set 
 * */
#define SWG  sWaitGrp    (ti,RFB,__LINE__) /*wait grp ready */
#define SSG  sSetGrp     (ti,RFB) /*set sub  ready */

#define GWS  gWaitSubs   (ti,RFB) /*wait sub thread*/
#define GSS  gSetSubs    (ti,RFB) /*set gmthread ok*/

#define GWM  gWaitMain   (ti,RFB) /*wait mmthread  */
#define GSM  gSetMain    (ti,RFB) /*set grp ok     */

#define MWG  mWaitGrps   (   RFB) /*mmt */
#define MSG  mSetGrps    (   RFB) /*mmt */

#define SWGR sWaitGrpr   (ti,RFB) /*wait grp ready */
#define GWSR gWaitSubsr  (ti,RFB) /*wait sub thread*/
#define GWMR gWaitMainr  (ti,RFB) /*wait mmthread  */

#define DBG  //printf("AA %d %d %2.2d %2.2d %d %d %d\n",__LINE__,mpi_id,ti->nstep,ti->igrp,ti->ind,RFB,RFG);
#define DBGA //printf("AA %d %d %2.2d %2.2d %d %d %d\n",__LINE__,mpi_id,ti->nstep,ti->igrp,ti->ind,RFB,RFG);
void prth(HTHREADINFO ti){
  char *hdg= "GT(s):id ig ind|   prop 1|      GWM|implsch 1|setspec 1|  accum 1|  ckrecv1|propagat2|implsch 2| setspec2|    GWS  | senddata|  accum 2|  GWS    |";       
  //         "SF(s): 0  4  0:|  0.01900| 13.73147|  0.19263|  0.00011|  0.00016| 10.82476|  0.00109|  0.00006|  0.00008| 15.51763|  0.00181|  0.00050|  0.00162|  0.00000|  0.00000|  0.00000 
  char *hds= "ST(s):id ig ind|   prop 1|      SWG|implsch 1|setspec 1|  accum 1|      SWG|propagat2|implsch 2| setspec2|    SSG  |      NUL|  accum 2|  SWG    |";       
  char *hdm= "MT(s):id ig ind|   RWIND |  checkb |Exchange |setbbound|      NUL|      MWG|SetWind  |implsch 2| setspec2|    SSG  |      NUL|  accum 2|  SWG    |";       
  fprintf(ti->fo,"%s\n",hds);
  fprintf(ti->fo,"%s\n",hdm);
  fprintf(ti->fo,"%s\n",hdg);
}
void _threadMain_(HTHREADINFO ti){
  int indg=ti->indg;
  int RFB=1,RFG=1;
  threadGroup*gi=ti->pg;
  double tb,td;
  double tbs[40],tds[40];
  bindcpu(ti->indg);
  //sleep(1000); exit(0);
  memset(tds,0,sizeof(tds));
  //ti->mstate=&gi->state;
  caliacbe(ti);	  //not need wait	, ready here
  RFB=1;
  ti->nstep=-1;
  if(ti->MainThread){
    ti->pg->nstep=-1;
    DBGA;
    GWM; 
    DBGA;
    g_initPropgats(gi); 
    DBGA;
    g_initImplsch(gi);
    DBGA;
    g_initSetspec(gi);
    DBGA;
    g_initOutput_cal(gi);
    DBGA;
    GSS; GSM;
    DBGA;
  }else{ 
    SWG; 
  }
  setspec   (gi,gi->wav_et,1,gi->nwpc,1)  ;
  ti->wav_et=gi->wav_et;
  ti->wav_ec=gi->wav_ec;
  RFB=10;
  if(ti->MainThread){ // forced sync
    DBG;
    GWS; DBG;
    GSM; DBG;
    GWM; DBG;
    GSS; DBG;
  }else{
    SSG; SWG;
  }
  if(ti->MainThread){
    RFG=1; DBG;
    gp_send_data(&gi->s_segs,sizeof(double)*mkj,16,++RFG);
    DBG;
    //gp_check_recv(&gi->s_segs,16,RFG);
    tscinit(ti);
    for(;;){
      ti->pg->nstep++;
      ti->nstep++;
      if(RFB>100000){// reset Sync ID
        RFB=10; GWSR; GSM; GWMR; GSS;
      }
      DBG;
      BT(0);propagats(gi,ti->wav_et,ti->wav_ec,ti->iacb1,ti->iace1 );
      DBG;
      //wait wind ok
      RFB++; GSM;
      DBG;
      EBT(1);GWM; GSS; //1
      DBG;
      EBT(2);implschs (gi,ti->wav_ec,ti->wav_et,ti->iacb1,ti->iace1 );
      EBT(3);setspec  (gi,ti->wav_et,ti->iacb1,ti->iace1,2)  ;
      EBT(4);accumea_ (&indg,ti->wav_et,&ti->iacb1,&ti->iace1);
      DBG;
      //wait boudary ok
      EBT(5);gp_check_recv(&gi->s_segs,16,RFG); //fixme:RFG
      DBG;
      RFB++; //GWM; wait by gp_check_recv
      GSS;GSM;//2
      DBGA;
      EBT(6);propagats(gi,ti->wav_et,ti->wav_ec,ti->iacb2,ti->iace2 );
      EBT(7);implschs (gi,ti->wav_ec,ti->wav_et,ti->iacb2,ti->iace2 );
      EBT(8);setspec  (gi,ti->wav_et,ti->iacb2,ti->iace2,2);
      DBG;
      RFB++;GSM; 
      DBG;
      EBT(9); GWS;//Send_data;//3
      DBG;
      EBT(10);gp_send_data(&gi->s_segs,sizeof(double)*mkj,16,++RFG);//fixme:RFG
      DBG;
      //gp_sdma_send_data(&gi->s_segs,sizeof(double)*mkj,16,++RFG);
      EBT(11);accumea_  (&indg,ti->wav_et,&ti->iacb2,&ti->iace2);
      DBG;
      EBT(12);RFB++;GSS;GWS;//4 for debug
      ET(12);
      DBG;
      if(md.hist_eot){
        //checkoutput_(ti,&indg,&RFB,ti->wav_et,&ti->iacb1,&ti->iace1,&ti->iacb2,&ti->iace2,&md.hist_eot);
      }
      if(ti->nstep==0){
        PT("gt");
      }
      if(CheckStop())break;
    }
    PH(); PT("GT");
    RFB+=1000;
    //GWS;
    GSM;
    //GWM; //if NEED Clean
    GSS;
  }else{
    tscinit(ti);
    for(;;){
      ti->nstep++;
      if(RFB>100000){// reset Sync ID
        RFB=10; SSG; SWGR;
      } 
      BT(0);propagats(gi,ti->wav_et,ti->wav_ec,ti->iacb1,ti->iace1 );
      //wait wind ok
      EBT(1);RFB++; SSG; SWG; //1
      EBT(2);implschs  (gi,ti->wav_ec,ti->wav_et,ti->iacb1,ti->iace1 );
      EBT(3);setspec   (gi,ti->wav_et,ti->iacb1,ti->iace1,2)  ;
      EBT(4);accumea_  (&indg,ti->wav_et,&ti->iacb1,&ti->iace1);
      //wait boudary ok
      EBT(5);RFB++;SWG; //2
      EBT(6);propagats(gi,ti->wav_et,ti->wav_ec,ti->iacb2,ti->iace2 );
      EBT(7);implschs  (gi,ti->wav_ec,ti->wav_et,ti->iacb2,ti->iace2 );
      EBT(8);setspec   (gi,ti->wav_et,ti->iacb2,ti->iace2,2);
      EBT(9);RFB++;SSG;  // 3
      EBT(10);
      EBT(11);accumea_  (&indg,ti->wav_et,&ti->iacb2,&ti->iace2);
      EBT(12);RFB++;SSG;SWG; //4
      ET(12);
      if(md.hist_eot){
        //checkoutput_(ti,&indg,&RFB,ti->wav_et,&ti->iacb1,&ti->iace1,&ti->iacb2,&ti->iace2,&md.hist_eot);
      }
      if(ti->nstep==0){
        PT("st");
      }
      if(CheckStop())break;
    }
    PT("ST");
    RFB+=1000;
    SSG;//not need wait
    //SWG;
  }
}
