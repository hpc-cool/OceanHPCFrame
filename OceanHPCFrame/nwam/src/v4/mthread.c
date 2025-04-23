#define _GNU_SOURCE
#include <sched.h>
#include <pthread.h>
#include "svars.h"
int bindcpu(int id){
  cpu_set_t mask;  //CPU核的集合
  CPU_ZERO(&mask);    //置空
  CPU_SET(id,&mask);   //设置亲和力值
  //printf("bindcpu %d\n",id);
  if (sched_setaffinity(0, sizeof(mask), &mask) == -1){//设置线程CPU亲和力
    printf("warning: could not set CPU affinity, continuing...\n");
    return -1;
  }
  return 0;
}
threadData md={0};
int bindcpu(int id);
void threadStart(int bi,int ei);
void threadEnd(int bi,int ei);
void initmd(){
  int i,j;
  if(md.PThreadInited)return ;
  if(NThPClu>MCOREPC){
    printf("NThPClu %d biger than surported MCOREPC %d\n",NThPClu,MCOREPC);
    exit(-1);
  }
  if(NGrpPNode>MCLUST){
    printf("NGrpPNode %d biger than surported MCLUST %d\n",NGrpPNode,MCLUST);
    exit(-1);
  }
  md.ngrp=NGrpPNode;
  for(i=0;i<NGrpPNode;i++){
    bindcpu(i*NCorePClu );
    threadGroup *p=md.grps[i]=HNEW(threadGroup );
    memset(p,0,sizeof(threadGroup));
    p->pmd=&md;
    p->ngrp=NGrpPNode;
    p->igrp=i;
    p->NSpecial=0;
    p->Nthreads=NThPClu;
    for(j=0;j<NThPClu;j++){
      struct _THREADINFO  *t=p->threads+j;
      t->pg=p;
      t->threads=p->threads;
      t->ind=j;
      t->Nthreads=NThPClu;
      t->MainThread=(j==0);
      t->igrp=i;
      t->indg=i*NCorePClu +j;
      //t->mstate=&p->state;
    }
  }
  bindcpu(ManageCoreId  );
  md.PThreadInited=1;
}
void setflag_(int *hist_eot,int *stop_now,int *rest_eot){
  md.hist_eot=*hist_eot;
  md.stop_now=*stop_now;
  md.rest_eot=*rest_eot;
  for(int i=0;i<NGrpPNode ;i++){
    md.grps[i]->hist_eot=md.hist_eot;
    md.grps[i]->stop_now=md.stop_now;
    md.grps[i]->rest_eot=md.rest_eot;
  }
}
void C_ctStart(){
  //printf("Start Threads %d x %d = %d\n",NThPClu,NGrpPNode,NGrpPNode *NThPClu);
  threadStart(0,NGrpPNode *NThPClu);
}
void C_ctEnd(){
  threadEnd(0,NGrpPNode *NThPClu);
}
int GetVInt(int volatile *volatile p);
__attribute__ ((optnone)) void ntdelay(int n){
  usleep(0);
//  static int vt=0; for(int i=0;i<n;i++) vt=(vt*1357+2581);
}
#define MAXCHECK 1000000
#define FWAIT(NM,COP) int Wait_##NM(volatile int *pv,int cv,int nl){ \
  for(int i=0;GetVInt(pv) COP cv;i++){ \
  if(i>=MAXCHECK){ printf("wait state timeout %d %s %d\n",*pv,#COP,cv);exit(0); }\
  ntdelay(10); } return 0;\
}
FWAIT(LNE,==);
FWAIT(LEQ,!=);
//FWAIT(LGE,<);
FWAIT(LLE,>);
int Wait_LGE(volatile int *pv,int cv,int nl){ 
  for(int i=0;GetVInt(pv) < cv;i++){ 
    if(i>=MAXCHECK){ 
      printf("wait state timeout %d:%p %d %s %d\n",nl,pv,GetVInt(pv),"<",cv);
      exit(0); 
    }
    ntdelay(10); 
  } 
  if(GetVInt(pv)<cv){
    printf("EEEEEEEEEE %d %d %d\n",GetVInt(pv),*pv,cv);
  }
  return 0;
}
//#define DBGSYNC
#define THLOG
#ifdef DBGSYNC
#define THLOG
#endif
void opentf(HTHREADINFO ti){
#ifdef THLOG
  //ti->fo=stdout;
  if(!ti->fo){
    char fn[256];
    sprintf(fn,"MTW_%2.2d_%2.2d.log",ti->igrp,ti->ind);
    //printf("Openfile %s\n",fn);
    ti->fo=fopen(fn,"wt");
  }
#endif
}
//called by sub threads
//thread wait gmain state,pg->state
void sWaitGrp(HTHREADINFO ti,int state,int nl){
#ifdef DBGSYNC
  fprintf(ti->fo,"SWG:%3.2d %3.2d %3.2d\n",md.nstep,state,ti->ind);
  fflush(ti->fo);
#endif
  Wait_LGE(&ti->pg->state,state,nl);
}
void sWaitGrpr(HTHREADINFO ti,int state){
#ifdef DBGSYNC
  fprintf(ti->fo,"SRG:%3.2d %3.2d %3.2d\n",md.nstep,state,ti->ind);
  fflush(ti->fo);
#endif
  Wait_LLE(&ti->pg->state,state,__LINE__);
}
//thread set state for gmain ,ti->state
void sSetGrp(HTHREADINFO ti,int state){
#ifdef DBGSYNC
  fprintf(ti->fo,"SSG:%3.2d %3.2d %3.2d %3.2d\n",md.nstep,ti->state,state,ti->ind);
  fflush(ti->fo);
#endif
  //ti->state=state;
  ti->pg->sstate[ti->ind*MSBG]=state;
}
//called by main group threads
//gmain thread wait sub threads,threads[x].state
void gWaitSubs(HTHREADINFO ti,int state){//mt
  int i;
  int ib=ti->pg->NSpecial+1;
#ifdef DBGSYNC
  if(!ti->fo) opentf(ti);
  fprintf(ti->fo,"GWS:%3.2d %3.2d:",md.nstep,state);
  for(i= ib;i<md.Nthreads;i++){
    fprintf(ti->fo,"%3.2d ",ti->threads[i].state);
  }
  fprintf(ti->fo,"\n ");
  fflush(ti->fo);
#endif
  for(i=ib;i<ti->Nthreads;i++){
#ifdef DBGSYNC
    printf("GWS:%3.2d %3.2d %d %d\n",md.nstep,state,i,ti->threads[i].state);
    fflush(stdout);
#endif
    //Wait_LGE(&ti->threads[i].state,state,__LINE__);
    Wait_LGE(&ti->pg->sstate[i*MSBG],state,__LINE__);
  }
}
void gWaitSubsr(HTHREADINFO ti,int state){//mt
  int i;
  int ib=ti->pg->NSpecial+1;
#ifdef DBGSYNC
  if(!ti->fo) opentf(ti);
  fprintf(ti->fo,"GRS:%3.2d %3.2d:",md.nstep,state);
  for(i=ib;i<md.Nthreads;i++){
    fprintf(ti->fo,"%3.2d ",ti->threads[i].state);
  }
  fprintf(ti->fo,"\n ");
  fflush(ti->fo);
#endif
  for(i=ib;i<ti->Nthreads;i++){
    Wait_LLE(&ti->pg->sstate[i*MSBG],state,__LINE__);
  }
}
//gmain thread set for subthread ,pg->state
void gSetSubs (HTHREADINFO ti,int state){
#ifdef DBGSYNC
  fprintf(ti->fo,"GSS:%3.2d %3.2d %3.2d %3.2d\n",md.nstep,ti->state,state,ti->ind);
  fflush(ti->fo);
#endif
  ti->pg->state=state;
}

//gmain thread wait mmthreads,pg->gstate
void gWaitMain(HTHREADINFO ti,int state){//mt
  int i;
#ifdef DBGSYNC
  if(!ti->fo) opentf(ti);
  fprintf(ti->fo,"GWM:%3.2d %3.2d %3.2d:\n",md.nstep,state,ti->igrp); fflush(ti->fo);
  printf("GWM:%3.2d %d %d %3.2d:\n",md.nstep,state,ti->pg->gstate,ti->igrp); fflush(stdout);
#endif
  Wait_LGE(&ti->pg->gstate,state,__LINE__);
#ifdef DBGSYNC
  //fprintf(ti->fo,"GWM:%3.2d %3.2d %d %3.2d:\n",md.nstep,state,ti->pg->gstate,ti->igrp); fflush(ti->fo);
  printf("GWME:%3.2d %3.2d %d %3.2d:\n",md.nstep,state,ti->pg->gstate,ti->igrp); fflush(stdout);
#endif
}
void gWaitMainr(HTHREADINFO ti,int state){//mt
  int i;
#ifdef DBGSYNC
  if(!ti->fo) opentf(ti);
  fprintf(ti->fo,"GRM:%3.2d %3.2d %3.2d:\n",md.nstep,state,ti->igrp);
  printf("GRM:%3.2d %3.2d %3.2d:\n",md.nstep,state,ti->igrp); fflush(stdout);
  fflush(ti->fo);
#endif
  Wait_LLE(&ti->pg->gstate,state,__LINE__);
}
//gmain thread set for mmthread ,md.gstate[(x)*MSB]
void gSetMain(HTHREADINFO ti,int state){
#ifdef DBGSYNC
  fprintf(ti->fo,"GSM:%3.2d %3.2d %3.2d %3.2d\n",md.nstep,md.gstate[(ti->igrp)*MSB],state,ti->igrp);
  fflush(ti->fo);
#endif
  md.gstate[(ti->igrp)*MSB]=state;
#ifdef DBGSYNC
  printf("GSME:%3.2d %3.2d %3.2d %3.2d\n",md.nstep,md.gstate[(ti->igrp)*MSB],state,ti->igrp);
  fflush(stdout);
#endif
}
//called by main main threads
//mmthread wait gmain thread,md.gstate[(x)*MSB]
void mWaitGrps(int state){//mmt
  int i;
#ifdef DBGSYNC
  printf("MWG:%d %d %d:",md.nstep,state,md.ngrp);
  for(i=0;i<md.ngrp;i++){ printf("%d ",md.gstate[(i)*MSB]); }
  printf("\n ");
  fflush(stdout);
#endif
  HTHREADINFO ti=md.grps[0]->threads;
  for(i=0;i<md.ngrp;i++){
    Wait_LGE(&md.gstate[(i)*MSB],state,__LINE__);
    //printf("VE %2d: %p %3.2d %d\n",i,&md.gstate[(i)*MSB],GetVInt(&md.gstate[(i)*MSB]),state); fflush(stdout);
  }
}
void mWaitGrpsr(int state){//mmt
  int i;
#ifdef DBGSYNC
  HTHREADINFO ti=md.grps[0]->threads;
  if(!ti->fo) opentf(ti);
  fprintf(ti->fo,"MRG:%3.2d %3.2d:",md.nstep,state);
  for(i=0;i<md.ngrp;i++){
    fprintf(ti->fo,"%3.2d ",md.gstate[(i)*MSB]);
  }
  fprintf(ti->fo,"\n ");
  fflush(ti->fo);
#endif
  for(i=0;i<md.ngrp;i++){
#ifdef DBGSYNC
    fprintf(ti->fo,"V %d:%3.2d \n",i,md.grps[i]->state); fflush(ti->fo);
#endif
    Wait_LLE(&md.gstate[(i)*MSB],state,__LINE__);
  }
#ifdef DBGSYNC
  fprintf(ti->fo,"\n ");
#endif
}
//mmthread set for gmain threads,md.grps[x].gstate
void mSetGrps(int state){// mmt
#ifdef DBGSYNC
  HTHREADINFO ti=md.grps[0]->threads;
  if(!ti->fo) opentf(ti);
  fprintf(ti->fo,"MSG:%3.2d %3.2d %3.2d\n",md.nstep,md.state,state);
  fflush(ti->fo);
  printf("MSG:%3.2d %3.2d %d\n",md.nstep,state,md.ngrp);
#endif
  md.state=state;
  for(int i=0;i<md.ngrp;i++){
    md.grps[i]->gstate=state;
  }
#ifdef DBGSYNC
  {
    HTHREADINFO ti=md.grps[0]->threads;
    printf("MSG:%3.2d %3.2d\n",md.nstep,state);
    for(int i=0;i<md.ngrp;i++){
      printf("%d:%3.2d ",i,md.grps[i]->gstate); 
    }
    printf("\n ");
    fflush(stdout);
  }
#endif
}
// Fortran interface
void swaitgrp_   (HTHREADINFO ti,int *state){ sWaitGrp   (ti,*state,__LINE__); }
void swaitgrpr_  (HTHREADINFO ti,int *state){ sWaitGrpr  (ti,*state); }
void ssetgrp_    (HTHREADINFO ti,int *state){ sSetGrp    (ti,*state); }
void gwaitsubs_  (HTHREADINFO ti,int *state){ gWaitSubs  (ti,*state); }
void gwaitsubsr_ (HTHREADINFO ti,int *state){ gWaitSubsr (ti,*state); }
void gsetsubs_   (HTHREADINFO ti,int *state){ gSetSubs   (ti,*state); }
void gwaitmain_  (HTHREADINFO ti,int *state){ gWaitMain  (ti,*state); }
void gwaitmainr_ (HTHREADINFO ti,int *state){ gWaitMainr (ti,*state); }
void gsetmain_   (HTHREADINFO ti,int *state){ gSetMain   (ti,*state); }
void waitstate_  (int *state    ,int *nstep){ md.nstep=*nstep; mWaitGrps (*state); }
void waitstater_ (int *state    ,int *nstep){ md.nstep=*nstep; mWaitGrpsr(*state); }
void setstate_   (int *state    ,int *nstep){ md.nstep=*nstep; mSetGrps  (*state); }

static inline uint64_t atsc(void) {
  uint64_t tsc;
  asm volatile("mrs %0, cntvct_el0" : "=r" (tsc));    //读取系统时间戳
  return tsc;
}
void tscinit(HTHREADINFO ti){
  memset(ti->tscds,0,sizeof(ti->tscds));
}
void tscb(HTHREADINFO ti,int id){
  ti->tscbs[id]=atsc();
}
void tsce(HTHREADINFO ti,int id){
  uint64_t d=atsc();
  ti->tscds[id]+=((d-ti->tscbs[id])&0xFFFFFFFFFFFFFFL);//most tsc is 56bits
}
void tsceb(HTHREADINFO ti,int id){
  uint64_t d=atsc();
  ti->tscds[id-1]+=((d-ti->tscbs[id-1])&0xFFFFFFFFFFFFFFL);//most tsc is 56bits
  ti->tscbs[id]=d;
}
void tscinit_(){
  HTHREADINFO ti=&md.tm;
  memset(ti->tscds,0,sizeof(ti->tscds));
}
void tscb_(int *id_){
  HTHREADINFO ti=&md.tm;int id=*id_;
  ti->tscbs[id]=atsc();
}
void tsce_(int *id_){
  HTHREADINFO ti=&md.tm;int id=*id_;
  uint64_t d=atsc();
  ti->tscds[id]+=((d-ti->tscbs[id])&0xFFFFFFFFFFFFFFL);//most tsc is 56bits
}
void tsceb_(int *id_){
  HTHREADINFO ti=&md.tm;int id=*id_;
  uint64_t d=atsc();
  ti->tscds[id-1]+=((d-ti->tscbs[id-1])&0xFFFFFFFFFFFFFFL);//most tsc is 56bits
  ti->tscbs[id]=d;
}
void prtsc(HTHREADINFO ti,const char*tag){
  char buf[1024]={0};
  for(int i=0;i<16;i++){
    sprintf(buf+i*10,"|%9.5f     Z",ti->tscds[i]*(1./(100*1000*1000)));
  }
  fprintf(ti->fo,"%2s(s):%2.2d %2.2d %2d:%s\n",tag,mpi_id,ti->igrp,ti->ind,buf);
  fflush(ti->fo);
}
void prtsc_(){
  prtsc(&md.tm,"MT");
}
static void ClearThread(HTHREADINFO ti){
    if(ti->threadid){
       void*vres;
       pthread_cancel(ti->threadid);
       pthread_join(ti->threadid,&vres);
       //ti->state=0;
       ti->pg->sstate[ti->ind*MSBG]=0;
    }
#if 0
    if(ti->nlocv){
      int i;
      for(i=0;i<ti->nlocv;i++){
        free(ti->locv[i]);ti->locv[i]=NULL;
      }
      free(ti->locv);ti->locv=NULL;
    }
    ti->nlocv=0;
#endif
    ti->threadid=0;
}
static void InitThread(HTHREADINFO ti,TFunc tfun,void*para,int ind,int detach,int nlocv,int Mainthread){
  ClearThread(ti);
  ti->tfun=tfun;
  //ti->para=para;
  ti->ind=ind%NThPClu;
  ti->igrp=ind/NThPClu;
  ti->indg=ti->igrp*NCorePClu+ti->ind ;
  ti->detach=detach;
  if(tfun==NULL)return;
  ti->Nthreads=ti->pg->Nthreads;
  ti->nlocv=nlocv;
  //ti->flog=NULL;
  if(Mainthread){
    if(ind!=ti->pg->idMainThread){
      if(ti->pg->idMainThread<ti->Nthreads)ti->threads[ti->pg->idMainThread].MainThread=0;
    }
    ti->pg->idMainThread=ind;
    ti->MainThread=1;
  }else ti->MainThread=0;
}
void zStartThreads(int bthread, int ethreads,TFunc tfun,void*para,int detach,int clear){
  int i,ib,ie,init0=0;
  initmd();
  ib=bthread; ie=ethreads;
  if(ib<=0){
    init0=1;//nthreads=nthreads+1;
    if(ManageCoreId ==0) ib=1;
  }
  //printf("threads :%d %d\n",ib,ie);
  if(bthread<=0){
    bindcpu(ManageCoreId  );
    //printf("MT Init\n");
    InitThread(&md.tm,NULL,NULL,0,0,0,1);
    md.tm.threadid=0;
    md.tm.igrp=-1;
    opentf(&md.tm);
    if(ib>0) md.grps[0]->NSpecial=1;
  }
  if(tfun){
    for(i=ib;i<ie;i++){
      int ig=i/NThPClu;
      int it=i%NThPClu;
      HTHREADINFO tis=md.grps[ig]->threads;
      int mt=0;
      if(it==0)mt=1;
      else if(ig==0&&ManageCoreId  ==0&&it==1)mt=1;
      InitThread(tis+it,tfun,para,i,detach,0,mt);
      opentf(tis+it);
      pthread_create(&tis[it].threadid,NULL,(TSFunc)tfun,&tis[it]);
      if(detach)pthread_detach(tis[it].threadid);
    }
    //printf("SSS End %d/%d\n",i,ie);
  }
}
void threadStart(int bi,int ei){
  zStartThreads(bi,ei,_threadMain_,0,1,1);
}
void threadEnd(int bi,int ei){
  zStartThreads(bi,ei,NULL,0,1,1);
}
