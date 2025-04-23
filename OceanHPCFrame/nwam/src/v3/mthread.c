#define _GNU_SOURCE
#include <sched.h>
#include <pthread.h>
#include "svars.h"
int bindcpu(int id){
  cpu_set_t mask;  //CPU核的集合
  CPU_ZERO(&mask);    //置空
  CPU_SET(id,&mask);   //设置亲和力值
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
  if(md.PThreadInited)return ;
  if(NThPClu>MCOREPC){
    printf("NCorePClu %d biger than surported MCOREPC %d\n",NThPClu,MCOREPC);
    exit(-1);
  }
  if(NProcPNode>MCLUST){
    printf("NProcPNode %d biger than surported MCLUST %d\n",NProcPNode ,MCLUST);
    exit(-1);
  }
  md.Nthreads=0;
  md.PThreadInited=1;
  md.threads=NULL;
}
void setflag_(int *hist_eot,int *stop_now,int *rest_eot){
  md.hist_eot=*hist_eot;
  md.stop_now=*stop_now;
  md.rest_eot=*rest_eot;
}
void exchage_ee(){
  double *pt;
  pt=gppar.wav_ec  ;gppar.wav_et  =gppar.wav_ec  ;gppar.wav_ec  =pt;
}
void C_ctStart(){
  threadStart(0,NThPClu);
}
void C_ctEnd(){
  threadEnd(0,NThPClu);
}
int GetVInt(int volatile *p);
void ntdelay(){
  usleep(1);
  //static int vt=0; for(int i=0;i<10;i++) vt=(vt*1357+2581);
}
#define FWAIT(N,O) void Wait_##N(volatile int *pv,int cv){ while(GetVInt(pv) O cv){ ntdelay(); } }
FWAIT(LNE,==);
FWAIT(LEQ,!=);
FWAIT(LGE,<);
FWAIT(LLE,>);
//#define DBGSYNC
void opentf(HTHREADINFO ti){
#ifdef DBGSYNC
  if(!ti->fo){
    char fn[256];
    sprintf(fn,"MTW_%d.log",ti->ind);
    ti->fo=fopen(fn,"wt");
  }
#endif
}
//called by sub threads
void sWaitState(HTHREADINFO ti,int state){
#ifdef DBGSYNC
  fprintf(ti->fo,"TWT:%3.2d %3.2d %3.2d\n",md.nstep,state,ti->ind);
  fflush(ti->fo);
#endif
  Wait_LGE(ti->mstate,state);
}
void sWaitStater(HTHREADINFO ti,int state){
#ifdef DBGSYNC
  fprintf(ti->fo,"TRT:%3.2d %3.2d %3.2d\n",md.nstep,state,ti->ind);
  fflush(ti->fo);
#endif
  Wait_LLE(ti->mstate,state);
}
void sSetState(HTHREADINFO ti,int state){
#ifdef DBGSYNC
  fprintf(ti->fo,"TST:%3.2d %3.2d %3.2d %3.2d\n",md.nstep,ti->state,state,ti->ind);
  fflush(ti->fo);
#endif
  ti->state=state;
}
//called by main thread
void WaitSubsr(int state){
  int i;
#ifdef DBGSYNC
  HTHREADINFO ti=md.threads;
  if(!ti->fo) opentf(ti);
  fprintf(ti->fo,"MWT:%3.2d %3.2d:",md.nstep,state);
  for(i=1;i<md.Nthreads;i++){
    fprintf(ti->fo,"%3.2d ",md.threads[i].state);
  }
  fprintf(ti->fo,"\n ");
  fflush(ti->fo);
#endif
  for(i=1;i<md.Nthreads;i++){
    Wait_LLE(&md.threads[i].state,state);
  }
}
void mWaitSubs(int state){
  int i;
#ifdef DBGSYNC
  HTHREADINFO ti=md.threads;
  if(!ti->fo) opentf(ti);
  fprintf(ti->fo,"MWT:%3.2d %3.2d:",md.nstep,state);
  for(i=1;i<md.Nthreads;i++){
    fprintf(ti->fo,"%3.2d ",md.threads[i].state);
  }
  fprintf(ti->fo,"\n ");
  fflush(ti->fo);
#endif
  for(i=1;i<md.Nthreads;i++){
#ifdef DBGSYNC
    fprintf(ti->fo,"V %d:%3.2d \n",i,md.threads[i].state); fflush(ti->fo);
#endif
    Wait_LGE(&md.threads[i].state,state);
  }
}
void mWaitSubsr(int state){
  int i;
#ifdef DBGSYNC
  HTHREADINFO ti=md.threads;
  if(!ti->fo) opentf(ti);
  fprintf(ti->fo,"MRT:%3.2d %3.2d:",md.nstep,state);
  for(i=1;i<md.Nthreads;i++){
    fprintf(ti->fo,"%3.2d ",md.threads[i].state);
  }
  fprintf(ti->fo,"\n ");
  fflush(ti->fo);
#endif
  for(i=1;i<md.Nthreads;i++){
#ifdef DBGSYNC
    fprintf(ti->fo,"V %d:%3.2d \n",i,md.threads[i].state); fflush(ti->fo);
#endif
    Wait_LLE(&md.threads[i].state,state);
  }
}

void mSetSubs(int state){
#ifdef DBGSYNC
  HTHREADINFO ti=md.threads;
  if(!ti->fo) opentf(ti);
  fprintf(ti->fo,"MST:%3.2d %3.2d %3.2d\n",md.nstep,md.state,*state);
  fflush(ti->fo);
#endif
  md.state=state;
}
void swaitstate_ (HTHREADINFO ti,int *state){ sWaitState (ti,*state); }
void swaitstater_(HTHREADINFO ti,int *state){ sWaitStater(ti,*state); }
void ssetstate_  (HTHREADINFO ti,int *state){ sSetState(ti,*state);   }

void waitstate_  (int *state,int *nstep)    { mWaitSubs (*state); md.nstep=*nstep; }
void waitstater_ (int *state,int *nstep)    { mWaitSubsr(*state); md.nstep=*nstep; }
void setstate_   (int *state,int *nstep)    { mSetSubs  (*state); md.nstep=*nstep; }

double dclktime();
double dclktimeinit();
static void InitPThread(){
  if(!md.PThreadInited){
    int i;
    dclktimeinit();
    initmd();
  }
}
static void ClearThread(HTHREADINFO ti){
    if(ti->threadid){
       void*vres;
       pthread_cancel(ti->threadid);
       pthread_join(ti->threadid,&vres);
       ti->state=0;
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
  ti->para=para;
  ti->ind=ind;
  ti->detach=detach;
  ti->Nthreads=md.Nthreads;
  ti->nlocv=nlocv;
  //ti->flog=NULL;
  if(Mainthread){
    if(ind!=md.idMainThread){
      if(md.idMainThread<md.Nthreads)md.threads[md.idMainThread].MainThread=0;
    }
    md.idMainThread=ind;
    ti->MainThread=1;
  }
}
void zStartThreads(int bthread, int ethreads,TFunc tfun,void*para,int detach,int clear){
  int i,ib,ie,init0=0;
  ib=bthread;ie=ethreads;
  HTHREADINFO tis;
  InitPThread();
  if(ib<=0){
    init0=1;
    ib=1;
  }
  tis=md.threads;
  md.NSpecial=1;
  //printf("thread :%3.2d %3.2d\n",ib,ie);
  if(tis&&clear){
    for(i=1;i<md.Nthreads;i++){
      ClearThread(tis+i);
    }
    memset(tis,0,(md.Nthreads-ib)*sizeof(THREADINFO));
  }else if(ib<md.Nthreads){
    for(i=ib;i<md.Nthreads;i++){
      ClearThread(tis+i);
    }
    memset(tis+ib,0,(md.Nthreads-ib)*sizeof(THREADINFO));
  }
  if(md.Nthreads<ie){
    md.threads=tis=(THREADINFO*)realloc(tis,ie*sizeof(THREADINFO));
    memset(tis+md.Nthreads,0,(ie-md.Nthreads)*sizeof(THREADINFO));
    md.Nthreads=ie;
  }
  for(int i=0;i<md.Nthreads;i++){
    md.threads[i].Nthreads=md.Nthreads;
  }
  md.threads=tis;
  if(bthread<=0){
    bindcpu((md.mpi_id%NProcPNode  )*NCorePClu);
    memset(tis,0,sizeof(THREADINFO));
    InitThread(tis+0,NULL,NULL,0,0,0,1);
    tis[0].threadid=0;
    opentf(tis);
  }
  if(tfun){
    for(i=ib;i<ie;i++){
      InitThread(tis+i,tfun,para,i,detach,0,0);
      opentf(tis+i);
      pthread_create(&tis[i].threadid,NULL,(TSFunc)tfun,&tis[i]);
      if(detach)pthread_detach(tis[i].threadid);
    }
  }
}
void threadStart(int bi,int ei){
  zStartThreads(bi,ei,_threadMain_,0,1,1);
}
void threadEnd(int bi,int ei){
  zStartThreads(bi,ei,NULL,0,1,1);
}
