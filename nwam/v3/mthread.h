#include <pthread.h>
#define MCOREPC 38
#define MCLUST  16
typedef struct _THREADINFO THREADINFO,*HTHREADINFO;
typedef void (*TFunc)(HTHREADINFO);
typedef void *(*TSFunc)(void*);
typedef struct _THREADINFO{
	short ind,Nthreads,MainThread,igrp;	   //4*2 2
	short nlocv,detach,indg  ;//3*2 2
  int volatile state;       //1
  int volatile *mstate;
  struct _THREADINFO  *threads;//2
  TFunc tfun; //8 2
	void*para;  //8 2
	void**locv; //8 2
	pthread_t threadid;//8 2
  FILE*fo;
  TUSERDECLARE ;
}THREADINFO;
#define MAXTHREADS 40
typedef struct _threadData{
  short Nthreads,NSpecial,PThreadInited,idMainThread ;
  THREADINFO *threads;
  int mpi_id,mpi_npe;
  int volatile state;
  NUSERDECLARE;
}threadData ;
__BEGIN_DECLS
extern threadData md;
void _threadMain_(HTHREADINFO ti);
//tread ctrl 
void threadStart(int bi,int ei);
void threadEnd(int bi,int ei);
int bindcpu(int id);
void opentf(HTHREADINFO ti);

void sSetState(HTHREADINFO ti,int state);
void sWaitState(HTHREADINFO ti,int state);
void sWaitStater(HTHREADINFO ti,int state);

void swaitstater_(HTHREADINFO ti,int *state);
void swaitstate_(HTHREADINFO ti,int *state);
void ssetstate_(HTHREADINFO ti,int *state);

void waitstate_(int *state,int *nstep);
void waitstater_(int *state,int *nstep);
void setstate_(int *state,int *nstep);

__END_DECLS
