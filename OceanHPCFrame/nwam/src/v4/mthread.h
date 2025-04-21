typedef struct _THREADINFO THREADINFO,*HTHREADINFO;
typedef void (*TFunc)(HTHREADINFO);
typedef void *(*TSFunc)(void*);

typedef struct _THREADINFO{
	short ind,Nthreads,MainThread,igrp;	   //4*2 2
	short nlocv,detach,indg  ;//3*2 2
  //int volatile state;       //1    tset,gwait
  uint64_t tscbs[16];
  uint64_t tscds[16];
  struct _threadGroup *pg;  //2
  struct _THREADINFO  *threads;//2
  TFunc tfun; //8 2
	void**locv; //8 2
	pthread_t threadid;//8 2
  FILE*fo;
  TUSERDECLARE ;
}THREADINFO;

#define MCOREPC 38
#define MCLUST  16
#define MSB 16
#define MSBG 16
typedef struct _threadGroup{
  short igrp,ngrp,mainGroup;
  short PThreadInited,idMainThread,Nthreads,NSpecial;
  int volatile state;  //twait,gset
  int volatile gstate;        //msetg gwaitm
  int volatile mstate;        //msetg gwaitm
  struct _threadData *pmd;
  int sstate[(MCOREPC)*MSBG];     //gsetm ,mwaitg
  THREADINFO threads[MCOREPC];
  GUSERDECLARE ;
}threadGroup;
typedef struct _threadData{
  short Nthreads,NSpecial,PThreadInited,idMainThread ,ngrp;
  threadGroup *grps[MCLUST];
  THREADINFO tm;
  int gstate[(MCLUST)*MSB];     //gsetm ,mwaitg
  int volatile state;     //msetg
  NUSERDECLARE;
}threadData ;
__BEGIN_DECLS
extern threadData md;
//tread main  function
void _threadMain_(HTHREADINFO ti);

//tread ctrl 
void threadStart(int bi,int ei);
void threadEnd(int bi,int ei);

//sync function
void sWaitGrp(HTHREADINFO ti,int state,int nl);
void sWaitGrpr(HTHREADINFO ti,int state);
void sSetGrp(HTHREADINFO ti,int state);

void gWaitSubs(HTHREADINFO ti,int state);
void gWaitSubsr(HTHREADINFO ti,int state);
void gSetSubs (HTHREADINFO ti,int state);

void gWaitMain(HTHREADINFO ti,int state);
void gWaitMainr(HTHREADINFO ti,int state);
void gSetMain(HTHREADINFO ti,int state);

void mWaitGrps(int state);
void mWaitGrpsr(int state);
void mSetGrps(int state);
int bindcpu(int id);
void opentf(HTHREADINFO ti);

// Fortran interface
void swaitgrp_   (HTHREADINFO ti,int *state);
void swaitgrpr_  (HTHREADINFO ti,int *state);
void ssetgrp_    (HTHREADINFO ti,int *state);
void gwaitsubs_  (HTHREADINFO ti,int *state);
void gwaitsubsr_ (HTHREADINFO ti,int *state);
void gsetsubs_   (HTHREADINFO ti,int *state);
void gwaitmain_  (HTHREADINFO ti,int *state);
void gwaitmainr_ (HTHREADINFO ti,int *state);
void gsetmain_   (HTHREADINFO ti,int *state);
void waitstate_  (int *state    ,int *nstep);
void waitstater_ (int *state    ,int *nstep);
void setstate_   (int *state    ,int *nstep);
void tscinit(HTHREADINFO ti);
void tscb(HTHREADINFO ti,int id);
void tsce(HTHREADINFO ti,int id);
void tsceb(HTHREADINFO ti,int id);
void prtsc(HTHREADINFO ti,const char*tag);
__END_DECLS
