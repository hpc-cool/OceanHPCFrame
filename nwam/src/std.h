#ifndef STD_H_INCLUDED
#define STD_H_INCLUDED
#include "wavedef.h"
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#define DBGINF printf("%s %d\n",__FILE__,__LINE__)
#define DBGO0 if(mpi_id==0){printf("%s %d %2.2d\n",__FILE__,__LINE__,mpi_id);fflush(stdout);}
#define DBG0  if(mpi_id==0)
#ifdef USEHBM 
#define USE_C_ALLOC 
#include <hbwmalloc.h>
#include <memkind.h>
// /pacific_fs/HPCKit/third/include/hbwmalloc.h
#define hmalloc(n) zhmalloc((n))
#define HMALLOC(n) zhmalloc((n))
#else
#define hmalloc(n) malloc(n)
#define HMALLOC(n) malloc((n))
#endif
#define dmalloc(n) malloc(n)
#define NEW(T)     (T*)malloc(sizeof(T)) 
#define NEWN(T,N)  (T*)malloc(((long)sizeof(T))*(N))
#define NEWPN(p,N) (typeof(p))malloc(((long)sizeof(*p))*(N))
#define NEWZN(p,N) p=(typeof(p))malloc(((long)sizeof(*p))*(N));if(p)memset(p,0,sizeof(*p)*(N))

#define HNEW(T)     (T*)hmalloc(sizeof(T)) 
#define HNEWN(T,N)  (T*)hmalloc(((long)sizeof(T))*(N))
#define HNEWPN(p,N) (  typeof(p))hmalloc(((long)sizeof(*p))*(N))
#define HNEWZN(p,N) p=(typeof(p))hmalloc(((long)sizeof(*p))*(N));if(p)memset(p,0,sizeof(*p)*(N))

#define CPYN(d,s,N)  memcpy(d,s,sizeof(*d)*(N))
#define FILLN(p,v,N) memset(p,v,sizeof(*p)*(N))
#define ZERON(p,N)   memset(p,0,sizeof(*p)*(N))

#define CONSTKGRID
#define LOGSCURR 0
#ifdef CONSTKGRID
#ifndef KLPAT
# define KLPAT NKP
#endif
# define kl CKL
# define klp (CKL+KLPAT)
# define kld CKL
# define jnthet   CJNTHET
# define mkj     (kl*jnthet)       //25*12 kl*jnthet
# define mkjp    (klp*jnthet)      //25*12 kl*jnthet
# define mkjd    (kld*jnthet)      //31*12 kl*jnthet
# define mkjpat  (klp*jnthet)      //40*12 kl*jnthet
# define MBVDEP   50
# define MBVDEPV (MBVDEP+3)/4
# define SIGMALVL 0
#define mkjpat1  (mkjpat+1)
#endif
struct fdim{
  long dimb;//1
  long ndim;//0x80
  long vf1;//1
  long vf2;//0
  long sdim;// 1
  long dime;//0x7f
};
struct fpointer{
  //void *p;//400026200010  71ef6c0
  long vt1; //0
  long vt2; //23
  long ndim;//2
  long vt4; //1b
  long vt5; //4
  long vt6; //10002 20010002
  long vt7; //sizeall 0x80
  long vt8; //sizeall 0x80
  long vt9; //0x81 -0xfe
  void *pb; // 71ef6c0
  long vt11;//400020075a40 400024765a40 400024765a40
  //struct fdim dims[1];//ndim
};
#define FP(n,ND)  struct fpointer FP##n;struct fdim FD_##n[ND]
struct AVHIST{
      int NN,itype,NTPF,mean,hist_ect;
      double nextTime;
      int oldID;
      int nea,itpf;
      int istime;
      int wrwa,ncidwa;
      int wrbv,ncidbv;
      int Logbv,LogWind,Loghs,Logtp,Logtz,Logth,Logau;
      long waoff,bvoff;
      char desc[256],ftimestr[256],wafn[256],bvfn[256];
      char OPTION[8],TAG[8];
      double *ea;FP(ea,3);
      float *aet ;FP(aet ,1);
      float *tpf ;FP(tpf ,1);
      float *h1_3;FP(h1_3,1);
      float *ape ;FP(ape ,1);
      float *bv  ;FP(bv  ,2);
};
typedef struct _windvs{
  float wx,wy,wv,wi;
}windvs;
typedef struct _ImplschP{
  double ws    [CKL]; //sqrt(g*wkk*tanhdk)  |256 aligned 128
  double iwsdk [CKL]; //(1/wsk)*dwk         |256 aligned 128 used mean2
  double ssbo  [CKL]; //ssbo                |256 aligned 128
  //double ccg   [CKL]; //cg                  |256 aligned 128
#if LOGSCURR==1
  double ccgd  [CKL]; //cg*wkk/wsk        |256 aligned 128; used sscu when use current
#endif
  windvs wvs;
  double depth,tRsd_tanLat;//32
  //float ucurc[8];
}ImplschP;
typedef struct _propinf{
  float dep;
  int   i8[9];
}propinf;
typedef struct {
  float pab[4];
  float pabp[4];
  short kpos[4],iquad,fill[3];
}p_interg;
typedef struct _prggpar{
  int nwps,nwpw,nwpc,nwpa;
  int LXB,LYB,LXN,LYN;
  int nbvdep;
  int mpi_id,mpi_comm_wav,mpi_npe;
  windvs constwind;
  int   uconstwind;
  double *wav_ec;
  double *wav_et;
  //double *wav_eco;
  //double *wav_eto;

  ImplschP*ip;

  double *rtcb;
  double *rtcl;
  p_interg *pis;  // (mkj,0:nwpc)
  int *ipos12; // (12,0:nwpc)
  propinf*ipos8; // ( 9,0:nwpc)
  int *nsp;
  int *ieind;
  struct Iepos*iepos;
  float *dep;
  windvs *pwvs;
  double*pvkdp;
  double*pvws;
  double*pvkdo, *pebdep;
  float*h1_3, *ape , *tpf , *aet , *bv  ;

}prggpar;
# define C_SetGPar            c_setgpar_
# define C_Allocate_Vee       c_allocate_vee_
__BEGIN_DECLS
extern prggpar gppar;
extern int mpi_id,NCorePClu ,NThPClu ,NGrpPNode ,NProcPNode,ManageCoreId;
//void C_SetGPar(FC_PARA *fcp,struct Iepos*iepos,int*ieind,int*nsp);
#ifdef C_CALCULATE
extern void InitCPropgats();
#endif
#ifdef USE_C_ALLOC //must define USE_C_ALLOC in platform_init_mod.F90
void C_Allocate_Vee(int *res);
#else
void c_setpointers_(double*eec,double*eet,ImplschP*ip,p_interg*pp_vs,int*ipos12,propinf*ipos8,windvs*pwvs);
#endif
void init_cpu_(int *mpi_id_,int *mpi_npe_,int*mpi_comm_,int * NCorePClu_ ,int *NThPClu_ ,int *NGrpPNode_ ,int *NProcPNode_,int *ManageCoreId_);
void*zhmalloc(size_t size);
__END_DECLS
#endif
