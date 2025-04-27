#define _GNU_SOURCE
#include "std.h"
#include <linux/types.h>
//#include <linux/atomic.h>
#include <mpi.h>
#include <sched.h>
#include <pthread.h>
prggpar gppar;
int indg,grpind,ngrps,mpi_id,mpi_npe,grp_npe;
int NCorePClu=38 ,NThPClu=38 ,NGrpPNode=16 ,NProcPNode=1,ManageCoreId=-1;
MPI_Comm mpi_comm;
#define atomic_t int
struct hbmstate{
  atomic_t lock;
  size_t hbmused,HbmMax;
};
static int hbmlock=0;
struct hbmstate hbms[16]={0};
static int getcpuid(){
  cpu_set_t mask;  //CPU核的集合
  CPU_ZERO(&mask);    //置空
  if (sched_getaffinity(0, sizeof(mask), &mask) == -1){//设置线程CPU亲和力
    printf("warning: could not get CPU affinity , continuing...\n");
    return -1;
  }
  for(int i=0;i<CPU_SETSIZE;i++){
    if(CPU_ISSET(i,&mask)){
      return i;
    }
  }
  return 0;
}
static int bindcpu(int id){
  cpu_set_t mask;  //CPU核的集合
  CPU_ZERO(&mask);    //置空
  CPU_SET(id,&mask);   //设置亲和力值
  printf("Bind CPU %d\n",id);
  if (sched_setaffinity(0, sizeof(mask), &mask) == -1){//设置线程CPU亲和力
    printf("warning: could not set CPU affinity %d, continuing...\n",id);
    return -1;
  }
  return 0;
}
static void ntdelay(int n){ static int vt=0; for(int i=0;i<n;i++) vt=(vt*1357+2581); }
void zunlock(atomic_t *p){ *p=0;}
void zlock(atomic_t *p){ while(__sync_val_compare_and_swap(p,0,1))ntdelay(100);}
#ifndef NO_MPI
#define DBG printf("%s %d %3.3d\n",__FILE__,__LINE__,mpi_id);fflush(stdout);MPI_Barrier(mpi_comm);
#else
#define DBG printf("%s %d %3.3d\n",__FILE__,__LINE__,mpi_id);fflush(stdout);
#endif
#define ADBG printf("========= HBM %d %3.3d %3.2d %d %d\n",__LINE__,mpi_id,indg,grpind,grp_npe);fflush(stdout)
void*zhmalloc(size_t size){
#ifdef USEHBM
  int cpuid=getcpuid();
  struct hbmstate *ph=&hbms[cpuid/NCorePClu];
  void*p=NULL;
#ifndef MMTHREAD
#ifndef NO_MPI
  MPI_Status status;
  size_t tta=0;
  //for node serial run
  if(grpind&&indg==0) MPI_Recv(&tta,1,MPI_LONG,(grpind-1)*NThPClu,10001,mpi_comm,&status);
  if(indg){
    MPI_Recv(&ph->hbmused,1,MPI_LONG,mpi_id-1,10001,mpi_comm,&status);
  }
#endif
#endif
#ifdef MTHREAD
  zlock(&hbmlock);
#endif
  if(ph->HbmMax==0){
    ph->HbmMax=0xFFB00000L;
    ph->hbmused=0;
  }
  if(ph->hbmused+size>ph->HbmMax){//fixme
    printf("HBM warning:nospace mpi_id:%d %lX %lX\n",mpi_id,ph->hbmused,size);
    p=malloc(size);
  }else{
    //printf("Bhmalloc %d %d %lx %lx \n",mpi_id,cpuid,ph->hbmused,size);
    p=memkind_malloc(MEMKIND_HBW,size);
    size_t st=memkind_malloc_usable_size(MEMKIND_HBW, p);
    ph->hbmused+=st;
    //printf("Ehmalloc %d %d %lx %lx %lx\n",mpi_id,cpuid,ph->hbmused,size,st);
  }
#ifdef MTHREAD
  zunlock(&hbmlock);
#endif
#ifndef MMTHREAD
#ifndef NO_MPI
  if(indg<grp_npe-1){
    MPI_Send(&ph->hbmused,1,MPI_LONG,mpi_id+1,10001,mpi_comm);
  }else{
    MPI_Send(&ph->hbmused,1,MPI_LONG,mpi_id-grp_npe+1,10001,mpi_comm);
  }
  if(indg==0){
    MPI_Recv(&ph->hbmused,1,MPI_LONG,mpi_id+grp_npe-1,10001,mpi_comm,&status);
    //for node serial run
    if(grpind<ngrps-1) MPI_Send(&tta,1,MPI_LONG,(grpind+1)*NThPClu,10001,mpi_comm);
  }
#endif
#endif
  //ADBG;
  return p;
#else 
  return malloc(size); 
#endif
}
void init_cpu_(int *mpi_id_,int *mpi_npe_,int*mpi_comm_,int * NCorePClu_ ,int *NThPClu_ ,int *NGrpPNode_ ,int *NProcPNode_,int *ManageCoreId_){
  mpi_id        =*mpi_id_;
  mpi_npe       =*mpi_npe_;
  mpi_comm      =MPI_COMM_WORLD;//*mpi_comm_;
  NCorePClu     =*NCorePClu_;
  NThPClu       =*NThPClu_;
  NGrpPNode     =*NGrpPNode_;
  NProcPNode    =*NProcPNode_;
  ManageCoreId  =*ManageCoreId_;
  if(ManageCoreId <=NGrpPNode*NCorePClu)ManageCoreId  =0; 
#ifndef MTHREAD
  int nthpn=NThPClu*NGrpPNode;
  int nid=mpi_id/nthpn;
  int mid=mpi_id%nthpn;
  grpind=mid/NThPClu;
  int nn=mpi_npe/nthpn;
  ngrps=NGrpPNode;
  if(nid>=nn){
    ngrps=(mpi_npe%nthpn+NThPClu-1)/NThPClu;
  }
  indg=mid%NThPClu;
  grp_npe=NThPClu;
  int nt=(mpi_npe/NThPClu)*NThPClu;
  if(nt<=mpi_id){
    grp_npe=mpi_npe%NThPClu;
  }
  int ic=indg+grpind*NCorePClu;
  bindcpu(ic);
#endif
}
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
void prtieind(struct Iepos*iepos,int*ieind,int*nsp);
__BEGIN_DECLS
#ifdef MMTHREAD
extern void InitThreadsGPar(FC_PARA *fcp,struct Iepos*iepos,int*ieind,int*nsp);
#endif
__END_DECLS
void C_SetGPar(FC_PARA *fcp,struct Iepos*iepos,int*ieind,int*nsp){
  int acwv;
  //printf("AAA %d %d\n",fcp->nwpc,fcp->nwpa);
  gppar.mpi_id=fcp->mpi_id;
  gppar.mpi_comm_wav=fcp->mpi_comm_wav;
  gppar.mpi_npe=fcp->mpi_npe;
  //pisp.mpi_id=fcp->mpi_id;
  gppar.nwps  =fcp->nwps  ;
  gppar.nwpc  =fcp->nwpc  ;
  gppar.nwpa  =fcp->nwpa  ;
  gppar.nbvdep=fcp->_NBVDEP;
  gppar.nwpw  =0;
  gppar.LXB   =fcp->LXB  ;
  gppar.LYB   =fcp->LYB  ;
  gppar.LXN   =fcp->LXN  ;
  gppar.LYN   =fcp->LYN  ;
  gppar.nsp   =nsp;
  gppar.ieind =ieind;
  gppar.iepos =iepos;

  gppar.constwind.wx=fcp->constwindx;
  gppar.constwind.wy=fcp->constwindy;
  gppar.constwind.wv=gppar.constwind.wx*gppar.constwind.wx+gppar.constwind.wy*gppar.constwind.wy;
  gppar.constwind.wi=0;
  acwv=gppar.constwind.wv*100;
  if(acwv>=10){
    gppar.uconstwind=2;
  }else if(acwv>=1){
    gppar.uconstwind=1;
  }else {
    gppar.uconstwind=0;
  }
#ifdef MMTHREAD
  InitThreadsGPar(fcp,iepos,ieind,nsp);
#endif
  //prtieind(iepos,ieind,nsp);
}

//#define USE_C_ALLOC
//#ifndef MMTHREAD
#ifdef USE_C_ALLOC //must define USE_C_ALLOC in platform_init_mod.F90
void f_setevar_  (double*,double*);
void f_setimpvar_(ImplschP*ip,double*pvws,windvs *pwvs);
void f_setprop_  (p_interg*pis,double*pvkdp,propinf*ipos8,int *ipos12);
void f_setoutput_(double*cpvkdo,double*cpebdep);
static char *fbuff=NULL;
void c_allocate_vee_(int *res){
  char*buf;
  long size;
  //printf("c_alloc %d %d\n",gppar.nwpc,gppar.nwpa);
#define IMDPSP  2
#define IMDPSO  4
#define BVIMDPS 1
#define VAVARS \
  VVAR(pebdep ,double  ,(gppar.nwpc+1)*kld*(BVIMDPS+gppar.nbvdep));\
  VVAR(wav_ec ,double  ,(gppar.nwpa+1)*mkj);\
  VVAR(wav_et ,double  ,(gppar.nwpa+1)*mkj);\
  VVAR(ip     ,ImplschP,(gppar.nwpc+1)    );\
  VVAR(pvws   ,double  ,(gppar.nwpc+1)*kl );\
  VVAR(pvkdp  ,double  ,(gppar.nwpc+1)*kl *IMDPSP);\
  VVAR(ipos12 ,int     ,(gppar.nwpc+1)*12 );\
  VVAR(ipos8  ,propinf ,(gppar.nwpc+1)    );\
  VVAR(pis    ,p_interg,(gppar.nwpc+1)*mkj);\
  VVAR(pwvs   ,windvs  ,(gppar.nwpc+1)    );\
  VVAR(pvkdo  ,double  ,(gppar.nwpc+1)*kld*IMDPSO);

  size=0;
#define  VVAR(vn,T,N) size+=sizeof(T)*(N)+256
  VAVARS
#undef VVAR
#ifndef MMTHREAD
  buf=fbuff=HNEWN(char,size);
#ifdef MTHREAD
  int ires=hbw_verify_memory_region(buf, size,  1);
  if(ires){
    printf(" Hbm alloc error %ld %p============ \n",size,buf);
  } 
  //else printf(" Hbm alloc ok %ld %p============ \n",size,buf);
#endif
#else
  buf=fbuff=NEWN(char,size);
  printf("C_alloc no hbm %d %lX\n",mpi_id,size);
#endif
  if(buf==NULL){*res=-1;return ;}
  buf=(char*)(( ((long)buf)+255)&(-256) );
#define  VVAR(vn,T,N) gppar.vn=(T*)buf;buf=(char*)( (((long)buf)+sizeof(T)*(N)+255) & (-256) )
  VAVARS
#undef VVAR
#undef VAVARS
  f_setevar_(gppar.wav_ec,gppar.wav_et);
  f_setimpvar_(gppar.ip ,gppar.pvws ,gppar.pwvs);
  f_setprop_(gppar.pis,gppar.pvkdp  ,gppar.ipos8,gppar.ipos12);
  f_setoutput_(gppar.pvkdo  ,gppar.pebdep );
  *res= 0;
}
#endif
void cpyeec_(int *mode){
#ifdef USE_C_ALLOC 
  if(*mode==0){
    memcpy(gppar.wav_ec+mkj,gppar.wav_ec+mkj,gppar.nwps*mkj*sizeof(double));
  }else{
    memcpy(gppar.wav_ec+mkj*(gppar.nwpc+1),gppar.wav_ec+mkj*(gppar.nwpc+1),(gppar.nwpa-gppar.nwpc)*mkj*sizeof(double));
  }
#endif
}
void c_setpointers_(double*eec,double*eet,ImplschP*ip,p_interg*pp_vs,int*ipos12,propinf*ipos8,windvs*pwvs  ){
#if 0
#ifdef USE_C_ALLOC 
  gppar.wav_eco=eec   ;
  gppar.wav_eto=eet   ;
#else
  gppar.wav_ec =eec   ;
  gppar.wav_et =eet   ;
#endif
#else
  gppar.wav_ec =eec   ;
  gppar.wav_et =eet   ;
#endif
  gppar.ip     =ip    ;
  gppar.pis    =pp_vs ;
  gppar.ipos12 =ipos12;
  gppar.ipos8  =ipos8 ;
  gppar.pwvs   =pwvs  ;
}
void initmtpar();
void c_init_distinct_(){
#ifdef MTHREAD
  initmtpar();
#endif
#ifdef C_CALCULATE
  InitCPropgats();
#endif
}
//#endif
int GetVInt(int volatile *volatile p){
	return *p;
}
