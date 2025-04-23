#ifndef SVARS_H_INCLUDED
#define SVARS_H_INCLUDED
#include "wavedef.h"
#include "std.h"
#include "zarm.h"

#define DBGLN printf("%s %d",__FILE__,__LINE__)
#define DBGL printf("%s %d\n",__FILE__,__LINE__)

typedef union _prgstate{
  double vlt;
  struct {
    long   state;
    long   ntimestep;
    long   slave_state;
    long   tstate;
  };
}prgstate;

#define ZBYTE signed char
#define  IMDPS   4
typedef struct _ImplschPar{
  double  wkal[4*12]  ;      //0x180
  long    kjps[CJNTHET*24];  //0x900 (  2,12,jn)
  ZBYTE    kjs[CJNTHET*48];  //0x900 (2,2,12,jn)
  double cwks17[CKL];        //0x100
  double grolim[CKL];        //0x100
  double wk[CKL],wkh[CKL];   //0x200
  double dwk[CKL],wkdk[CKL],wkibdk[CKL];     //0x300
#if SIGMALVL==2
  double expdep[CKL*(MBVDEP)];
#endif
  double  cosths[CJNTHET],sinths[CJNTHET]  ;     //0xc0
  double  lmdpd,lmdmd,lmdpd20,lmdmd20,Dm2,D1;
  float   beta10,beta11;        //0x88 beta11
  float   ads,brkd1,brkd2,deltts;
  float spdeltts;
  int   ntsplit;

  double  sds,awk,enh;
  double  zpi;


  float   D1m30,F44_24778,F16_0,F0_467;//mean3
  float   F18_20832,F0_001930368,F0_5,F1_0;//mean2
  float   F0_75,F5_5,F0_833,Fm1_25;//enh
  float   F28_,F0_80,F0_065,F0_001;
  float   tztp,tzzpi;

  ImplschP*pip;
  //p_interg *pis;                    //     (mkj,0:nwpc)
  //int *ipos12;                      //0x10 (12,0:nwpc)
  //propinf *ipos8;                       //0x18 (9,0:nwpc);
  //double *wav_ec,*wav_et;           //0x0
  //windvs *pwvs;
  //float   *wx,*wy,*wiv,*winc; //0x28
  int iacb,iace;
  int mpi_id,tid,dbglvl,dmalog;     //0x138
  int nbvlvl,nwps,nwpc,nwpa;     //0x148
  //int iacb1,iace1,iacb2,iace2,iacb3,iace3;
  //int fill0[2] ;//0x12d8 aligned to 0x80,size 0x1300
}ImplschPar;

typedef struct _gridpara{
  long _kl,_jnthet,_kld,_NBVDEP,_SIGMALVL;
}gridpara;

typedef struct _ImplschLV{
  double se[mkj];
  double ee[mkj];//aligned 256 for dma
  double eer[mkj];//aligned 256 for dma
  ImplschP pi; //256*3+32 aligned 128 for dma
  void*plv ,*plvn,*pee ,*pip ,*peet;
  long iac;
  long fill[10];
  //int  reply_gp,reply_ge,reply_pe;
#ifdef IMPDBGDATA
  //  double vv[mkj];//aligned 256 for dma
#endif
}ImplschLV;


# define NPPGBLOCK 32
# ifdef P_NVEC
#   undef  P_NVEC
# endif
typedef struct _PPGBUF{
#if P_NVEC==1
  double pe[mkj];  //0xc00
#else
  double  pe[mkj]; //11
#endif
  propinf ipos8;
  //int l_pos[12];
  int fill[22];
}PPGBUF;
typedef struct PPG_LV{
#if P_NVEC==1
  double pet[10*mkj];
#else
  double  pet[10*mkj]; //11
#endif
  p_interg pib; //0x2700
  PPGBUF pbuf[2];
}PPG_LV;

typedef struct _OutputLV{
  double ee[mkj];//aligned 256 for dma
  float ape,tpf,aet,h1_3;
  float bv[MBVDEP];
}OutputLV;

#ifdef MTHREAD
# define C_ctStart            c_ctstart_
# define C_ctEnd              c_ctend_
//# define C_waitslave          c_waitslave_
//# define C_setstate           c_setstate_
//# define C_Wait               c_wait_
#endif
//void checksp(gdata *pgpp,int ll);
#define CHECKESP  //checksp(pgpp,__LINE__);
#define CHECKMEM  //CheckMem(pgpp,__LINE__);

typedef struct _Boundary{
  double cosths[CJNTHET],sinths[CJNTHET];
  double rwk[CKL];
  int nwpc;
  double rg,g,gama,zpi,windfield,xj0,xj_,arlfa_,wsj_;
  double rsigma1,rsigma2,wv0;
  double F_1,F_m1d25,F_m0d5,F_m0d4,F_m0d33;
  double* pvws;
  int*nsp;
}Boundary;
struct v_interg{
  double pp0[8],pp1[8],pp2[8],pp3[8];//64*4
  double ps0[8],ps1[8],ps2[8],ps3[8];//64*4
  short  kpos[4][8];                 //16*4=64
  int    offset[8];                  // 8*4=32
  int    iquad,flag,ipos[3],fill[3]; // 8*4=32
};
struct g2lInd{
  int ipart,iac;
};
#include "gpart.h"
#ifdef MTHREAD
#define TUSERDECLARE struct{ \
  double *wav_ec, *wav_et;\
  int iacb1,iace1,iacb2,iace2;\
  int nstep;\
}
#define NUSERDECLARE struct{\
  int hist_eot,stop_now,rest_eot;\
  int nwps,nwpw,nwpc,nwpa;\
  int LXB,LYB,LXN,LYN;\
  int mpi_id,mpi_npe,mpi_comm_wav;\
  struct g2lInd*g2l;\
  windvs constwind;\
  int   uconstwind;\
  double *wav_ec, *wav_et;\
  int nstep;\
  int *ieind;\
  struct Iepos*iepos;\
  sendcomm s_segs;\
  int *nsp;\
  float *dep;\
  windvs *pwvs;\
  Boundary bdy;\
  p_interg *pis;\
  /*struct v_interg *pvis;*/\
  /**/\
  int *ipos12; \
  propinf*ipos8;\
  ImplschPar pisp;\
  int ndep,npdep;\
  double *pvkd,*pebdep;\
  float *ape , *tpf , *aet , *h1_3, *bv;\
}
#define GUSERDECLARE struct{ \
  int hist_eot,stop_now,rest_eot;\
  int nwpt,nwps,nwpc,nwpo,nwpa;\
  int NPART,block_id,LXB,LYB,LXN,LYN;\
  Rect recti,recto;\
  int *l2gind;\
  char*fbuff;\
  sendcomm s_segs;\
  int mpi_id,mpi_npe,mpi_comm_wav;\
  double *wav_ec, *wav_et;\
  int nstep;\
  /*inited at g_initSetspec*/\
  char *nsp;\
  /*float *dep;*/\
  Boundary bdy;\
  /*inited at g_initPropgats OK*/\
  /*p_interg *pis; */\
  struct v_interg *pvis;\
  /*int *ipos12; */\
  /*propinf*ipos8;*/\
  /*inited at C_Setwind OK*/\
  /*inited at g_initImplsch OK*/\
  ImplschPar pisp;\
  /*inited at g_initOutput_cal OK */\
  int ndep,npdep;\
  double *pvkd,*pebdep;\
}

#include "mthread.h"
void initmd();
#endif
__BEGIN_DECLS
void c_fillsendcomm_();
// in propagat.cpp
extern void f_initGPar(threadGroup *gi);
extern void g_initPropgats(threadGroup *gi);
extern void g_initImplsch(threadGroup *gi);
extern void g_initSetspec(threadGroup *gi);
// in svars.c
void g_setoutput_(threadGroup *gi);
extern ImplschPar pisp;

# define C_SetGPar            c_setgpar_
# define C_InitImplsch        c_initimplsch_
# define C_SetPropInterg      c_setpropinterg_
# define C_Setwind            c_setwind_

//void C_SetGPar(FC_PARA *fcp);
void C_InitImplsch(ImplschPar *pisp_,ImplschP*pip,float*pwvs);
void C_SetPropInterg(int *iac_,int*j_,int*k_,int*iquad_,int*kpos,float*pab,float*pabp);
void C_Setwind();
extern void setspec(threadGroup *gi,double *e,int nwpb,int nwpe,int n) ;
extern void propagats(threadGroup *gi,double *wav_e,double *wav_ee,int iacb,int iace) ;
extern void implschs (threadGroup *gi,double *eec,double *eet,int iacb,int iace) ;
extern void implschs_(threadGroup *gi,double *eec,double *eet,int *piacb,int *piace) ;
extern void accumea_(int *ind,double*ee,int*nwpb,int *nwpe);
extern void checkoutput_(HTHREADINFO ti,int*ind,int*RFB,double*ee,int*iacb1,int*iace1,int*iacb2,int*iace2,int*nout);
extern void t_endoutput_cal_(HTHREADINFO ti,int*iacb1,int*iace1);
extern void g_initOutput_cal(threadGroup *gi);
__END_DECLS
#endif //SVARS_H_INCLUDED
