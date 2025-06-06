#include "svars.h"
// Main Core Code
extern "C" void offsetinit();

static inline uint64_t arm64_cntvct(void) {
	uint64_t tsc;
    asm volatile("mrs %0, cntvct_el0" : "=r" (tsc));
    return tsc;
}

static inline uint64_t arm64_cntfrq(void) {
	uint64_t freq;
	asm volatile("mrs %0, cntfrq_el0" : "=r" (freq));
	return freq;
}

static inline uint64_t rdtsc(void) {
	return arm64_cntvct();
}
void ppa();
void ppb();
//static struct v_interg *pvis=NULL;
//static double*ppab,*ppabp,*ppabt=NULL;
/*
 *svfloat64_t {p0,p1,p2,p3} ...
 * not surport load from float32_t to svfloat64_t ,
 * speed for cache
 * */
//#define DBGINF if(iac>=292)printf("%s,%d,%d,%d,%d\n",__FILE__,__LINE__,iac,iace,kj)
#ifdef DBGINF
#undef DBGINF
#endif
#define DBGINF //printf("%s %4d %04d ",__FILE__,__LINE__,pisp.mpi_id)
#define DBGINFE //printf("%s %4d %04d\n",__FILE__,__LINE__,pisp.mpi_id)
void InitCPropgats(){
}
#define DBG //fprintf(gi->threads[0].fo,"CC %d %2.2d \n",__LINE__,gi->igrp);fflush(gi->threads[0].fo);
#define DBGA //fprintf(gi->threads[0].fo,"CCC %d %2.2d %4.3d/%4d %2.2d\n",__LINE__,gi->igrp,iac,gi->nwpc,kj);fflush(gi->threads[0].fo);
void g_initPropgats(threadGroup *gi){
  //ppa();
  //ppb();
  int iac;
  struct v_interg *pvis;
  iac=0;
  //printf("CCA %d %2.2d %4d/%4d %ld\n",__LINE__,gi->igrp,iac,gi->nwpc,sizeof(struct v_interg)*(gi->nwpc+1)*mkj/8);fflush(stdout);
  gi->pvis=pvis=HNEWN(struct v_interg ,(gi->nwpc+1)*mkj/8);

  for(iac=1;iac<gi->nwpc;iac++){
    int it=(iac*mkj),*ii;
    int kj,iquad=-1;
    int giac=gi->l2gind[iac];
    p_interg *pi=md.pis+giac*mkj;  // (mkj,0:nwpc)
    int *gl_pos=md.ipos12+giac*12;
    int l_pos[12]={0};
    for(kj=0;kj<12;kj++){
      int mt=l_pos[kj]=md.g2l[gl_pos[kj]].iac;// assume part is ok
      if(mt<0||mt>gi->nwpa){
        l_pos[kj]=0;//fixme
      }
    }
    kj=0;
    struct v_interg*pv=pvis+it/8;
    for(kj=0;kj<mkj/8;kj++,pi+=8){
      int i,k;
      int idf=0,igf=0;
      iquad=pi->iquad;
      if(1){
        int *ii=l_pos+pi->iquad*3;
        int i,j,it[4],ie[4];
        ie[0]=iac;
        for(i=0;i<4;i++){
          it[i]=pi[0].kpos[i];
          int *ii=l_pos+pi[0].iquad*3;
          if(i)ie[i]=ii[i-1];
          for(j=1;j<8;j++){
            int *ii=l_pos+pi[j].iquad*3;
            int ig=0, id=pi[j].kpos[i]-it[i]-j;
            if(i) ig=ii[i-1]-ie[i];
            if(id){
              if(pi[j].pab[i]<1e-5)id=0;
            }
            idf|=id;igf|=ig;
          }
        }
        pv[kj].flag= (idf||igf);
      }else pv[kj].flag= 0;
      pv[kj].iquad  =pi[0].iquad;
      pv[kj].ipos[0]=l_pos[iquad*3+0];
      pv[kj].ipos[1]=l_pos[iquad*3+1];
      pv[kj].ipos[2]=l_pos[iquad*3+2];
      {
        for(int im=0;im<3;im++){
          int mt=pv[kj].ipos[2];
          if(mt<0||mt>gi->nwpa){
            printf("EEEEEEEEE %d %d %5d %5d\n",mpi_id,gi->igrp,mt,gi->nwpa);
            exit(0);
          }
        }
      }
      for(i=0;i<8;i++){
        pv[kj].pp0[i]=pi[i].pab [0];
        pv[kj].pp1[i]=pi[i].pab [1];
        pv[kj].pp2[i]=pi[i].pab [2];
        pv[kj].pp3[i]=pi[i].pab [3];
        pv[kj].ps0[i]=pi[i].pabp[0];
        pv[kj].ps1[i]=pi[i].pabp[1];
        pv[kj].ps2[i]=pi[i].pabp[2];
        pv[kj].ps3[i]=pi[i].pabp[3];
        pv[kj].kpos[0][i]=pi[i].kpos[0];
        pv[kj].kpos[1][i]=pi[i].kpos[1];
        pv[kj].kpos[2][i]=pi[i].kpos[2];
        pv[kj].kpos[3][i]=pi[i].kpos[3];
        if(idf||igf){//fixme
          //pv[kj].offset[]  =0;
        }
      }
    }
    //printf(" %d %d %d\n",iac,kj,pvis[384/8].iquad[0]);
  }
}
#define LDO(P,l) svld1_f64(pm,P+kpos[l*8])
void propagats(threadGroup *gi,double *wav_e,double *wav_ee,int iacb,int iace){
  int iac;
  wav_ee--; //kpos start from 1,fortran
  svbool_t pm=svptrue_b64();
  for(iac=iacb;iac<iace;iac++){
    if(gi->nsp[iac]!=1)continue;
    int it=(iac*mkj);
    int kj,iquad=-1;
    struct v_interg *pi=gi->pvis+it/8;  // (mkj,0:nwpc)
    double *pe0,*pe1,*pe2,*pe3,*pe;
    pe0=wav_ee+it;
    pe=wav_e+it;
    for(kj=0;kj<mkj;kj+=8,pi++){
      if(!pi->flag){
        if(iquad!=pi->iquad){
          pe1=wav_ee+pi->ipos[0]*mkj;
          pe2=wav_ee+pi->ipos[1]*mkj;
          pe3=wav_ee+pi->ipos[2]*mkj;
        }
        short *kpos=&pi->kpos[0][0];
        VDOUBLE et,et1;
        VDOUBLE  p0,p1,p2,p3;
        p0=LD(pi->pp0 ); p1=LD(pi->pp1 ); p2=LD(pi->pp2 ); p3=LD(pi->pp3 );
        et =MUL(   LD(pi->ps0 ),(MLA(MLA(MLA(MUL(p0,LDO(pe0,0)),p1,LDO(pe0,1)),p2,LDO(pe0,2)),p3,LDO(pe0,3))));
        et =MLA(et,LD(pi->ps1 ),(MLA(MLA(MLA(MUL(p0,LDO(pe1,0)),p1,LDO(pe1,1)),p2,LDO(pe1,2)),p3,LDO(pe1,3))));
        et1=LDO(pe2,0);
        et1=MUL(p0,et1);
        et1=MLA(et1,p1,LDO(pe2,1));
        et1=MLA(et1,p2,LDO(pe2,2));
        et1=MLA(et1,p3,LDO(pe2,3));
        et =MLA(et,LD(pi->ps2 ),(et1));
        //et =MLA(et,LD(pi->ps2 ),(MLA(MLA(MLA(MUL(p0,LDO(pe2,0)),p1,LDO(pe2,1)),p2,LDO(pe2,2)),p3,LDO(pe2,3))));
        et1=LDO(pe3,0);
        et1=MUL(p0,et1);
        et1=MLA(et1,p1,LDO(pe3,1));
        et1=MLA(et1,p2,LDO(pe3,2));
        et1=MLA(et1,p3,LDO(pe3,3));
        et =MLA(et,LD(pi->ps3 ),(et1));
        //et =MLA(et,LD(pi->ps3 ),(MLA(MLA(MLA(MUL(p0,LDO(pe3,0)),p1,LDO(pe3,1)),p2,LDO(pe3,2)),p3,LDO(pe3,3))));
        VVD(pe+kj)=et;
      }else{
        // fixme...
      }
#undef LDO
    }
  }
}

void propagats_(threadGroup *gi,double *wav_e,double *wav_ee,int *piacb,int *piace){
  propagats(gi,wav_e,wav_ee,*piacb,*piace);
}
#if 0
void propagats_naive(double *wav_e,double *wav_ee,int *piacb,int *piace){
  int iac,iacb=*piacb,iace=*piace;
  /*printf("iacb=%d,iace=%d\n",iacb,iace);*/
  wav_ee--; // kpos is index from 1
  for(iac=iacb;iac<iace;iac++){
    int it=(iac*mkj);
    int kj,iquad=-1;
    p_interg *pi=md.pis+it;  // (mkj,0:nwpc)
  	int *l_pos=md.ipos12+iac*12;
		double *pe0,*pe1,*pe2,*pe3,*pe;
    pe0=wav_ee+it;
    pe=wav_e+it;
    for(kj=0;kj<mkj;kj++,pi++){
      if(iquad!=pi->iquad){
        int iqt=(iquad=pi->iquad)*3;
        int *ii=l_pos+iqt;
        pe1=wav_ee+ii[0]*mkj;
        pe2=wav_ee+ii[1]*mkj;
        pe3=wav_ee+ii[2]*mkj;
      }
      {
				short *kpos=pi->kpos;
				float*pab=pi->pab,*pabp=pi->pabp;
				double et;
        et =pabp[0]*(pab[0]*pe0[kpos[0]]+pab[1]*pe0[kpos[1]]+pab[2]*pe0[kpos[2]]+pab[3]*pe0[kpos[3]]);
        et+=pabp[1]*(pab[0]*pe1[kpos[0]]+pab[1]*pe1[kpos[1]]+pab[2]*pe1[kpos[2]]+pab[3]*pe1[kpos[3]]);
        et+=pabp[2]*(pab[0]*pe2[kpos[0]]+pab[1]*pe2[kpos[1]]+pab[2]*pe2[kpos[2]]+pab[3]*pe2[kpos[3]]);
        et+=pabp[3]*(pab[0]*pe3[kpos[0]]+pab[1]*pe3[kpos[1]]+pab[2]*pe3[kpos[2]]+pab[3]*pe3[kpos[3]]);
        pe[kj]=et;
      }
    }
  }
}
#endif
#if 0
void ppa(){
  int iac,iacb=1,iace=md.nwpc;
  FILE*fl=fopen("ppa.dat","wt");
  fprintf(fl,"ZAAAAAA %d %d \n",iacb,iace);
  int nnvec=0;
  for(iac=iacb;iac<iace;iac++){
    int it=(iac*mkj);
    int kj,iquad=-1;
    int *l_pos=md.ipos12+iac*12;
    p_interg *pi=md.pis+it;  // (mkj,0:nwpc)
    for(kj=0;kj<mkj;kj+=8,pi+=8){
      iquad=pi->iquad;
      int *ii=l_pos+pi->iquad*3;
      int i,j,it[4],ie[4];
      int idf,igf,iie,ije;
      iie=-1,ije=-1;
      idf=0;igf=0;
      ie[0]=iac;
      for(i=0;i<4;i++){
        it[i]=pi[0].kpos[i];
        int *ii=l_pos+pi[0].iquad*3;
        if(i)ie[i]=ii[i-1];
        for(j=1;j<8;j++){
          int *ii=l_pos+pi[j].iquad*3;
          int ig=0, id=pi[j].kpos[i]-it[i]-j;
          if(i) ig=ii[i-1]-ie[i];
          if(id){
            if(pi[j].pab[i]<1e-5)id=0;
          }
          idf|=id;igf|=ig;
          if(id){
            iie=i;ije=j;
          }
        }
      }
      if(idf){
        fprintf(fl,"%d %d %d A: %d %d\n",iac,kj,i,iie,ije);
        for(i=0;i<4;i++){
          fprintf(fl,"%d %d %d B: ",iac,kj,i);
          for(j=0;j<8;j++){
            fprintf(fl,"%8d ",pi[j].kpos[i]);
          }
          fprintf(fl,"\n");
          fprintf(fl,"%d %d %d C: ",iac,kj,i);
          for(j=0;j<8;j++){
            fprintf(fl,"%8.5f ",pi[j].pab[i]);
          }
          fprintf(fl,"\n");
        }
      }
      if(igf){
        for(i=1;i<4;i++){
          fprintf(fl,"%d %d %d D:%d ",iac,kj,i,ie[0]);
          for(j=1;j<8;j++){
            int *ii=l_pos+pi[j].iquad*3;
            fprintf(fl,"%d ",ii[i-1]);
          }
          fprintf(fl,"\n");
        }
      }
      if(idf||igf)nnvec++;
    }
  }
  fclose(fl);
  printf("no vec %d/%d\n",nnvec,(iace-iacb)*mkj/8);
}
void ppb(){
  int iac,iacb=1,iace=md.nwpc;
  FILE*fl=fopen("ppb.dat","wt");
  fprintf(fl,"ZAAAAAA %d %d \n",iacb,iace);
  int nnvec=0;
  for(iac=iacb;iac<iace;iac++){
    int it=(iac*mkj);
    int kj,iquad=-1;
    int *l_pos=md.ipos12+iac*12;
    p_interg *pi=md.pis+it;  // (mkj,0:nwpc)
    for(kj=0;kj<mkj;kj+=8,pi+=8){
      iquad=pi->iquad;
      int *ii=l_pos+pi->iquad*3;
      int i,j;
      fprintf(fl,"%d %d :\n",iac,kj);
      for(i=0;i<4;i++){
        fprintf(fl,"%d %d %d B: ",iac,kj,i);
        for(j=0;j<8;j++){
          fprintf(fl,"%8d ",pi[j].kpos[i]);
        }
        fprintf(fl,"\n");
        fprintf(fl,"%d %d %d C: ",iac,kj,i);
        for(j=0;j<8;j++){
          fprintf(fl,"%8.5f ",pi[j].pab[i]);
        }
        fprintf(fl,"\n");
      }
      for(i=0;i<4;i++){
        fprintf(fl,"%d %d %d D: ",iac,kj,i);
        for(j=0;j<8;j++){
          int *ii=l_pos+pi[j].iquad*3;
          if(i>0) fprintf(fl,"%d ",ii[i-1]);
        }
        fprintf(fl,"\n");
        fprintf(fl,"%d %d %d E: ",iac,kj,i);
        for(j=0;j<8;j++){
          fprintf(fl,"%8.5f ",pi[j].pabp[i]);
        }
        fprintf(fl,"\n");
      }
    }
  }
  fclose(fl);
}
#endif
