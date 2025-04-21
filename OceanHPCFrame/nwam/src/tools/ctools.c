/*
	some tools write in c code,used by fortran code
*/
#define _GNU_SOURCE
#include <stdio.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

#define DBGPINT(a) //printf a
#if ( defined FORTRANCAPS )
#define msleep    MSLEEP
#define IntVolatile  INTVOLATILE
#elif ( defined FORTRANUNDERSCORE )
#define msleep    msleep_
#define IntVolatile  intvolatile_
#elif ( defined FORTRANDOUBLEUNDERSCORE )
#define msleep    msleep__
#define IntVolatile  intvolatile__
#endif
#define PENV(name) printf("ENV %s=%s\n",#name,getenv(#name))
extern char **environ;

static long dclkb=0;
double dclktime(){
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME,&ts);
  if(dclkb==0)dclkb=ts.tv_sec;
  return (ts.tv_sec-dclkb)+ts.tv_nsec/1000000000.;
}
double dclktimeinit(){
  return dclkb=dclktime();
}
void checkenv_(){
	char *s;
	int i;
	PENV(TM_API);
	//for(i=0;environ[i];i++){
	//	printf("ENV %s\n",environ[i]);
	//	fflush(stdout);
	//}
	//sleep(10000);
	//execl("/bin/sh","tenv",0);
}
void msleep(int *ms){
  usleep(*ms*1000);
}
void setmask1_(int *p,int *nn){
	int i;
	for(i=0;i<*nn;i++)p[i]=p[i]>0?1:0;
}
void setgrid_(int *p,int *nn,float*pvb,float *pv0,float *pv1){
	int i;
	float vb,v0,v1,v,vz;
	vb=*pvb;v0=*pv0;v1=*pv1;vz=0.;
	//printf("%p %p %d\n",p,p+1,*nn)
	for(i=0;i<*nn;i++){
		v=p[i];
		p[i]=v>vb?(v<v0?v0:(v>v1?v1:v)):vz;
	}
}
void countrect_(int*p,int *pm,int *pim,int*pjm,int *prect,int*pnc){
	int i,j;
	int im,jm,nc,m,i1,i2,j1,j2,i0,j0;
	im=*pim;jm=*pjm;m=*pm;nc=0;
	i1=i2=j1=j2=-1;
	for(j=1;j<=jm;j++){
		for(i=1;i<=im;i++,p++){
			if(*p!=m)continue;
			i0=i1=i2=i;j0=j1=j2=j;			
			nc++;
			break;
		}
		if(i1>0)break;
	}
	for(j=j0;j<=jm;j++){
		for(i=i0;i<=im;i++,p++){
			if(*p!=m)continue;
			nc++;
			i1=i1<=i?i1:i;
			i2=i2>=i?i2:i;
			j1=j1<=j?j1:j;
			j2=j2>=j?j2:j;
		}
		i0=1;
	}
	prect[0]=i1;	prect[1]=j1;	prect[2]=i2;	prect[3]=j2;
	*pnc=nc;
}
/*
	fortran interface :
	return int value at given pointer,
	for volatile variable,avoid optimize
*/
int IntVolatile(void*p){
	return *(int *)p;
}

//*====================================================================*//
#ifdef CP_SWF90
/*  in SW set cpu mode,for underflow
*/
union long_double {
        unsigned long long ldata;        double ddata;
};
void set_fpcr(unsigned long long a){
        union long_double fpcr;fpcr.ldata = a; asm("wfpcr %0"::"f"(fpcr.ddata));
}

unsigned long long get_fpcr(){
        union long_double fpcr; asm("rfpcr %0":"=f"(fpcr.ddata));return fpcr.ldata;
}
void syslimit_(long value){
        struct rlimit rlim;
        if (value == 0) rlim.rlim_cur = RLIM_INFINITY;
        else    rlim.rlim_cur = value;
        rlim.rlim_max = RLIM_INFINITY;
        if (setrlimit(RLIMIT_STACK, &rlim) == -1) {
                perror("setrlimit A:");
                exit(-1);
        }
        rlim.rlim_cur = RLIM_INFINITY;
        if (setrlimit(RLIMIT_CORE, &rlim) == -1) {
                perror("setrlimit B:");
                exit(-1);
        }
        return;
}

void set_rndmode_(){
        union long_double fpcr; fpcr.ldata = get_fpcr();fpcr.ldata |= (0x1ULL<<48)   ;set_fpcr(fpcr.ldata);
}
void reset_rndmode_(){
        union long_double fpcr; fpcr.ldata = get_fpcr();fpcr.ldata &= (~(0x1ULL<<48));set_fpcr(fpcr.ldata);
}

void initcpu_sw_(){  syslimit_(0);  set_rndmode_();}
void initcpu_(){  syslimit_(0);  set_rndmode_();}

void rpcc_(unsigned long *c)
{
  unsigned long a;
  asm("rtc %0": "=r" (a) : );
  *c=a;
};

#elif defined( INTEL)
void initcpu_(){  }
void rpcc_(unsigned long *c)
{
unsigned long a, d;
 __asm__ __volatile__ ("rdtsc" : "=a" (a), "=d" (d));
*c = a | (d << 32);
};

#else
void initcpu_(){  }
void rpcc_(unsigned long *c)
{
unsigned long a, d;
printf("unknown CPU not surport rpcc\not");
*c=0;
};
#endif
void meminf(int mpiid){
	FILE*fi;
	char fn[256];
	long rss;
	long pid=getpid();
	sprintf(fn,"/proc/%ld/status",pid);
	fi=fopen(fn,"rt");
	rss=0;
	while(!feof(fi)){
		char*ir=fgets(fn,256,fi);
		if(strncmp(fn,"VmRSS",5)==0){
			//printf("MEMINFA: %d %d %s AAA\n",mpiid,pid,fn);
			rss=atol(fn+6);break;
		}
		//if(strncmp(fn,"VMRSS")==0)rss=atol(fn+6);
	}
	fclose(fi);
	printf("MEMINFB: %d %ld %ld kb\n",mpiid,pid,rss);
}
void meminf_(int *mpiid){
	meminf(*mpiid);
}

#include <stdarg.h>
#define ferrrpt stderr
#define finfrpt stderr


#ifndef finfrpt
INMIC FILE* finfrpt=NULL;
static void vopen(FILE**f,char*fnt){
	char fn[128];
	if(*f)return;
	sprintf(fn,"%s_%d.txt",fnt,md.si.selfid);	
	*f=fopen(fn,"wt");
}
#endif
//#define MID md.si.selfid
#define MID 0
void errrpt(const char*s,int id,const char *fmt,...) {
	va_list args;
	int err;
	err=errno;
  va_start(args, fmt);
  fflush(stdout);
  fprintf(ferrrpt," %s %d %d",s,id,MID );
  vfprintf(ferrrpt,fmt, args);
  fprintf(ferrrpt," %d %s\n",err,strerror(err));
  va_end(args);	
  fflush(ferrrpt);
}

void infrpt(const char*s,int id,const char *fmt,...) {
	va_list args;
  va_start(args, fmt);
  fflush(stdout);
#ifndef finfrpt
  vopen(&ferrrpt,"Infrpt");
#endif
  fprintf(finfrpt," %s %d %d",s,id,MID);
  vfprintf(finfrpt,fmt, args);
  fprintf(finfrpt,"\n");
  va_end(args);	
  fflush(finfrpt);
}
#include <math.h>
int ccheckee(double *ee,int n){
  int i ,err=0;
  for(i=0;i<n;i++){
    if(isnan(ee[i])||isinf(ee[i])){
      printf("ee err:%d %lf \n",i,ee[i]);
      err=i+1;
      return err;
    }
  }
  return 0;
}
void ccheckee_(double *ee,int *n,int *err){
  *err=ccheckee(ee,*n);
}

#include <mpi.h>
void wav_abort(){
#ifdef NO_MPI
  exit(0);
#else
  MPI_Abort(MPI_COMM_WORLD,0);
#endif

}
