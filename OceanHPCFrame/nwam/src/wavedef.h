#ifndef WAVEDEF_H_INCLUDED
#define WAVEDEF_H_INCLUDED
#if 0 
  === Include File for Fortran and C Code
  === do not Define COMMENT__
  === comment Lines Must in #ifdef COMMENT__
  ==== #include "wavedef.h" =======
#  define TEST_PART
#  define CP_SWF90
#  define NO_MPI
#endif
#ifdef MMTHREAD
#define MTHREAD
#endif
#ifdef C_CALCULATE
#define C_PROPAGAT 
#define C_IMPLSCH 
#endif
#ifndef LOGSCURR
# define LOGSCURR 0
#else
# define LOGSCURR 1
#endif
#ifdef USEHBM
#define USE_C_ALLOC 
#endif

# define MBVDEP   50
#define SIGMALVL 0
#define CKL 32
#define CJNTHET 12
#define CMKJ (CJNTHET*CKL)
#define NKP 8
#define KLPAT NKP
#define NTSPLIT 10
#define MAXAVHIST 16
#ifdef _LFORTRAN 
! fortran code
!#define NODBGINFO
# define REALD real(8)
#  define DBGO0(lvl,format)  if(mpi_id==0.and.lvl<=DBGLvl)write(*,format)
#  define DBGO(lvl,format)  if(lvl<=DBGLvl)write(*,format)
# ifdef NODBGINFO
#  define OUT_ALLOCSIZE(v)   !
#  define DBGINFN(n) !
#  define DBGINF0  !
#  define DBGINF   !
#  define DBGINFA  !
#  define SDBGINF  !
#  define FFO      !
#  define FFO0     !
#  define DBG(lvl) !
#  define DBARRIER
# else
#  define OUT_ALLOCSIZE(v) !write(6,'(a,a8,a,i8,a,18i6)')"ALLOCATE:size(","v",")=",size(v),__FILE__,__LINE__,mpi_id;
#  define DBGINFN(n)  if(mpi_id==n)write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id,runstate%nTimeStep
#  define DBGINF0  if(mpi_id==0.and.DbgLvl>8)write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id,runstate%nTimeStep
#  define DBGINFA0 if(mpi_id==0.and.DbgLvl>8)write(*,*)__FILE__,__LINE__,mpi_id,runstate%nTimeStep,";"
#  define DBGINF   if(DbgLvl>8)write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id,runstate%nTimeStep
#  define DBGINFA  if(DbgLvl>8)write(*,*)__FILE__,__LINE__,mpi_id,runstate%nTimeStep,";"
#  define SDBGINF  if(DbgLvl>8)write(*,'(a,3i6,";",a,18i8)')__FILE__,__LINE__,mpi_id,runstate%nTimeStep
#  define FFO0     if(mpi_id==0)call flush6
#  define FFO      call flush6
#  define DBG(lvl) if(lvl<=DBGLvl)
#  define DBARRIER call wav_mpi_barrier
# endif
#else 
// C code
#define DBGINF printf("%s %d\n",__FILE__,__LINE__)
#endif
#endif
