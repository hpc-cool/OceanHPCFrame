#include "wavedef.h"
#ifdef DBGINF
#undef DBGINF
#endif
#ifdef DEBUG
# define CKEE(et) call checkee(et,nwpa+1,__FILE__,__LINE__)
#ifdef NO_MPI
# define DBGINF   write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id
#else
# define DBGINF   call wav_mpi_barrier("AA");write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id
#endif
#else
# define CKEE(et) !call checkee(et,nwpa+1,__FILE__,__LINE__)
# define DBGINF
#endif
Module wavemdl_mod
  use varcommon_mod
  use implsch_mod
  use propagat_mod
  use boundary_mod
  use windin_mod
  use output_mod
  use output_cal_mod
  use partition_mod
  use restart_mod
  use waveinit_mod
  use platform_init_mod
  IMPLICIT NONE
  private
  public::RunWaveMdl
  integer ::curexgid=0
contains
  ! Main program
  SUBROUTINE RunWaveMdl !{
    integer key,ith,ko,iid,srcid,i
    real(8),allocatable:: irtcfs(:)
    real(8) rttt
    integer RFB
    REALD,pointer:: ett(:,:,:)
    integer n,dnwp,cnwp,nnwp,exg_state,nwpb,ncheck
    !program Init
    call InitWaveMdl   ! wavemdl Init
    DBGO0(4,'(a)')trim(startType)
    call Init_Startup
    call Resettimer
    call starttimer(1)
#ifdef MTHREAD
#define MWGR call waitstater(RFB,runstate%nTimeStep)
#define MWG  call waitstate(RFB,runstate%nTimeStep)
#define MSG  call setstate (RFB,runstate%nTimeStep)
#else
#define MWGR 
#define MWG  
#define MSG  
#endif
#ifdef MTHREAD
    RFB=1; if(startType/='STARTUP') RFB=2;
#endif
    call c_setwind
    MSG; MWG;
    curexgid=1;
    !wait calend
    RFB=10; MSG;MWG;
    DO   !{
      if(RFB>100000)then
        RFB=10; MSG;MWGR;
      endif
      runstate%nTimeStep=runstate%nTimeStep+1
      !print*,runstate%nTimeStep,Gettimer(3)
      call SetModelTime
      call starttimer(3)
      curexgid=2-(curexgid-1) ! 1=>2 2=>1
      CKEE(eec)
#ifdef MTHREAD
      !set wind ok
      call setflag(runstate%hist_eot,runstate%stop_now,runstate%rest_eot)
      !wait calend
      call c_setwind
      call starttimer( 9);key=IDataIO(4);call endtimer( 9)
      RFB=RFB+1; MSG;MWG; ! 1
#endif
#ifndef NO_MPI
      call starttimer( 4);call exchange_boundary_start(eec,curexgid);call endtimer(4)
      call starttimer( 4);call exchange_boundary_end()       ;call endtimer(4)
#endif
#ifdef MTHREAD
      !set boundary ok
      RFB=RFB+1; MSG; ! 2
#else
      if(1)then
        !eec => eet
        call starttimer( 5);call propagats(eet,eec,1,nwpc)     ;call endtimer(5)
        !eec <= eet
        call starttimer( 6);call implschs(eec,eet,1,nwpc)          ;call endtimer(6)
        call starttimer( 7);call setspec(eec,1,nwpc,2)  ;call endtimer( 7)      !set ee !      set water boundary
        !call starttimer( 7);call Monitor(eec,-1,0);call endtimer( 7)
        call starttimer( 8);call ACCUMEA(eec,1,nwpc) ;call endtimer( 8)
        call starttimer( 9);key=IDataIO(4);call endtimer( 9)
        call c_setwind
      endif
#endif
      if(runstate%iPreCalT>0)then
        call DecIPrecalT
        if(runstate%iPreCalT<=0)then
          call Resettimer
          call starttimer(1)
        endif
        cycle
      endif
#ifdef MTHREAD
      if(runstate%hist_eot>0)then ! check output have wait
        !call starttimer(10);call CheckOutPut(0,-1,RFB,eec,1,0,1,0,runstate%hist_eot);call endtimer(10)
      endif
      if(runstate%rest_eot/=0)then
        !call starttimer(11);call CheckRestart;call endtimer(11)
      endif
#else
      if(runstate%hist_eot>0)then
        !call starttimer(10);call CheckOutPut(0,-1,RFB,eec,1,nwpc,1,0,runstate%hist_eot);call endtimer(10)
      endif
      if(runstate%rest_eot/=0)then
        !call starttimer(11);call CheckRestart;call endtimer(11)
      endif
#endif
      if(runstate%stop_now/=0)exit
      call endtimer(3)
#ifdef MTHREAD
      !set IO OK
      !call starttimer( 9);call setstate(35,runstate%nTimeStep);call endtimer( 9)
#endif
    end do !}
    !DBGO0(0,*) "end wave model A",runstate%nTimeStep
    call flush6
#ifdef MTHREAD
    RFB=RFB+1000-10
    MSG;MWG;
#endif
    call EndWaveMdl
  end SUBROUTINE RunWaveMdl !}
end Module wavemdl_mod
