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
#undef DBGINF
# define DBGINF   !write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id,runstate%nTimeStep,RFB,RFG;call flush(6)
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
#define MWG  call waitstate(RFB,runstate%nTimeStep)
#define MSG  call setstate (RFB,runstate%nTimeStep)
#define MWGR call waitstater(RFB,runstate%nTimeStep)

!#define DTIMES
#ifdef DTIMES
#define BT(i)  call tscb(i)
#define ET(i)  call tsce(i)
#define EBT(i) call tsceb(i);
#define PT()   call prtsc();
#else
#define BT(i)
#define ET(i)
#define EBT(i)
#define PT()
#endif
  SUBROUTINE RunWaveMdl !{
    integer key,ith,ko,iid,srcid,i
    real(8),allocatable:: irtcfs(:)
    real(8) rttt
    integer RFB,RFG
    REALD,pointer:: ett(:,:,:)
    integer n,dnwp,cnwp,nnwp,exg_state,nwpb,ncheck
    !program Init
    RFG=1;
    call InitWaveMdl   ! wavemdl Init
    DBGO0(4,'(a)')trim(startType)
    DBGINF
    call Init_Startup
    DBGINF
    call Resettimer
    call starttimer(1)
    !call sleep(1000);
#ifdef MTHREAD
    RFB=1; !if(startType/='STARTUP') RFB=2;
    DBGINF
    MSG; MWG;
    DBGINF
#endif
    call c_setwind
    !write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id,runstate%nTimeStep,RFB
    curexgid=1;
    ! Main Loop
    !wait threads ready
    ! W 1   25 L  r4 25
    ! S 1/2 25 L  15 25
    !TW 1   25 L r15 25
    !TS 1   25 L   4 25
    DBGINF
    RFB=10; MSG;MWG;
    DBGINF
    RFG=RFG+1;
    DBGINF
    call c_checkboundary(RFG) !MWG
    DBGINF
#ifndef NO_MPI
    call exchange_boundary_start(eec,curexgid);
    DBGINF
    call exchange_boundary_end()       ;
#endif
    DBGINF
    call c_sendboundary(RFG)
    DBGINF
    DBGO0(0,*)'model Begin ',Difftimer(2),mpi_npe
    call Resettimer
    call SetModelTime
    DO   !{
      BT(0);
      if(RFB>100000)then
        RFB=10; MSG;MWGR;
      endif
      DBGINF
      runstate%nTimeStep=runstate%nTimeStep+1
      call SetModelTime
      call starttimer(3)
      curexgid=2-(curexgid-1) ! 1=>2 2=>1
      CKEE(eec)
      DBGINF
#ifdef MTHREAD
      !set wind ok
      call setflag(runstate%hist_eot,runstate%stop_now,runstate%rest_eot)
      RFB=RFB+1; MSG; ! 1
      DBGINF
      ! nexttime wind
      key=IDataIO(4);
#endif
      DBGINF
      EBT(1);RFG=RFG+1;call c_checkboundary(RFG) !MWG//fixme:RFG
      DBGINF
#ifndef NO_MPI
      EBT(2);call exchange_boundary_start(eec,curexgid);
      call exchange_boundary_end()       
#endif
      DBGINF
#ifdef MTHREAD
      !set boundary ok
      EBT(3);call c_sendboundary(RFG) ; !MSB//fixme:RFG
      EBT(4);
      RFB=RFB+1;  !2 
      EBT(5);RFB=RFB+1; MSG; MWG ! 3 Next IDataIo must wait calend
#else
      if(1)then
        !eec => eet
        CKEE(eec)
        call propagats(eet,eec,1,nwpc)     
        CKEE(eet)
        !eec <= eet
        call implschs(eec,eet,1,nwpc)     
        CKEE(eec)
        call setspec(eec,1,nwpc,2)   !set ee !      set water boundary
        CKEE(eec)
        !call Monitor(eec,-1,0)
        call ACCUMEA(eec,1,nwpc) 
        key=IDataIO(4)
        call c_setwind
      endif
#endif
      DBGINF
      !wait calend
      EBT(6);call c_setwind

      DBGINF
      if(runstate%iPreCalT>0)then
        call DecIPrecalT
        if(runstate%iPreCalT<=0)then
          call Resettimer
          call starttimer(1)
        endif
        cycle
      endif
      EBT(7);

      DBGINF
#ifdef MTHREAD
      RFB=RFB+1;  ! 4
      if(runstate%hist_eot>0)then
        !MWG by CheckOutPut
        !call CheckOutPut(0,-1,RFB,eec,1,0,1,0,runstate%hist_eot);
      endif
      if(runstate%rest_eot/=0)then
        !call CheckRestart;
      endif
#else
      if(runstate%hist_eot>0)then
        !call CheckOutPut(0,-1,RFB,eec,1,nwpc,1,0,runstate%hist_eot);
      endif
      if(runstate%rest_eot/=0)then
        !call CheckRestart
      endif
#endif
      if(runstate%nTimeStep==0)then
        PT();
      endif
      if(runstate%stop_now/=0)exit
      DBGINF
      call endtimer(3)
      ET(7);
    end do !}
#ifdef MTHREAD
    DBGINF
    RFB=RFB+1000;MSG; !MWG;
    DBGINF
#endif
    PT();
    DBGINF
    call flush6
    call EndWaveMdl
  end SUBROUTINE RunWaveMdl !}
end Module wavemdl_mod
