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
# define DBGINF0   !if(mpi_id==0)write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id
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
    REALD,pointer:: ett(:,:,:)
    integer n,dnwp,cnwp,nnwp,exg_state,nwpb,ncheck
    !program Init
    call InitWaveMdl   ! wavemdl Init
    DBGO0(4,'(a)')trim(startType)
    call Init_Startup
    !call wav_mpi_barrier("CC");
    call Resettimer
    call starttimer(1)
    call c_setwind
    curexgid=1;
    ! Main Loop
    DO   !{
      runstate%nTimeStep=runstate%nTimeStep+1
      DBGINF0
      call SetModelTime
      call starttimer(3)
      curexgid=2-(curexgid-1) ! 1=>2 2=>1
      !call C_setstate(2,runstate%nTimeStep); ! if comment this line ,cal &mpi not orlay
      CKEE(eec)
      DBGINF0
#ifndef NO_MPI
      !call cpyeec(0);
      call starttimer( 4);call exchange_boundary_start(eec,curexgid);call endtimer(4)
      call starttimer( 4);call exchange_boundary_end()       ;call endtimer(4)
      !call cpyeec(1);
#endif
      DBGINF0
      !eec => eet
      CKEE(eec)
      call starttimer( 5);call propagats(eet,eec,1,nwpc)     ;call endtimer(5)
      CKEE(eet)
      DBGINF0
      !eet => eec
      call starttimer( 6);call implschs(eec,eet,1,nwpc)          ;call endtimer(6)
      DBGINF0
      CKEE(eec)
      call starttimer( 7);call setspec(eec,1,nwpc,2)  ;call endtimer( 7)      !set ee !      set water boundary
      DBGINF0
      CKEE(eec)
      !call starttimer( 7);call Monitor(eec,-1,0);call endtimer( 7)
      call starttimer( 8);call ACCUMEA(eec,1,nwpc) ;call endtimer( 8)
      call starttimer( 9);key=IDataIO(4);call endtimer( 9)
      DBGINF0
      call c_setwind
      DBGINF0
      if(runstate%iPreCalT>0)then
        call DecIPrecalT
        if(runstate%iPreCalT<=0)then
          call Resettimer
          call starttimer(1)
        endif
        cycle
      endif
      DBGINF0

      if(runstate%hist_eot>0)then
        !call starttimer(10);call CheckOutPut(eec,1,nwpc,runstate%hist_eot);call endtimer(10)
      endif
      DBGINF0
      if(runstate%rest_eot/=0)then
        !call starttimer(11);call CheckRestart;call endtimer(11)
      endif
      if(runstate%stop_now/=0)exit
      call endtimer(3)
    end do !}
    !DBGO0(0,*) "end wave model A",runstate%nTimeStep
    call flush6
    call EndWaveMdl
  end SUBROUTINE RunWaveMdl !}
end Module wavemdl_mod
