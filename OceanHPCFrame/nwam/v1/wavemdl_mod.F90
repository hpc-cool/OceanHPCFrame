#include "wavedef.h"
#ifdef DEBUG
!#  define DBGINF
#  define CKEE(e)
#else
!#  define DBGINF   write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id
#  define CKEE(e)
#endif
# define DBGINF0   !if(mpi_id==0)write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id
Module WaveMdl_mod
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
    integer n,dnwp,cnwp,nnwp,exg_state,nwpb,ncheck
    call InitWaveMdl   ! wavemdl Init
    DBGO0(4,'(a)')trim(startType)
    call Init_Startup
    call Resettimer
    call starttimer(1)
    curexgid=1;
    ! Main Loop
    DO   !{
      call starttimer(3)
      runstate%nTimeStep=runstate%nTimeStep+1
      !set Model time ,not
      call starttimer( 7);call SetModelTime;call endtimer( 7)
      curexgid=2-(curexgid-1) ! 1=>2 2=>1
      ! exchange data in MPI
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
      !eet => eec
      call starttimer( 6);call implschs(eec,eet,1,nwpc)          ;call endtimer(6)
      !call starttimer( 7);call Monitor(eet,-1,0);call endtimer( 7)
      CKEE(eec)
      ! precalc not
      if(runstate%iPreCalT>0)then
        call DecIPrecalT
        if(runstate%iPreCalT<=0)then
          call Resettimer
          call starttimer(1)
        endif
        cycle
      endif
      !for output ,not
      call starttimer( 8);call ACCUMEA(eec,1,nwpc) ;call endtimer( 8)
      !get next wind ,not
      call starttimer( 9);key=IDataIO(4);call endtimer( 9)
      !call c_setwind
      if(runstate%hist_eot>0)then
        !call starttimer(10);call CheckOutPut(eec,1,nwpc,runstate%hist_eot);call endtimer(10)
      endif
      if(runstate%rest_eot/=0)then
        !call starttimer(11);call CheckRestart;call endtimer(11)
      endif
      call endtimer(3)
      if(runstate%stop_now/=0)exit
    end do !}
    call EndWaveMdl
  end SUBROUTINE RunWaveMdl !}

end Module WaveMdl_mod
