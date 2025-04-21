!#define USE_C_ALLOC
#include "wavedef.h"
module platform_init_mod
#ifndef NO_MPI
  use wav_mpi_mod
#endif
  use varcommon_mod
  use partition_mod
  use windin_mod
  use output_cal_mod
  use output_mod
  use restart_mod
  use propagat_mod
  use boundary_mod
  use implsch_mod
  !use,intrinsic::iso_c_binding
  !use wavemdl_mod
  IMPLICIT NONE

  public::Init_Distinct,InitAfterPartition,InitAfterMpi
  type FC_PARA
    integer(4) kl_,jnthet_,kld_,NBVDEP_,sigmalvl_,logscurr_
    integer(4) nwps,nwpc,nwpa;
    integer LXB,LYB,LXN,LYN
    integer(4) mpi_id,mpi_comm_wav,mpi_npe;
    real   (4) :: constwindx,constwindy
  end type FC_PARA
  real(8),allocatable::nfpcs(:)
  real(8),external::c_getfpcd
contains

  SUBROUTINE   InitAfterMpi
    call Init_cpu(mpi_id,mpi_npe,mpi_comm_wav, NCorePClu ,NThPClu ,NGrpPNode,NProcPNode ,ManageCoreId);
  end SUBROUTINE   InitAfterMpi
  SUBROUTINE End_Distinct
    !call C_ctEnd
  end SUBROUTINE End_Distinct
  real(8) function GetNfpc(nt)
    integer i,nt
  end function GetNfpc

#undef DBARRIER 
#undef FFO      
#undef DBGINF0  
#  define DBARRIER call wav_mpi_barrier
#  define FFO      call flush6
#  define DBGINF0  if(mpi_id==0.and.DbgLvl>8)write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id,runstate%nTimeStep
#define DBGZ DBARRIER;DBGINF0;FFO
  SUBROUTINE InitAfterPartition
    type(FC_PARA) gd
    integer sigmlvl
    integer res;
    !call c_checkstructs ;
    nullify(eec,eet,ip,pp_vs,ipos12,ipos8);
    sigmlvl=0
    !if(vdep[1]<0)  sigmlvl=1
    gd%kl_        =kl
    gd%jnthet_    =jnthet
    gd%kld_       =kld
    gd%NBVDEP_    =ndep
    gd%sigmalvl_  =sigmlvl
    gd%logscurr_  =logscurr
    gd%mpi_id     =mpi_id
    gd%mpi_npe    =mpi_npe
    gd%mpi_comm_wav=mpi_comm_wav
    gd%nwps       =nwps
    gd%nwpc       =nwpc
    gd%nwpa       =nwpa
    gd%lxb       =lxb
    gd%lyb       =lyb
    gd%lxn       =lxn
    gd%lyn       =lyn
    gd%constwindx =constwindx
    gd%constwindy =constwindy
    !call C_CheckcGrid(gd)

    call C_SetGPar(gd,iepos,ieind,nsp)
#ifdef USE_C_ALLOC
    call C_Allocate_Vee(res)! Set Gvar same time
    if(res<0)then
      print*,'C_Allocate_Vee allocate mem error,',mpi_id
      call InitMpi(-1)      ! Mpi Deinit
      stop
    endif
#endif
    if(.not.associated(eet))  ALLOCATE(eet(kl,jnthet,0:nwpa))
    if(.not.associated(eec))  ALLOCATE(eec(kl,jnthet,0:nwpa));
    eec=0;eet=0;
    allocate(nfpcs(mpi_npe*2))
    nfpcs=0
  end SUBROUTINE InitAfterPartition
  SUBROUTINE Init_Distinct
    call c_setpointers(eec,eet,ip,pp_vs,ipos12,ipos8,wxy)
    call c_Init_Distinct
    !call C_ctStart
  end SUBROUTINE Init_Distinct
  SUBROUTINE Init_Startup
    integer ir
    if(startType=='STARTUP')then
      ir=IDataIO(2)          !Get First Data for cold start
      call setspec(eec,1,nwpc,1)     ! cpp Not surport setspec
      call InitCheckOutPut(1)
    else
      call read_restart(eec)
      runstate%iPreCalT=0
      ir=IDataIO(3)        !Get First Data for warm start
      !call setspec(eec,1,nwpc,2)    ! now surport global only
    end if 
  end SUBROUTINE Init_Startup
end module platform_init_mod
!===============================
