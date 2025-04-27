#include "wavedef.h"
!#define ALLLOG
!# define DBGINF   !call wav_mpi_barrier("AA");write(*,'(a,i6,i4.4,";",18i8)')__FILE__,__LINE__,mpi_id

Module waveinit_mod
#ifndef NO_MPI
  use wav_mpi_mod
#endif
  use varcommon_mod
  use implsch_mod
  use propagat_mod
  use boundary_mod
  use windin_mod
  use output_mod
  use partition_mod
  use restart_mod
  use platform_init_mod
IMPLICIT NONE
private
  integer::hist,rest
  public::InitWaveMdl,EndWaveMdl,DecIPrecalT,SetModelTime
  integer ::ndstep,oldday=-1
  real(8)  dfpc
contains
#define RINFO(TAG) if(mpi_id==0)write(*,'(i8," ",a,f8.3," ",i8)')iwalltime(),TAG,Difftimer(2),mpi_npe;  call flush(6);
#define DBGA call wav_mpi_barrier;write(*,'("DBGA",i4,":",i4.3)')__LINE__,mpi_id;  call flush(6);
#define DBGB write(*,'("DBGA",i4,":",4i4.3)')__LINE__,mpi_id,mpi_npe;  call flush(6);
  SUBROUTINE InitWaveMdl
    integer key,ir,out
    real*8 tt
    !ir=ieee_flags('clear', 'exception', 'all',out )
    call Resettimer
    call starttimer(1)
    call starttimer(2)
    if(IDataIO(0)<0)then   ! IGetWind(0) Init MPI msg etc,must before partition
      !DBGO(0,*)"Init MPI msg  Error"
      stop
    endif
    call InitMpi(0)      ! Mpi Init
    call zf_setpid(mpi_id) ! Set Mpi_id for zfile.c
    RINFO('InitMpi End')
#ifndef NO_MPI
    !call wav_mpi_gather(nodeid,NodeIds)
    !RINFO('Gather nodeid End')
    call gathercheck(__LINE__)
    RINFO('gathercheck End')
#endif
    !call testtimes
    !call zftest;         ! call zftest;         ! stop
    Call InitPara      !varcommon_mod Init
    call InitAfterMpi
    DBGLVL=20
    ! call wav_chStdOut
    if(mpi_id==0) then
      RINFO('InitPara End')
    else
      if(DBGLvl<5)DBGLvl=0
    endif
    call setwave
    RINFO('setwave End')
    !call nlweight;    stop
    RINFO('ALLOCATE depg,nspg begin')
    ALLOCATE(depg(gixl,giyl),nspg(gixl,giyl));
    RINFO('ALLOCATE depg,nspg End')
    call ReadTopog
    RINFO('ReadTopog End')
#ifndef NO_MPI
    call gathercheck(__LINE__)
    RINFO('gathercheck End')
#endif
    !partition !!!!!
    call partition
    RINFO('partition End')
    !call ReadTopog;    RINFO('ReadTopog End')
    !call InitMonitor
    !RINFO('InitMonitor End')
    if(mpi_id==0)then
      DBGO(0,*) 'KL       = ',kl
      DBGO(0,*) 'DBGLVL   = ',DBGLVL
      DBGO(0,*) 'DELTTS   = ',DELTTS
      DBGO(0,*) 'GXLON0   = ',GXLON0
      DBGO(0,*) 'GYLAT0   = ',GYLAT0
      DBGO(0,*) 'GXLON1   = ',GXLON1
      DBGO(0,*) 'GYLAT1   = ',GYLAT1
      DBGO(0,*) 'GRDSZX   = ',GRDSZX
      DBGO(0,*) 'GRDSZY   = ',GRDSZY
      DBGO(0,*) 'GCIRCLE  = ',GCIRCLE
      DBGO(0,*) 'DEPTHMOD = ',DEPTHMOD
      DBGO(0,*) 'LOGSCURR = ',LOGSCURR
      DBGO(0,*) 'LOGSDISS = ',LOGSDISS
      DBGO(0,*) 'nwpa = ',nwpa,nwpc,nwps,gnwpc,(gnwpc+0.)/mpi_npe,gixl,giyl
      DBGO(4,'(10i9)')nwpcs
      call flush(6);
    endif
    ! sum case Needed allocate Mem use C Program for mem align
    ! if Not ,InitAfterPartition  is null
    !print*,__FILE__,__LINE__,mpi_id
    !call wav_mpi_barrier();
    RINFO('InitAfterPartition start')
    call flush6

    call InitAfterPartition
    RINFO('InitAfterPartition')
    call set_neighbor
    RINFO('set_neighbor')
    call InitOutPut(0)   !Init OutPut
    RINFO('InitOutPut End')
    key=IDataIO(1)            !must Init before InitImplsch
    call InitImplsch
    RINFO('InitImplsch End')
    call InitPropagat
    RINFO('InitPropagat End')
    call InitSetspec
    RINFO('InitSetspec End')
    call endtimer(1)
    !    call wav_mpi_barrier();
    ! For every case ,sum Distinct Init here
    call Init_Distinct
    ! Init run
    runstate%nTimeStep=-1
    runstate%edaycur=dayinit
    runstate%timecur=timeinit
    runstate%ndays=0
    call InitCheckOutPut(0)
    !if(iCheckRestart>0)then
    !  call InitCheckOutPut(1)
    !  runstate%edaycur=dayinit
    !  runstate%timecur=timeinit
    !  runstate%ndays=0
    !  call setmodeltime
    !  call TestRestart(iCheckRestart,eec)
    !  RINFO('TestRestart End')
    !  stop
    !endif
    if(mpi_id==0)then
      RINFO( "wave model Init End")
    endif

    eec=0;eet=0;

    call NextTime(resttime,rest_option,rest_N,0)
    ndstep=3600/DELTTS
    if(ndstep<=1)ndstep=1
  end SUBROUTINE InitWaveMdl !}

  SUBROUTINE DecIPrecalT
    integer nn
    if( mod(runstate%iPreCalT,ndstep)==0)then
      !DBGO0(0,'("p",$)')
    endif
    nn=3600*4/DELTTS;if(nn<1)nn=1
    runstate%iPreCalT=runstate%iPreCalT-1
    if( mod(runstate%iPreCalT,nn)==0)then
       DBGO0(1,'(a,f10.2,a,i8,a,i8,a)')'runstate%iPreCalT use',Gettimer(1),"s",runstate%iPreCalT,' leave ',int(runstate%iPreCalT*DELTTS/3600),'hours'
       !call Monitor(eec,-1,0)

    endif
  end SUBROUTINE DecIPrecalT

  SUBROUTINE EndWaveMdl
    integer ir
    call End_Distinct
    call InitOutPut(-1)  ! Close  OutPut
    ir=IDataIO(-1)
    call InitMpi(-1)      ! Mpi Deinit
    DBGO0(0,*) "end wave model A",runstate%nTimeStep
    call flush6
  end SUBROUTINE EndWaveMdl

  subroutine SetModelTime
    real*8 dt,adt,tt,nfpca,ta
    integer  eday
    if(runstate%timecur>1-delltday*0.01)then
      runstate%timecur=runstate%timecur-1
      runstate%edaycur=runstate%edaycur+1
    endif
    call day2cdatetime(runstate%edaycur+runstate%timecur,runstate%cdatecur,runstate%ctimecur)
      !nfpca=GetNfpc(runstate%nTimeStep)*1e-9
      !if(mpi_id==0)then
      ! ta=Gettimer(3)+1e-300
      !  DBGO(1,'(i8.8,i7.6,i5,2f12.5,f15.5,":",f17.4,"G",f8.2,"% |T:",i,7f10.5)')runstate%cdatecur,runstate%ctimecur,runstate%ndays-1 ,dt,adt, &
      !  ta,nfpca,nfpca*100/(ta*dfpc),iwalltime()  !,Gettimer(4),Gettimer(5),Gettimer(6),Gettimer(7),Gettimer(8),Gettimer(9),Gettimer(10)
      !endif
     ! DBGINF,noOutPut,runstate%iPreCalT,runstate%hist_eot
    if(noOutPut==0 .and.runstate%iPreCalT<=0)then  !
      call CheckNeedOutPut(runstate%hist_eot)
      !call CheckNeedRestart(runstate%rest_eot)
    endif
    if(oldday/=runstate%edaycur)then
      runstate%ndays=runstate%ndays+1
      oldday=runstate%edaycur
      dt =Difftimer(1)
      call endtimer(1);
      if(runstate%ndays>1)then
        adt=Gettimer(1)/(runstate%ndays-1)
      else
        if(runstate%ndays==1)call Resettimer(1)
        adt=dt;
      endif
      if(runstate%ndays>RunDays)runstate%stop_now=1
      !print*,runstate%edaycur,runstate%ndays,oldday,runstate%nTimeStep
#ifdef ALLLOG
        write(1000+mpi_id,'(a)')" pid    date   time Days                  dtime    adtime  alltime exchange  propgat implschs    other   accuma    Input   output"
        write(1000+mpi_id,'(i4,i9.8," ",i6,i5,3i5,2f9.5,14f9.2)')mpi_id,runstate%cdatecur,runstate%ctimecur,runstate%ndays-1 ,nwps,nwpa,nwpc,dt,adt, &
        Gettimer(3),Gettimer(4),Gettimer(5),Gettimer(6),Gettimer(7),Gettimer(8),Gettimer(9),Gettimer(10)
        write(1000+mpi_id,'(14f9.3)')Gettimer(11),Gettimer(12),Gettimer(13),Gettimer(14),Gettimer(24),Gettimer(25),Gettimer(26),Gettimer(27),Gettimer(28)
#endif
      if(mpi_id==0)then

        if(runstate%ndays==0.or.mod(runstate%ndays,50)==1)then
        write(1,*)"   date   time Days       dtime      adtime     alltime" !  exchange   propgat   implsch     other    accuma     Input    output"
        write(*,*)"   date   time Days       dtime      adtime     alltime" !  exchange   propgat   implsch     other    accuma     Input    output"
        endif
        write(1,'(i8.8,i7.6,i5,2f12.5,f12.5,":","T:",i8,7f10.5)')runstate%cdatecur,runstate%ctimecur,runstate%ndays-1 ,dt,adt, &
        ta,iwalltime()  !,Gettimer(4),Gettimer(5),Gettimer(6),Gettimer(7),Gettimer(8),Gettimer(9),Gettimer(10)
        write(*,'(i8.8,i7.6,i5,2f12.5,f12.5,":","T:",i8,7f10.5)')runstate%cdatecur,runstate%ctimecur,runstate%ndays-1 ,dt,adt, &
        ta,iwalltime()  !,Gettimer(4),Gettimer(5),Gettimer(6),Gettimer(7),Gettimer(8),Gettimer(9),Gettimer(10)

        !Gettimer(3),0. !,Gettimer(4),Gettimer(5),Gettimer(6),Gettimer(7),Gettimer(8),Gettimer(9),Gettimer(10)
        !DBGO(1,'(14f9.3)')Gettimer(11),Gettimer(12),Gettimer(13),Gettimer(14),Gettimer(15),Gettimer(24),Gettimer(25),Gettimer(26),Gettimer(27),Gettimer(28)
        call flush6
      endif
      call  starttimer(1)
    endif
    if(0.and.mpi_id==0)then
      dt =Difftimer(1)
      write(*,'(i8.8,i7.6,i5,2f12.5,f12.5,":","S:",i8,7f10.5)')runstate%cdatecur,runstate%ctimecur,runstate%ndays-1 ,dt,adt, &
        ta,iwalltime()  !,Gettimer(4),Gettimer(5),Gettimer(6),Gettimer(7),Gettimer(8),Gettimer(9),Gettimer(10)
    endif
  END subroutine SetModelTime
end Module waveinit_mod

! called by C code for oprate Fortran variable

#ifdef USE_C_ALLOC
subroutine F_SetEVar(cec,cet)
  use varcommon_mod
IMPLICIT NONE
  REALD,target:: cec (kl,jnthet,0:nwpa)
  REALD,target:: cet (kl,jnthet,0:nwpa)
  eec=>cec;   eet=>cet;
  eec=0;eet=0;
end subroutine F_SetEVar
subroutine F_SetImpVar(cip,cpvws ,cwxy);
  use varcommon_mod
  use implsch_mod
  use boundary_mod
IMPLICIT NONE
  type(implsch_pack_type),target::cip(0:nwpc)
  real(8),target::cpvws(kl,0:nwpc)
  real(4),target::cwxy(4,0:nwpc)
  ip=>cip;
  pvws   =>cpvws
  wxy    =>cwxy;
  pvws=0;
  wxy=0.;
end subroutine F_SetImpVar
#define IMDPSP  2
#define IMDPSO  4
#define BVIMDPS 1
subroutine F_SetOutput(cpvkdo,cpebdep)
  use varcommon_mod
  use output_cal_mod
  IMPLICIT NONE
  real(8),target::cpvkdo(kld,IMDPSO,0:nwpc)
  real(8),target::cpebdep(kld,BVIMDPS+ndep,0:nwpc)
  pvkdo  =>cpvkdo
  pebdep =>cpebdep
  pvkdo=0;
  pebdep=0;

end subroutine F_SetOutput
subroutine F_SetProp(cpis,cpvkdp,cipos8,cipos12)
  use varcommon_mod
  use propagat_mod
  IMPLICIT NONE
  type(propagat_pre_type),target::cpis(kl*jnthet,0:nwpc)
  real(8),target::cpvkdp(kl,IMDPSP,0:nwpc)
  integer(4),target::cipos12(12,0:nwpc)
  type (propinf_type),target::cipos8(0:nwpc)

  pp_vs=>cpis;
  pvkdp  =>cpvkdp
  ipos12=>cipos12
  ipos8=>cipos8
  !pp_vs=0; ipos8=0;
  pvkdp=0; ipos12=0;
end subroutine F_SetProp
#endif
