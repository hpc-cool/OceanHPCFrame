Module debughlp_mod
#ifndef NO_MPI
  use wav_mpi_mod
#endif
  IMPLICIT NONE
public
  integer ::OutDType,NTYPE
  integer::dbglvl=10
#ifdef NO_MPI
   integer,public:: mpi_comm_wav=-1,mpi_id=0,mpi_npe=1
#endif
#define MAXTIMER 100
  real*8:: rtcb(0:MAXTIMER+1),rtcl(0:MAXTIMER+1)
!#define HAVE_CPU_TIME
#define HAVE_SYSTEM_CLOCK
contains
#  define DBGINF   write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__
#ifndef HAVE_CPU_TIME
# ifdef HAVE_SYSTEM_CLOCK
    subroutine cpu_time(rtime)
      real(8) rtime
      integer COUNT
      integer,save::last_count=-1,COUNT_RATE, COUNT_MAX
      real(8),save::dtime,rcr=1;
      if(last_count/=-1)then
        CALL SYSTEM_CLOCK(COUNT)
        if(COUNT<last_count)then
          dtime=dtime+dble(COUNT_MAX)/float(COUNT_RATE)
        endif
      else
        CALL SYSTEM_CLOCK(COUNT, COUNT_RATE, COUNT_MAX)
        dtime=0;rcr=1.D0/dble(COUNT_RATE)
        !rcr=0.001
      endif
      last_count=COUNT
      rtime=dtime +COUNT*rcr
    end subroutine cpu_time
# else
    subroutine cpu_time(rtime)
      real(8) rtime
      integer vals(8)
      integer,save::oday=-1
      real(8),save::dtime=0
      DBGINF
      call date_and_time(VALUES=vals)
      if(vals(3)/=oday)dtime=dtime+86400.
      oday=vals(3)
      rtime=dtime+vals(5)*3600+vals(6)*60+vals(7)+vals(8)*0.001
    end subroutine cpu_time
# endif
		
#endif
		integer function iwalltime() 
      integer vals(8)
      call date_and_time(VALUES=vals)
      iwalltime=vals(5)*10000+vals(6)*100+vals(7)
		end function iwalltime
  SUBROUTINE Resettimer(iid)
    integer,optional:: iid
    real(8) rtime
    rtime=0;
    call cpu_time(rtime)
    if(present( iid))then
      rtcl(iid)=0;rtcb(iid)=rtime
    else
      rtcl=0;rtcb=rtime
    endif
!    call C_SetTimerPos(rtcb,rtcl)
  end SUBROUTINE Resettimer

  SUBROUTINE endstarttimer(iide,iids)
    integer iide,iids
    real*8 rtcc
    call cpu_time(rtcc)
    if(iide>=1)rtcl(iide)=rtcl(iide)+(rtcc-rtcb(iide))
    if(iids>=1)rtcb(iids)=rtcc
  end SUBROUTINE endstarttimer

  SUBROUTINE starttimer(iid)
    integer iid
    call cpu_time(rtcb(iid))
  end SUBROUTINE starttimer

  SUBROUTINE endtimer(iid)
    integer iid
    real(8) rtime
    integer COUNT
    call cpu_time(rtime)
    rtcl(iid)=rtcl(iid)+(rtime-rtcb(iid))
  end SUBROUTINE endtimer

  real*8 function Difftimer(iid)
    integer iid
    real(8) rtime
    call cpu_time(rtime)
    Difftimer=(rtime-rtcb(iid)) !+rtcl(iid)
  end function Difftimer

  real*8 function Gettimer(iid)
    integer iid
    Gettimer=rtcl(iid)
  end function Gettimer
!==========================================================================
  subroutine wav_stdio()
#ifdef USE_SHR_MODS
      !call shr_msg_chdir   ('wav') ! changes cwd
      call shr_msg_stdio('wav')
#endif
  end subroutine wav_stdio
  SUBROUTINE wav_chStdOut
    character*256 pid
    call flush6
    if(mpi_id>=0)then
      close(6)
      if(mpi_id==0)then
        write(pid,'("wav.log.",i6.6)')mpi_npe
        open(unit=6,file=trim(pid),status='unknown')
        write(6,*)mpi_id,'Redirect Out to ',trim(pid)
      else
        if(DBGLvl<=5)then
          DBGLvl=0
        else
          ! write(pid,'("logs/log",i6.6,"/wav.log.",i6.6)')mpi_npe,mpi_id
          write(pid,'("wav.log.",i6.6)')mpi_id
          open(unit=6,file=trim(pid),status='unknown')
          write(6,*)mpi_id,'Redirect Out to ',trim(pid)
        endif
      endif
      call flush6
    endif
  END SUBROUTINE wav_chStdOut

  subroutine wav_abort(msg)
    character*(*) msg
#ifndef NO_MPI
    call wav_mpi_abort(msg,-1)
#else
    write(*,*)msg
    stop
#endif
  end subroutine wav_abort

  subroutine WaitPrev
    integer it
#ifndef NO_MPI
    if(mpi_id>0)CALL wav_mpi_RECV(it,mpi_id-1,100)
#endif
  end subroutine WaitPrev
  subroutine NEXTCPU
    integer it
#ifndef NO_MPI
    if(mpi_id<mpi_npe-1)call wav_mpi_send(it,mpi_id+1,100)
#endif
  end subroutine NEXTCPU
  subroutine Zbarrier
    integer it
    integer,allocatable,save::vi(:)
#ifndef NO_MPI
  if(.not.allocated(vi))allocate(vi(mpi_npe))
  call wav_mpi_gather(it,vi,0)
  call wav_mpi_scatter(vi,it,0)
!  deallocate(vi)
#endif
  end subroutine Zbarrier
  subroutine flush6
    call flush (6)
    !call flush (0)
  end subroutine flush6
  subroutine flusht(ift)
  integer ift
    call flush (ift)
  end subroutine flusht



end Module debughlp_mod
