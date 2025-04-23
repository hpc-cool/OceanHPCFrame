#include "wavedef.h"


module split_mod
  use partition_mod
  use wav_mpi_mod
  implicit none

  type segdef_type
    integer(2) :: id ! :=0,1,...,nBLK  ! 0 means global
    integer(2) :: ty ! 1:write to mpi for send 
    ! 2:write to mpi inner;
    ! 3:read from mpi from recv
    ! 5:mpi read from numa for mpi send ;
    ! 6:mpi read from numa for mpi inner;
    ! 7:mpi write to numa from mpi recv
    ! 8:write 2 neihbor ; 
    ! 9:read from neighbor
    integer(4) :: sl ! start in local.
    integer(4) :: sr ! start in remote.
    integer(4) :: np ! number of points.
  end type segdef_type

  type BLKpmap_type
    integer(4) :: nwps ! 1-->nwps: MPI send
    integer(4) :: nwpw ! nwps+1-->nwpw: BLK write
    integer(4) :: nwpc ! nwpw-->nwpc: pure inner points. End of calculate points.
    integer(4) :: nwpr ! nwpc-->nwpr: recv from other blocks.
    integer(4) :: nwpa ! all points.
    integer(4) :: nseg ! number of segs.
    integer(4) :: fill(2)
    type(segdef_type) :: segs(1000)
  end type BLKpmap_type
  integer(4),pointer::ipos12l(:,:)

  type dnb8_type
    integer(4) :: np
    integer,pointer :: nb8(:,:) ! (8,np)
  end type dnb8_type

  type dieind
    integer lxn,lyn,lxb,lyb
    integer,pointer::ieind(:,:) !(lxn,lyn)
    !integer,pointer::
  end type dieind

  contains

#define C_DBGINF call C_DBGL(0,__LINE__,mpi_id,0)
  subroutine subsplit
    real(4) numamips(100)
    integer nnuma
    type(dnb8_type),allocatable:: nb8s(:)
    type(BLKpmap_type),allocatable:: pmap(:)
    type(dieind),allocatable::ieinds(:)
    integer i;
    call C_GetnumaInfo(nnuma,numamips)
    !write(6,*)"MPI",mpi_id," NodeName=",trim(NodeName)," numaS=",nnuma 
    if(.not. associated(ipos12l))then
      ALLOCATE(ipos12l(12,0:nwpc))
    endif
    allocate(pmap(0:nnuma),nb8s(nnuma),ieinds(0:nnuma))
    call set_list_for_BLK(nnuma,numamips,nb8s,pmap)

    !DBGINF,pmap(:)%nwpc
    call C_SetSegs(nnuma,pmap)

    do i=1,nnuma
      !call C_SetSegNb8(i-1,nb8s(i)%np,nb8s(i)%nb8)
      !call C_SetIEInd(i-1,ieinds(i)%lxn,ieinds(i)%lyn,ieinds(i)%ieind)
      deallocate(nb8s(i)%nb8)
    enddo
    deallocate(nb8s,pmap)
  end subroutine  subsplit

  subroutine set_list_for_BLK(nBLK,BLKmips,nb8,pmap)
!  subroutine set_list_for_BLK(ieind,ipos8,BLKmips,nBLK,lxn,lyn,nwpc,nwpa,nwps)
    ! Input:
    integer :: nBLK       ! - Number of BLKs & main CPU (few BLK and one main CPU).
!    integer :: LXN        ! - Size of local matrix in i-direction.
!    integer :: LYN        ! - Size of local matrix in j-direction.
!    integer :: nwpc       ! - Number of calculation points: sending points & inner points
!    integer :: nwpa       ! - All points, including inner points and outer points.
!    integer :: nwps       ! - Sending of sending points.
!    integer :: ieind(:,:) ! - Global table of ij2s from IRRP (LXN,LYN), 0 is for land.
!-------------------------------------------------------------------------------------------------
! - Order of points
!      0,1,2,...,nwps,nwps+1,...,nwpc,nwpc+1,...,nwpa
! - Oreder of ipos8
!      4---3---2
!      5   0---1
!      6---7---8
!-------------------------------------------------------------------------------------------------
!    integer :: ipos8(:,:) ! - Neighbors of inner points.(0:8,0:nwpc)
    real(4) :: BLKmips(nBLK) ! - Speed of each BLKs or main CPU.(nBLK)
    ! Output:
    type(dnb8_type) :: nb8(nBLK)
    type(BLKpmap_type) :: pmap(0:nBLK) ! 0 is for MPI (global)
!-------------------------------------------------------------------------------------------------
!          m(W)        D       U          I          R     d       u        R     M
!      |-----------|-------|--------|------------|------|-------|-------|------|------|
!     01          nwps             nwpw         nwpc                          nwpr   nwpa
!     * Calculate PNTs order:        m(W) --->  D ---> U --->  I
!     * Outer PNTs needed:    0 --->  M   --->  u ---> d ---> mdw
!     * Local PNTs order:     0 ---> m(W) --->  D ----> U ---->  I  ----> u ------ d ------ M
!-------------------------------------------------------------------------------------------------
    integer,allocatable :: ij2sBLK(:)
    integer,allocatable :: ij2ptyp(:,:)
    integer,allocatable :: ptable(:,:),recs(:)
    integer,pointer :: ip8(:)
    integer :: i,j,s,iBLK,iBLK1,iBLK2,npts,ii,snp,s1,it
    real(4) smips,tmips
    integer :: nsl
    do iBLK=0,nBLK
    pmap(iBLK)%nwps=0
    pmap(iBLK)%nwpw=0
    pmap(iBLK)%nwpc=0
    pmap(iBLK)%nwpr=0
    pmap(iBLK)%nwpa=0
    pmap(iBLK)%nseg=0
    enddo
    pmap(0)%nwps=nwps
    pmap(0)%nwpw=nwps
    pmap(0)%nwpc=nwpc
    pmap(0)%nwpr=nwpc
    pmap(0)%nwpa=nwpa
!    if(mpi_id==0)then
!      open(11,file='ipos8')
!      do s=1,nwpc
!        write(11,'(10i8)')s,ipos8(:,s)
!      enddo
!      close(11)
!    endif
!#ifdef cc
!-------------------------------------------------------------------------------------------------
! bl,br,nn,id=(0,1,2,3...)
!
!  Local PNTs order:     0 ---> m(W) --->  D ---> U --->  I  ---> (R)d ---> u(R) ---> M
!                            ---nwps ------------nwpw----nwpc---------------nwpr-----mwpa
!
!                            |-----------|-------|-----|-------|----|------|----|---|------|
!
!                                m(w)        D      U      I      R    d     u    R    M
!                            |-----------|-------|-----|-------|----|------|----|---|------|
!                                       nwps          nwpw    nwpc                 nwpr----nwpa
!
!
!
!-------------------------------------------------------------------------------------------------
!      (1) Keep the order of MPI send & recv
!      (2) For the points read by other BLK, keep them continue, but can be aR-Ranged by some parts.
!      (3) For the points read from the other BLK, keep them continue, but can be aR-Ranged in some parts.
!      (4) Keep the same order in the other BLK.
!-------------------------------------------------------------------------------------------------
! - Symbles defination.
!      I    - Pure inner points;
!      M    - MPI recv points;
!      m    - MPI sending points;
!      d(R) - BLK reading points from down side.
!      u(R) - BLK reading points from up side.
!      D(W) - BLK writing points by down side.
!      U(W) - BLK writing points by up side.
!      R    - MPI sending points in the inner cornner. It will be read by neighbor.
!      r    - BLK reading points and also MPI sending points from other BLK parts.
! - Location of different point type.
!      M | R-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-R | M
!      --+-----------------------------------------------------------------------------------+--
!      M | W-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-U-W | M
!      M | m-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-m | M
!      M | m-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-m | M
!      M | m-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-m | M
!      M | m-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-m | M
!      M | m-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-m | M
!      M | m-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-m | M
!      M | m-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-m | M
!      M | m-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-m | M
!      M | m-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-m | M
!      M | m-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-m | M
!      M | m-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-m | M
!      M | m-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-m | M
!      M | m-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-m | M
!      M | m-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-m | M
!      M | m-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-m | M
!      M | m-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-m | M
!      M | m-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-I-m | M
!      M | W-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-D-W | M
!      --+-----------------------------------------------------------------------------------+--
!      M | R-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-u-R | M
!-------------------------------------------------------------------------------------------------
! - Calculating order:
!     (1) Calculate Inner Points (I)
!     (2) Prepare Reading Points (R) & (a,b,c,d), then calculate Write points (W),
!     (3) Prepare (M), then calculate Sending points (m), and conner points (A,B,C,D).
!     Calculate MPI sending points (mW), pure inner points (I), BLK writing points (UDW),
!     Calculate PNTs order:        m(W) --->  D ----> U ---->  I
!     Outer PNTs needed:    0 --->  M   --->  u ----> d ----> mdw
!     Local PNTs order:     0 ---> m(W) --->  D ----> U ---->  I  ----> u ------ d ------ M
!-------------------------------------------------------------------------------------------------
!! C Code Needed
!type segdef
!integer(4) id; !i=-1,0,1...nBLK-1; -1 mean MPI DATA
!integer(4) nl,nr,np;
!end type segdef
!type BLKrunPara_type
! !Not needed in C Code
!integer(4) nwps,nwpw,nwpc,nwpr,nwpa;
!integer(4) nseg ;
!type(segdef)::segs(100)
!end type BLKrunPara_type
!type BLKtemp_type
!integer(4) bwpc,ewpc;
!end type BLKtemp_type
!type (BLKrunPara_type) ::BLK (-1:100)
!type (BLKrunData_type) ::BLKd(-1:100)
!! segs Need Map relationships
!! calc nodes read from neighbor
!! calc nodes write to neighbor
!! calc nodes read from MPIData
!! calc nodes write to MPIData
!! MPI Node write to calc nodes
!! MPI Node readfrom to calc nodes
!! End C Code Needed
!#endif
!-------------------------------------------------------------------------------------------------
    allocate(ij2sBLK(0:nwpa))
    !----------------------------------------------------------------------------------
    ! --- Set for index of BLK for each calculate point.
    smips=sum(BLKmips)
    tmips=0;snp=0;
    do iBLK=1,nBLK
    tmips=tmips+BLKmips(iBLK)
    npts=nwpc*tmips/smips
    if(iBLK==nBLK)npts=nwpc
    pmap(iBLK)%nwpc=npts-snp;
    snp=npts
    enddo
    !write(*,*)'nwpc=',nwpc
    ! SetBLK
    !do iBLK=1,nBLK
    !  write(6,*)'nwpc=',iBLK,pmap(iBLK)%nwpc
    !enddo
    iBLK=1; npts=0; ij2sBLK=0;
    !DBGINF,LXN,LYN
    do j=1,LYN
    do i=1,LXN
    s=ieind(i,j);
    if(s<=0)cycle
    if(s<=nwpc)then
    ij2sBLK(s)=iBLK;
    npts=npts+1
  endif
    if(npts==pmap(iBLK)%nwpc.and.iBLK/=nBLK)then
    iBLK=iBLK+1; npts=0
  endif
  enddo
    enddo
    !DBGINF,LXN,LYN
    !----------------------------------------------------------------------------------
    ! --- Set point type.
    allocate(ij2ptyp(0:nwpa,nBLK));ij2ptyp=1000
    do s=1,nwpc
    iBLK=ij2sBLK(s)
    if(s<=nwps)then
    ij2ptyp(s,iBLK)=100     ! - For MPI send points (m)
  elseif(s<=nwpc)then
    ij2ptyp(s,iBLK)=300     ! - For MPI inner points (I)
  endif
    enddo
 !   DBGINF,LXN,LYN
    do s=1,nwpc
    ip8=>ipos8(:,s); iBLK=ij2sBLK(s)
!      write(*,*)lbound(ip8),ubound(ip8),s,ip8;      stop
!      write(*,*)__LINE__,lbound(ip8),ubound(ip8),s
!      write(*,*)__LINE__,lbound(ipos8(:,s)),ubound(ipos8(:,s)),s
!      write(*,*)__LINE__,ip8(1:5)
!      write(*,*)__LINE__,ipos8(1:5,s)
!      write(*,*)__LINE__,ipos8(0:5,s);      stop
!      *** Turn out: {
!         268           1           9           1
!         269           1           9           1
!         270           1        2601        2602        2472        2473
!         271        2601        2602        2472        2473           2
!         272           1        2601        2602        2472        2473
!           2
!      *** Turn out: }
    if(iBLK==0)cycle
    do ii=1,9
    s1=ip8(ii);if(s1==0)cycle;
    iBLK1=ij2sBLK(s1);
    if(s1>nwpc)then
    ij2ptyp(s1,iBLK )=500 ! - MPI recv pnts
    cycle
  endif
    if(iBLK1==iBLK .or. iBLK1==0)cycle
    if(iBLK1>iBLK)then
    if(s<nwps)then
    ij2ptyp(s ,iBLK )=211 ! - For WU-Point
  else
    ij2ptyp(s ,iBLK )=210 ! - For U-Point
  endif
    if(s1<nwps)then
    ij2ptyp(s1,iBLK )=411 ! - For Rd-Point
  else
    ij2ptyp(s1,iBLK )=410 ! - For d-Point
  endif
  endif
    if(iBLK1<iBLK)then
    if(s<nwps)then
    ij2ptyp(s ,iBLK )=221 ! - For WD-Point
  else
    ij2ptyp(s ,iBLK )=220 ! - For D-Point
  endif
    if(s1<nwps)then
    ij2ptyp(s1,iBLK )=421 ! - For Ru-Point
  else
    ij2ptyp(s1,iBLK )=420 ! - For u-Point
  endif
  endif
  enddo
    nullify(ip8)
    enddo
    !----------------------------------------------------------------------------------
    !count points
    do iBLK=1,nBLK
    pmap(iBLK)%nwps=0;pmap(iBLK)%nwpw=0;
    pmap(iBLK)%nwpc=0;pmap(iBLK)%nwpr=0;pmap(iBLK)%nwpa=0
    do s=1,nwpa
    it=ij2ptyp(s,iBLK);
    if(it==100 .or. it==211 .or. it==221)then
    pmap(iBLK)%nwps=pmap(iBLK)%nwps+1     ! --- nwps: 100m,211WUm,221WDm
  else if(it==220 .or. it==210)then
    pmap(iBLK)%nwpw=pmap(iBLK)%nwpw+1     ! --- nwpw: 220D, 210U
  else if(it==300)then
    pmap(iBLK)%nwpc=pmap(iBLK)%nwpc+1     ! --- nwpc: 300I
  else if(it==420 .or. it==410 .or. it==421 .or. it==411)then
    pmap(iBLK)%nwpr=pmap(iBLK)%nwpr+1     ! --- nwpr: 420u,421Ru,410d,411Rd
  else if(it==500)then
    pmap(iBLK)%nwpa=pmap(iBLK)%nwpa+1     ! --- nwpa: 500M
  endif
  enddo
    pmap(iBLK)%nwpw=pmap(iBLK)%nwpw+pmap(iBLK)%nwps
    pmap(iBLK)%nwpc=pmap(iBLK)%nwpc+pmap(iBLK)%nwpw
    pmap(iBLK)%nwpr=pmap(iBLK)%nwpr+pmap(iBLK)%nwpc
    pmap(iBLK)%nwpa=pmap(iBLK)%nwpa+pmap(iBLK)%nwpr
    !write(6,*)'CCC',iBLK,pmap(iBLK)%nwps,pmap(iBLK)%nwpw,pmap(iBLK)%nwpc,pmap(iBLK)%nwpr,pmap(iBLK)%nwpa
    enddo
    !DBGINF,pmap(:)%nwpc
    !----------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------!
    !          m(W)        D       U          I          R     d       u        R     M               !
    !      |-----------|-------|--------|------------|------|-------|-------|------|------|           !
    !     01          nwps             nwpw         nwpc                          nwpr   nwpa         !
    ! GLB 01          nwps                                                        nwpc   nwpa         !
    !     * Calculate PNTs order:        m(W) --->  D ---> U --->  I                                  !
    !     * Outer PNTs needed:    0 --->  M   --->  u ---> d ---> mdw                                 !
    !     * Local PNTs order:     0 ---> m(W) --->  D ----> U ---->  I  ----> u ----> d ----> M       !
    !                                 100m       220D    210U     300I     420u    410d    500M       !
    !                               211WUm                                421Ru   411Rd               !
    !                               221WDm                                                            !
    !-------------------------------------------------------------------------------------------------!
    allocate(ptable(0:nwpa,nBLK));ptable=0
    allocate(recs(nBLK));recs=0
    !if(mpi_id==0)then
    !	s=723
    !	do iBLK=1,nBLK
    !		write(*,'(a,20i6)')'AAAAAZ',__LINE__,iBLK,s,nwps,nwpc,ptable(s,iBLK),pmap(iBLK)%nwps,pmap(iBLK)%nwpc
    !  enddo
    !endif
    do s=1,nwps
    iBLK=ij2sBLK(s);it=ij2ptyp(s,iBLK);
    !if(it==100 .or. it==211 .or. it==221)then ! must be sure
    recs(iBLK)=recs(iBLK)+1;ptable(s,iBLK)=recs(iBLK)   ! --- Mpi Send Pnts
    !endif
    enddo
    do s=nwps+1,nwpc
    iBLK=ij2sBLK(s);it=ij2ptyp(s,iBLK);
    if(it==220 .or. it==210)then
    recs(iBLK)=recs(iBLK)+1;ptable(s,iBLK)=recs(iBLK)   ! --- write to neighbor points
  endif
    enddo
    do s=nwps+1,nwpc
    iBLK=ij2sBLK(s);it=ij2ptyp(s,iBLK);
    if(it==300)then
    recs(iBLK)=recs(iBLK)+1;ptable(s,iBLK)=recs(iBLK)   ! ---  inner points
  endif
    enddo
    do s=1,nwpc
    do iBLK=1,nBLK
    it=ij2ptyp(s,iBLK);
    if(it==420 .or. it==421 .or. it==410 .or. it==411)then
    recs(iBLK)=recs(iBLK)+1;ptable(s,iBLK)=recs(iBLK) ! --- readfrom neighbor points
  endif
  enddo
    enddo
    do s=nwpc+1,nwpa
    do iBLK=1,nBLK
    it=ij2ptyp(s,iBLK);
    if(it==500)then
    recs(iBLK)=recs(iBLK)+1;ptable(s,iBLK)=recs(iBLK) ! --- MPI recv points
  endif
  enddo
    enddo
    !if(mpi_id==0)then
    !  open(11,file='aaaa')
    !  do j=1,LYN
    !    do i=1,LXN
    !      s=ieind(i,j);
    !      if(s<=0)cycle
    !      write(11,'(16i10)')i,j,ij2sBLK(s),ij2ptyp(s,:),ptable(s,:),s
    !    enddo
    !  enddo
    !  close(11)
    !endif
    !----------------------------------------------------------------------------------
    ! nb8
    !  type dnb8_type
    !    integer(4) :: np
    !    integer,pointer :: nb8(:,:) ! (8,np)
    !  end type dnb8_type
    !  type(dnb8_type) :: nb8(nBLK)
    !----------------------------------------------------------------------------------
    do iBLK=1,nBLK
    nb8(iBLK)%np=pmap(iBLK)%nwpc
    allocate(nb8(iBLK)%nb8(8,0:pmap(iBLK)%nwpc));nb8(iBLK)%nb8=0;
    enddo
    do s=1,nwpc
    iBLK=ij2sBLK(s);i=ptable(s,iBLK)
    if(i==0)cycle ! yinxq 2015-4-5 9:32:15
    do j=1,8
    s1=ipos8(j,s);nb8(iBLK)%nb8(j,i)=ptable(s1,iBLK) !ipos8(j,s) ! in global seqs
  enddo
    enddo
    !----------------------------------------------------------------------------------
    ! segs
    pmap(0)%nseg=0;
    do iBLK=1,nBLK
    j=0;pmap(iBLK)%nseg=0
    do s=1,nwps                                  ! mpi send
    i=ptable(s,iBLK);iBLK2=ij2sBLK(s);
    if(iBLK==iBLK2.and.i>0.and.i<=pmap(iBLK)%nwps)then
    j=i+AddSeg(pmap,nBLK,iBLK,0,i/=j,1,5,i,s)
  else
    j=0
  endif
  enddo
    j=0
    do s=nwps+1,nwpc                             ! mpi inner
    i=ptable(s,iBLK)
    if(i>0.and.i<=pmap(iBLK)%nwpc)then
    j=i+AddSeg(pmap,nBLK,iBLK,0,i/=j,2,6,i,s)
  else
    j=0
  endif
  enddo
    j=0
    do s=nwpc+1,nwpa                             ! mpi recv
    i=ptable(s,iBLK);iBLK2=ij2sBLK(s);
    if(i>0)then
    j=i+AddSeg(pmap,nBLK,iBLK,0,i/=j,3,7,i,s)
  else
    j=0
  endif
  enddo
    enddo
    do iBLK=1,nBLK
    j=0;iBLK2=0;
    do s=1,nwpc
    i=ptable(s,iBLK)
    iBLK1=ij2sBLK(s)
    ii=ptable(s,iBLK1)
    if(iBLK1/=iBLK .and. i>0 .and. ii>0)then
    j=i+AddSeg(pmap,nBLK,iBLK,iBLK1,iBLK1/=iBLK2.or.i/=j,8,9,i,ii)
  else
    j=0
  endif
    iBLK2=iBLK1
  enddo
    enddo
    !DBGINF,nwpc,shape(ipos12),shape(ipos12l);
      !do s1=nwpc+1,nwpa
      !	DBGINF,s1,ptable(s1,1),ij2sBLK(s1),ij2ptyp(s1,1)
      !enddo
    do s=1,nwpc
      iBLK=ij2sBLK(s);
      do j=1,12
        s1=ipos12(j,s)
        ipos12l(j,s)=ptable(s1,iBLK);    		    	
        !if(mpi_id==0)then
        !	if(iBLK==1)then
        !  	DBGINF,s,j,iBLK,s1,ptable(s1,1),ij2sBLK(s1),ij2ptyp(s1,1)
        !  endif
        !endif
      enddo	    
    enddo
    !if(mpi_id==0)then
    !  write(*,'(a10,10i6)')'SDSDSD-',0,0,0,nwps,nwpc,nwpa
    !  do iBLK=0,nBLK
    !    write(*,'(a10,10i6)')'SDSDSD-',iBLK,0,pmap(iBLK)%nseg, &
    !             pmap(iBLK)%nwps,pmap(iBLK)%nwpw,pmap(iBLK)%nwpc,pmap(iBLK)%nwpr,pmap(iBLK)%nwpa
    !    do nsl=1,pmap(iBLK)%nseg
    !      write(*,'(a10,10i6)')'SDSDSD-',iBLK,1,pmap(iBLK)%segs(nsl)%id,&
    !                pmap(iBLK)%segs(nsl)%ty,&
    !                pmap(iBLK)%segs(nsl)%sl,&
    !                pmap(iBLK)%segs(nsl)%sr,&
    !                pmap(iBLK)%segs(nsl)%np
    !    enddo
    !  enddo
    !endif
  end subroutine set_list_for_BLK

  integer function AddSeg(pmap,nBLK,iBLK,iBLK1,new,tyl,tyr,il,ir)
    integer nBLK,iBLK,iBLK1,tyl,tyr,il,ir
    type(BLKpmap_type) :: pmap(0:nBLK) ! 0 is for MPI (global)
    logical new
    integer nsl,nsr
    nsl=pmap(iBLK)%nseg         ;  nsr=pmap(iBLK1)%nseg
    if(new)then
  nsl=nsl+1;                      nsr=nsr+1;
  pmap(iBLK)%nseg        =nsl  ;  pmap(iBLK1)%nseg        =nsr ;
  pmap(iBLK)%segs(nsl)%sl=il   ;  pmap(iBLK1)%segs(nsr)%sl=ir  ;
  pmap(iBLK)%segs(nsl)%sr=ir   ;  pmap(iBLK1)%segs(nsr)%sr=il  ;
  pmap(iBLK)%segs(nsl)%id=iBLK1;  pmap(iBLK1)%segs(nsr)%id=iBLK;
  pmap(iBLK)%segs(nsl)%ty=tyl  ;  pmap(iBLK1)%segs(nsr)%ty=tyr ;
  pmap(iBLK)%segs(nsl)%np=1    ;  pmap(iBLK1)%segs(nsr)%np=1   ;
  AddSeg=1
    else
  pmap(iBLK )%segs(nsl)%np=pmap(iBLK )%segs(nsl)%np+1
  pmap(iBLK1)%segs(nsr)%np=pmap(iBLK1)%segs(nsr)%np+1
  AddSeg=1
    endif
  end function AddSeg
  end module split_mod
