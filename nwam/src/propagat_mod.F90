#include "wavedef.h"
!#define USE_SMOOTH8
module propagat_mod
  use varcommon_mod
  use partition_mod
IMPLICIT NONE
  real   (8),allocatable::wkd(:)
  real   (8),allocatable::dxym(:,:)
  real (8),pointer::pvkdp(:,:,:)
  type  propagat_pre_type
    real(4) pab(4),pabp(4)
    integer*2 kpos(4),iquad,fill(3)
  end type propagat_pre_type
  type(propagat_pre_type),pointer:: pp_vs(:,:)
  integer mkj,iact,ia,ic
!NSMOOTH =4 | 8
#define NSMOOTH 4

 integer*1 LID(1)
 integer kjt,ixx,jyy,ixx1,jyy1
#define ICCG   1
#define IDSIDD 2
#define IMDPSP 2
! ICCG   = 1, cg                  | used implsch_mod & propagat_mod.F90
! IDSIDD = 2, dsidd               | used propagat_mod
contains
  SUBROUTINE InitPropagat !{
    integer k,iac,ia,ic,iad,icd,iaa,ica,i
    real(8) dxmt,dymt
    real(8) wkk,dk,wsk,cg,tanhdk
    ALLOCATE(wkd(kld+1)) ;
    mkj=kl*jnthet
    DO  k=1,kld !{
      if(k>1)wkd(k-1)=(wk(k)-wk(k-1))
    end do !}
    
    if(.not. associated(pp_vs))then
      ALLOCATE(pp_vs(kl*jnthet,0:nwpc))
    endif
    ALLOCATE(dxym(8,0:nwpc))
    if(.not.associated(pvkdp))ALLOCATE(pvkdp(kl,IMDPSP,0:nwpc))

    do iac=1,nwpc
      IF(nsp(iac)/=1) CYCLE
      call id2pos(iac,ia,ic)
      if(ia<=0.or.ic<=0)then
          write(6,*)mpi_id,'InitPropagat Error:',iac,nsp(iac),ia,ic
          stop
          call wav_abort('InitPropagat Error:1')
      endif
      !iad=adjix(ia-1);icd=adjiy(ic-1)
      !iaa=adjix(ia+1);ica=adjiy(ic+1)
      !if(icd<rectc(2).or.icd>rectc(4))icd=ic; !ZZZEERR
      !if(ica<rectc(2).or.ica>rectc(4))ica=ic; !ZZZEERR
      !   4---3---2
      !   |   |   |
      !   5---+---1
      !   |   |   |
      !   6---7---8
      !i  2  3  4
      !132534178576

      dymt=deg2m*grdszy
      dxmt=deg2m*grdszx
      dxym(1,iac)=  dxm(iac);  dxym(2,iac)=dym(iac)
      dxym(3,iac)=  dxm(iac);  dxym(4,iac)=dym(iac)
      dxym(5,iac)=  dxm(iac);  dxym(6,iac)=dym(iac)
      dxym(7,iac)=  dxm(iac);  dxym(8,iac)=dym(iac)
      DO  k=1,kl !{
        wkk=wk(k)
        dk=dep(iac)*wkk
        if(dk>40)then  !avoid float overflow
                       !  sinh(80)   =2.87e+34
                       !  cosh(40)**2=1.38e+34,
                       !1-tanh(40)   =3.61e-35
          wsk=sqrt(g*wkk)
          pvkdp(k,ICCG  ,iac)=0.5*wsk/wkk
          pvkdp(k,IDSIDD,iac)=0
        else
          wsk=sqrt(g*wkk*tanh(dk))
          pvkdp(k,ICCG  ,iac)=0.5*wsk*(1.+2.*dk/sinh(2.*dk))/wkk
          pvkdp(k,IDSIDD,iac)=0.5*g/cosh(dk)*wkk**2/wsk/cosh(dk)
        endif
      enddo
    end do
    call propagats_Pre
  end SUBROUTINE InitPropagat !}

#ifndef C_PROPAGAT
  SUBROUTINE propagats(e,ee,iacb,iace) !{
    REALD e(kl*jnthet,0:nwpc)
    REALD ee(kl*jnthet,0:nwpc)
    integer iacb,iace
    real(8) pab(4),pabp(4),et(6)
    integer kj,iac,kpos(4),ppos(3),iquad,iquad_pre
    type(propagat_pre_type),pointer:: ppv
    real(8) edd,eddm
    do iac=iacb,iace
      iquad_pre=-1
      DO  kj=1,jnthet*kl !{
        ppv=>pp_vs(kj,iac)
        iquad=ppv%iquad*3+1;
        ppos=ipos12(iquad:iquad+2,iac);
        et(1)=ppv%pab(1)*ee(ppv%kpos(1),iac    )+ppv%pab(2)*ee(ppv%kpos(2) ,iac    )+ppv%pab(3)*ee(ppv%kpos(3),iac    )+ppv%pab(4)*ee(ppv%kpos(4),iac    )
        et(2)=ppv%pab(1)*ee(ppv%kpos(1),ppos(1))+ppv%pab(2)*ee(ppv%kpos(2) ,ppos(1))+ppv%pab(3)*ee(ppv%kpos(3),ppos(1))+ppv%pab(4)*ee(ppv%kpos(4),ppos(1))
        et(3)=ppv%pab(1)*ee(ppv%kpos(1),ppos(2))+ppv%pab(2)*ee(ppv%kpos(2) ,ppos(2))+ppv%pab(3)*ee(ppv%kpos(3),ppos(2))+ppv%pab(4)*ee(ppv%kpos(4),ppos(2))
        et(4)=ppv%pab(1)*ee(ppv%kpos(1),ppos(3))+ppv%pab(2)*ee(ppv%kpos(2) ,ppos(3))+ppv%pab(3)*ee(ppv%kpos(3),ppos(3))+ppv%pab(4)*ee(ppv%kpos(4),ppos(3))
        iquad_pre=iquad
        et(5)=ppv%pabp(1)*et(1)+ppv%pabp(2)*et(2)+ppv%pabp(3)*et(3)+ppv%pabp(4)*et(4)
        if(et(5)<0)et(5)=0
        e(kj,iac)=et(5)
      end do
    end do
    call diffract_smooth(e,ee,iacb,iace)
  return
  end SUBROUTINE propagats !}
#endif
  SUBROUTINE GetIKJ(k,j,k_d,th_d,kpos,pab) !{
    integer k,j
    real(8) k_d,th_d
    real(8) wks,ths_
    real(4) pab(4)
    integer kpos(4)
    integer iwk,iwk1,jth,jth1
    real(8) p,q,r,ths
    integer ik,kj
    kj=(j-1)*kl+k
    if(abs(th_d)<1e-10)th_d=0
    ths=thet(j)+th_d;
    if(ths>=zpi)then
      ths=ths-zpi
    else if(ths<0)then
      ths=ths+zpi
    endif
    jth=idint(ths/deltth)+1;
    if(ths<thet(jth))then
      jth=jth-1
    endif
    jth1=jth+1;
    if(jth>=jnthet)then
      jth=jnthet
      jth1=1
    endif
    r=(ths-thet(jth))*ddeltth
    if(r>1)r=1
    if(abs(k_d)<1e-10)k_d=0
    wks=wk(k)+k_d
    if(k_d>=0)then
      if(wks>=wk(kld))then
        iwk=kl;iwk1=kl;
        p=wkh(kld-kl+1)
        q=0
      else
        iwk=k+1;    
        do while(wks>=wk(iwk))
          iwk=iwk+1;
        enddo
        iwk=iwk-1
        if(iwk<kl)then
          iwk1=iwk+1;
          q=(wks-wk(iwk))/wkd(iwk)
          p=1-q
        else if(iwk<kld)then
          q=(wks-wk(iwk))/wkd(iwk)
          p=(1-q)*wkh(iwk-kl+1)+q*wkh(iwk-kl+1+1)
          iwk=kl;iwk1=kl;
          q=0
        end if
      endif
    else
      if(wks<wk(1))then
        iwk=1;iwk1=1;
        if(wks>0)then
          q=wks/wkmin
          p=0;
        else
          p=0;q=0
        endif
      else
        iwk=k-1
        
        do while(iwk>0 .and. wks<wk(iwk))
          iwk=iwk-1;
        enddo
        q=(wks-wk(iwk))/wkd(iwk)
        iwk1=iwk+1;p=1-q
      endif
    endif
    pab(1)=p*(1-r);  pab(2)=q*(1-r)
    pab(3)=p*r    ;  pab(4)=q*r
    kpos(1)=(jth -1)*kl+iwk;  kpos(2)=(jth -1)*kl+iwk1;
    kpos(3)=(jth1-1)*kl+iwk;  kpos(4)=(jth1-1)*kl+iwk1;
    if(minval(pab)<-1e-10)then
      write(6,*)"K",k_d,th_d,q,r,iwk,iwk1,jth,jth1,wkd(iwk),wks,thet(jth),ths,zpi
    endif

  end SUBROUTINE GetIKJ !}
  SUBROUTINE GetIPOS(x_d,y_d,iac,iquad,pabp) !{
    real(8) x_d,y_d
    real(4) pabp(4)
    integer iac,iquad,ia,ic
    real(8) q,r
      !   4---3---2
      !   | 1 | 0 |
      !   5---+---1
      !   | 3 | 2 |
      !   6---7---8
    if(abs(x_d)<1e-5)x_d=0
    if(abs(y_d)<1e-5)y_d=0
    if(y_d>=0)then !{
      if(x_d>=0)then !{
        iquad=0
        q= x_d/dxm(iac)
        r= y_d/dym(iac)
      else     !}{
        iquad=1
        q=-x_d/dxm(iac)
        r= y_d/dym(iac)
      end if
    else
      if(x_d>=0)then !{
        iquad=2
        q= x_d/dxm(iac)
        r=-y_d/dym(iac)
      else
        iquad=3
        q=-x_d/dxm(iac)
        r=-y_d/dym(iac)
      end if !}
    end if !}
    pabp(1)=(1-q)*(1-r)
    pabp(2)=(q)*(1-r)
    pabp(3)=(1-q)*(r)
    pabp(4)=(q)*(r)
    if(minval(pabp)<-1e-10)then
      call id2pos(iac,ia,ic)
      write(6,*)"P",x_d,y_d,iquad,r,q,dxym(5,iac),dxym(6,iac),ia,ic
      call flush(6)
      stop
    endif
  end SUBROUTINE GetIPOS !}
  SUBROUTINE propagats_pre
  ! Local
    real(8) d0,dddx0,dddy0
    real(8) wk0,dk0,cg,cgx,cgy
    real(8) dsidd,ssrwk,ssrth,ssrwkd,ssrthd,ssrwkbc,ssrthbc
    real(4) pab(4),pabp(4)
    real(8) rt0,tRsd_tanLat
    integer iac,ia,ic,iquad,kpos(4),ppos(3)
    real(8)::x_d,y_d,k_d,th_d,rsec
    integer j,k,kj,isec
    integer nnt(10),nt
#if LOGSCURR==1
    real(8) duxdx0,duxdy0,duydx0,duydy0,ux,uy,dudsl,dudnl
    ux=0;uy=0;duxdx0=0;duxdy0=0;duydx0=0;duydy0=0
#endif
    nnt=0
    pab=0;pabp=0
    do iac=1,nwpc
      nt=0;
      iact=iac
      tRsd_tanLat=Rsd_tanLat(iac)
      d0   =dep(iac)
      dddx0=dddx(iac)
      dddy0=dddy(iac)
#if LOGSCURR==1
      if(logscurr/=0)then
        ux    =ucur(1,iac)
        uy    =ucur(2,iac)
        duxdx0=ucur(3,iac)
        duxdy0=ucur(4,iac)
        duydx0=ucur(5,iac)
        duydy0=ucur(6,iac)
      end if
#endif

      DO  j=1,jnthet !{
        ssrwkd=-( dddx0*cosths(j)+dddy0*sinths(j))
        ssrthd=-(-dddx0*sinths(j)+dddy0*cosths(j))
        rt0=cosths(j)*tRsd_tanLat !全球模式，大圆 系数  ! Big Circle
        ssrthbc=0 ;       ssrwkbc=0

#if LOGSCURR==1
          dudsl=duxdx0*cosths2(j)+(duxdy0+duydx0)*sincosths(j)+duydy0*sinths2(j)
          dudnl=(duydy0-duxdx0)*sincosths(j)+duxdy0*cosths2(j)-duydx0*sinths2(j)
          !波流 大圆  ! Big Circle
          ! !1:ux;2:uy;3:uxx;4:uxy;5:uyx;6:uyy
          ssrwkbc=(ux*sinths(j)-uy*cosths(j))*rt0
          ssrthbc=(ux*cosths(j)+uy*sinths(j))*rt0  ! + uy*tRsd_tanLat
#endif
        DO  k=1,kl !{
          wk0=wk(k)
          cg=pvkdp(k,ICCG,iac)
          ssrth=0 +cg*rt0 !大圆  ! Big Circle
          ssrwk=0;
          cgx=cg*cosths(j);cgy=cg*sinths(j)
#if LOGSCURR==1
            !******  1.  "the calculation of wave engery-current spreading"
            cgx=cgx+ux;cgy=cgy+uy
#endif
          !======================================================================c
          !******  2.  "the effect of refraction caused by topography and current"
            dsidd=pvkdp(k,IDSIDD,iac)
            ssrwk=ssrwk+ssrwkd*dsidd
            ssrth=ssrth+ssrthd*dsidd/wk0
          !if(iact==iacp(1))then
          !! write(3163,'(2i5,20e26.13)')k,j,ssrwk,dddx0,dddy0,dsidd
          !  !write(6,'(6i5," A",8e13.7)')k,j,kpos,wk0+k_d,pab
          !endif

#if LOGSCURR==1
          ssrwk=ssrwk-(dudsl*wk0)
          ssrth=ssrth-(dudnl) !*wk0/wk0
          ssrwk=ssrwk+ssrwkbc
          ssrth=ssrth+ssrthbc
#endif

          x_d=-deltts*cgx;y_d=-deltts*cgy
          k_d=-deltts*ssrwk;th_d=-deltts*ssrth

          call id2pos(iac,ia,ic)


          IF(abs(x_d)>dxm(iac).or.abs(y_d)>dym(iac).or.abs(th_d)>pi/2) THEN !{
            call id2pos(iac,ia,ic)
            DBGO(0,*) 'deltts too big, be careful !!! program have to stop !!!'
            DBGO(0,*) 'decrease deltts in nwamctl.ini please '
            DBGO(0,*) 'Note:86400( seconds of 1day ) MUST Divide Exactly deltts  '
            !86400=128*27*25  ;3600 =16*9*25 1440=32*9*5
            ! 1 2 4 8 16 32
            ! 3 6 12 24 48 96
            ! 9 18 36 72 144 288
            ! 5 10 20 40 80 160
            ! 15 30 60 120 240 480
            ! 45 90 180 360 720 1440
            ! 1 2 3 4 5 6 10 12 15 20 30 60 min
            !   8   9 16 18 24 32 36 40 45 48 72 80 90 96 120 144 160 180 240 288 360 480 720 1440 min
            ! 180 160 90 80 60 45 40 36 32 30
            ! 1 2 3 4 5 6 8 9 12 15 18,20,24,25,27,30,36,40,
            ! 32
            DBGO(0,*) 'at mpi_id:',mpi_id,iac,kj
            DBGO(0,*) 'ia=',ia,'ic=',ic,'k=',k,'j=',j,'iac=',iac,'x=',sxcord(iac),'y=',sycord(iac)
            DBGO(0,*) 'abs(x_d):abs(',x_d,')must <',dxm(iac)
            DBGO(0,*) 'abs(y_d):abs(',y_d,')must <',dym(iac)
            DBGO(0,*) 'abs(th_d):abs(',th_d,')Must <',pi/2
            rsec=abs(dxm(iac)/cgx);
            rsec=min(rsec,abs(dym(iac)/cgy));
            rsec=min(rsec,abs(pi/ssrth));
            isec=rsec;
            do while(mod(86400,isec)/=0)
              isec=isec-1
            enddo
            DBGO(0,*) 'try deltts from ',deltts,' to ',isec,'Sec',rsec,'sec'
            call wav_abort('deltts too big')
            stop
          end if !}
          call id2pos(iac,ia,ic)
          call GetIPOS(x_d,y_d,iac,iquad,pabp)
          call GetIKJ(k,j,k_d,th_d,kpos,pab)
          kj=(j-1)*kl+k
          if(mod(k,4)/=1)THEN
            if(kpos(1)-pp_vs(kj-1,iac)%kpos(1)/=1)THEN
              nt=nt+1;nnt(1)=nnt(1)+1
            endif
          endif
          pp_vs(kj,iac)%iquad=iquad ;
          pp_vs(kj,iac)%kpos =kpos  ! not -1 for fortran
          pp_vs(kj,iac)%pab  =pab   ;
          pp_vs(kj,iac)%pabp =pabp  ;
          if(minval(pab)<-1e-10.or.minval(pabp)<-1e-10)THEN
            write(6,'(4i5," A",16e11.3)')mpi_id,k,j,iac,k_d,th_d,thet(j)+th_d,x_d,y_d,pab,pabp
            call flush(6)
            stop
          endif
        end do
      end do
      if(nt/=0)THEN
        nnt(2)=nnt(2)+1
      endif
    end do
    !write(*,*)"propagats_Pre NO ALign:",mpi_id,nnt(1:2)
    call InitSmooth
  end SUBROUTINE propagats_Pre !}
  SUBROUTINE InitSmooth
    smooth_p0=1./16./8.
  end SUBROUTINE InitSmooth
  SUBROUTINE diffract_smooth(e,ee,nwpb,nwpe) !{
    REALD e(kl,jnthet,0:nwpc)
    REALD ee(kl,jnthet,0:nwpc)
    integer nwpb,nwpe
    real(8) et,ev
    integer iac,j,k,iac1,iac3,iac5,iac7
    real(8) sp1,sp3,sp5,sp7

#ifdef USE_SMOOTH8
    integer iac2,iac4,iac6,iac8
    real(8) sp2,sp4,sp6,sp8
#endif
      !   4---3---2
      !   |   |   |
      !   5---+---1
      !   |   |   |
      !   6---7---8
    !e_c_n=e_c_n*(1-nn*p0) +sum(e_o_t)*p0 +(e_c_n-e_c_o)*nn*p0
    !e_c_n=e_c_n +(sum(e_o_t)-e_c_o*nn)*p0
    !e_c_n=e_c_n +(sum(e_o_t-e_c_o))*p0
    do iac=nwpb,nwpe
      iac1=ipos8(iac)%i8(1);if(iac1>0)then;sp1=smooth_p0;else;sp1=0.;endif
      iac3=ipos8(iac)%i8(3);if(iac3>0)then;sp3=smooth_p0;else;sp3=0.;endif
      iac5=ipos8(iac)%i8(5);if(iac5>0)then;sp5=smooth_p0;else;sp5=0.;endif
      iac7=ipos8(iac)%i8(7);if(iac7>0)then;sp7=smooth_p0;else;sp7=0.;endif
#ifdef USE_SMOOTH8
      iac2=ipos8(iac)%i8(2);if(iac2>0)then;sp2=smooth_p0;else;sp2=0.;endif
      iac4=ipos8(iac)%i8(4);if(iac4>0)then;sp4=smooth_p0;else;sp4=0.;endif
      iac6=ipos8(iac)%i8(6);if(iac6>0)then;sp6=smooth_p0;else;sp6=0.;endif
      iac8=ipos8(iac)%i8(8);if(iac8>0)then;sp8=smooth_p0;else;sp8=0.;endif
#endif

      DO  j=1,jnthet !{
        DO  k=1,kl !{
            et=0;ev=ee(k,j,iac)
            et=et+(ee(k,j,iac1)-ev)*sp1
            et=et+(ee(k,j,iac3)-ev)*sp3
            et=et+(ee(k,j,iac5)-ev)*sp5
            et=et+(ee(k,j,iac7)-ev)*sp7
#ifdef USE_SMOOTH8
            et=et+(ee(k,j,iac2)-ev)*sp2
            et=et+(ee(k,j,iac4)-ev)*sp4
            et=et+(ee(k,j,iac6)-ev)*sp6
            et=et+(ee(k,j,iac8)-ev)*sp8
#endif
            e(k,j,iac)=e(k,j,iac)+et
      enddo
      enddo
    enddo
  end SUBROUTINE diffract_smooth !}
end module propagat_mod
