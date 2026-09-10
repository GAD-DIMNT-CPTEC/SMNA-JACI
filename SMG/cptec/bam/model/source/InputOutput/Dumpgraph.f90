MODULE Dumpgraph

  USE FieldsDynamics, ONLY : fgtmpm,fgum,fgvm,fgqm,fgdivm,fglnpm,fgzs,fgtmp,fgu,fgv,fgdiv,fglnps, &
                             fgqp, omg, fgq, fgyu, fgyv, fgtd, fgvdlnp, fgqd, fgprsl, fgphil, fgpsp, fgtmpp
    USE FieldsPhysics, ONLY:  PBL_CoefKm, PBL_CoefKh
  USE Utils,     ONLY : rcs2, colrad, lati
  USE Sizes,     ONLY : imax, jmax, kmax, ibmax,jbmax,ibPerIJ, jbPerIJ, eta
  USE Constants, ONLY :  pihalf, r8, r4
  USE Communications, ONLY : Collect_Grid_Full
  USE Parallelism, ONLY: myid 
   IMPLICIT NONE
  SAVE       

  PRIVATE
  PUBLIC :: dumpgra, Writectl
  CONTAINS

 ! dump fields for graphics
 ! 

  SUBROUTINE Writectl(jgraf)
  integer :: m,n,jgraf,k,mi,jini
  real :: h,ini,pi
  pi = acos(-1.)
  h = 360. / imax
  ini = 0.
  open(unit=98,file='gr.ctl')
  write(98,1000)'DSET   ^graphout'
  write(98,1000)'UNDEF  -2.56E33'
  write(98,1002)'*'
  write(98,1001)'OPTIONS SEQUENTIAL YREV BIG_ENDIAN '
  write(98,1002)'*'
  write(98,1001)'TITLE 1 Days of Sample Model Output'
  write(98,1002)'*'
  write(98,1003) imax, ini, h
  write(98,1002)'*'
  write(98,1008) jmax
  mi = 1
  do
     if (mi.gt.jmax) exit
     write(98,1010) (real(lati(k)/pi*180-90),k=mi,min(mi+9,jmax))
     mi = mi+10
  enddo
  write(98,1002)'*'
  write(98,1041) kmax
  mi = kmax
  do
     if (mi.lt.1) exit
     write(98,1009) (1000.*real(eta(k)),k=mi,max(mi-6,1),-1)
     mi = mi-7
  enddo
  write(98,1002)'*'
  write(98,1005)jgraf
  write(98,1002)'*'
  write(98,1006)'VARS 6'
  write(98,1020)kmax
  write(98,1021)kmax
  write(98,1022)kmax
  write(98,1023)kmax
  write(98,1001)'ps     0  99   Surface Pressure    '
  write(98,1001)'zs     0  99   Topography          '
  write(98,1007)'ENDVARS'
 1000  format(15A)
 1001  format(35A)
 1002  format(1A)
 1003  format('XDEF    ',I4,'  LINEAR     ', F10.5, F8.5)
 1005  format('TDEF ', I3, ' LINEAR 02JAN1987 1DY')
 1006  format(6A)
 1007  format(7A)
 1008  format('YDEF    ',I4, ' LEVELS  ')
 1009  format(1X,32F11.4)
 1010  format(1X,32F10.5)
 1020  format('U      ',I3,'  99   Zonal velocity')
 1021  format('V      ',I3,'  99   Meridional velocity')
 1022  format('T      ',I3,'  99   Temperature        ')
 1023  format('GEO    ',I3,'  99   Geopotential       ')
 1024  format('PRES   ',I3,'  99   Pressure           ')
 1025  format('Q      ',I3,'  99   Humidity           ')
 1041  format('ZDEF ',I4,' LEVELS    ')
       close(98)
  END SUBROUTINE Writectl
  SUBROUTINE dumpgra
  !
  INTEGER :: i,j,k,ib,jb,jhalf
  real(kind=r8) :: ff(imax,jmax,kmax+1)
  !
  CHARACTER*2 :: file(0:32)
  real(kind=r8) :: cosphi
  !
  DATA file  /'00','01','02','03','04','05','06','07','08', &
              '09','10','11','12','13','14','15', &
              '16','17','18','19','20','21','22', &
              '23','24','25','26','27','28','29', &
              '30','31','32'/
 !
   ff = 0.
   Call Collect_Grid_Full(fgu,ff,kmax,0)
!  Call Collect_Grid_Full(fgyu,ff,kmax,0)
   if (myid.eq.0) then
    DO j = 1,jMax
       cosphi = sin(colrad(j))
       ff(:,j,:) = ff(:,j,:)/cosphi
    ENDDO
    DO k=kmax,1,-1
       Call Writefield(imax*jmax,ff(:,:,k))
    ENDDO
   endif
   ff = 0.
   Call Collect_Grid_Full(fgv,ff,kmax,0)
!  Call Collect_Grid_Full(fgyv,ff,kmax,0)
   if (myid.eq.0) then
    DO j = 1,jMax
       cosphi = sin(colrad(j))
       ff(:,j,:) = ff(:,j,:)/cosphi
    ENDDO
    DO k=kmax,1,-1
       Call Writefield(imax*jmax,ff(:,:,k))
    ENDDO
   endif
   ff = 0.
   Call Collect_Grid_Full(fgtmp,ff,kmax,0)
!  Call Collect_Grid_Full(fgtd,ff,kmax,0)
   if (myid.eq.0) then
    DO k=kmax,1,-1
       Call Writefield(imax*jmax,ff(:,:,k))
    ENDDO
   endif
   ff = 0.
   Call Collect_Grid_Full(fgphil,ff,kmax,0)
!  Call Collect_Grid_Full(fgqd,ff,kmax,0)
   if (myid.eq.0) then
    DO k=kmax,1,-1
       Call Writefield(imax*jmax,ff(:,:,k))
    ENDDO
   endif
!  ff = 0.
!  Call Collect_Grid_Full(omg,ff,kmax,0)
!  Call Collect_Grid_Full(fgqd,ff,kmax,0)
!  if (myid.eq.0) then
!   DO k=kmax,1,-1
!      Call Writefield(imax*jmax,ff(:,:,k))
!   ENDDO
!  endif
!  ff = 0.
!  Call Collect_Grid_Full(fgqp,ff,kmax,0)
!  Call Collect_Grid_Full(fgqd,ff,kmax,0)
!  if (myid.eq.0) then
!   DO k=kmax,1,-1
!      Call Writefield(imax*jmax,ff(:,:,k))
!   ENDDO
!  endif
   ff = 0.
   Call Collect_Grid_Full(fgpsp,ff(:,:,1),1,0)
!  Call Collect_Grid_Full(fgvdlnp,ff(:,:,1),1,0)
   if (myid.eq.0) then
     Call Writefield(imax*jmax,ff(:,:,1))
   endif
   ff = 0.
   Call Collect_Grid_Full(fgzs,ff(:,:,1),1,0)
   if (myid.eq.0) then
     Call Writefield(imax*jmax,ff(:,:,1))
   endif
 !
!   WRITE(99) &
!      (((ffgu(i,j,k),i=1,imax),j=1,jmax), k=1,kmax),&
!      (((ffgv(i,j,k),i=1,imax),j=1,jmax), k=1,kmax),&
!      (((ft(i,j,k),i=1,imax),j=1,jmax), k=1,kmax),&
!      (((frot(i,j,k),i=1,imax),j=1,jmax), k=1,kmax),&
!      ((fp(i,j),i=1,imax),j=1,jmax)
 !
END SUBROUTINE dumpgra

SUBROUTINE WriteField (ndim, bfr)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ndim
  REAL(KIND=r8),  INTENT(IN) :: bfr(ndim)
  REAL(KIND=r4) :: bfr4(ndim)
  bfr4=REAL(bfr,r4)
  
  WRITE (UNIT=199) bfr4

END SUBROUTINE WriteField
 !
END MODULE Dumpgraph
