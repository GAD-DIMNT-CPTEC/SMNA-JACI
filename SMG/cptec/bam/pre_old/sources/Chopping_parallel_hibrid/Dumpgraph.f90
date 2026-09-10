MODULE Dumpgraph

  USE InputArrays, ONLY : glnPsInp, glnPsOut, gtopoinp, gtopoout,gUvelOut,gVvelOut,gTvirOut,gPresOut
  USE InputParameters, ONLY: i4, r4, r8, KmaxInp, KmaxInpp,kmaxout
  USE Sizes,     ONLY : imax, jmax, kmax, ibmax,jbmax,ibPerIJ, jbPerIJ
  USE Utils,     ONLY : lati
  USE Communications, ONLY : Collect_Grid_Full
  USE Parallelism, ONLY: myid 
  IMPLICIT NONE
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
     write(98,1010) (real(lati(jmax-k+1)),k=mi,min(mi+9,jmax))
     mi = mi+10
  enddo
  write(98,1002)'*'
  write(98,1041) 1 , 0.0
  write(98,1002)'*'
  write(98,1005)jgraf
  write(98,1002)'*'
  write(98,1006)'VARS 6'
  write(98,1023) kmaxout
  write(98,1022) kmaxout
  write(98,1001)'tpinp  0  99   topography          '
  write(98,1001)'tpout  0  99   topography          '
  write(98,1001)'psinp  0  99   Surface Pressure    '
  write(98,1001)'psout  0  99   Surface Pressure    '
  write(98,1007)'ENDVARS'
 1000  format(15A)
 1001  format(35A)
 1002  format(1A)
 1003  format('XDEF    ',I4,'  LINEAR     ', F10.5, F8.5)
 1005  format('TDEF ', I3, ' LINEAR 02JAN1987 1DY')
 1006  format(6A)
 1007  format(7A)
 1008  format('YDEF    ',I4, ' LEVELS  ')
 1009  format(1X,32F8.5)
 1010  format(1X,32F10.5)
 1020  format('U      ',I3,'  99   Zonal velocity')
 1021  format('V      ',I3,'  99   Meridional velocity')
 1022  format('T      ',I3,'  99   Temperature        ')
 1023  format('PRES   ',I3,'  99   Pressure           ')
 1041  format('ZDEF ',I4,' LEVELS    ',F8.5)
       close(98)
  END SUBROUTINE Writectl
  SUBROUTINE dumpgra
  !
  INTEGER :: i,j,k,ib,jb,jhalf
  real(kind=r8) :: ff(imax,jmax,kmaxout)
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
   Call Collect_Grid_Full(gpresOut,ff,kmaxout,0)
   if (myid.eq.0) then
     do k=1,kmaxout
     Call Writefield(imax*jmax,ff(:,:,1))
     enddo
   endif
   ff = 0.
   Call Collect_Grid_Full(gtvirOut,ff,kmaxout,0)
   if (myid.eq.0) then
     do k=1,kmaxout
     Call Writefield(imax*jmax,ff(:,:,1))
     enddo
   endif
   ff = 0.
   Call Collect_Grid_Full(gtopoinp,ff(:,:,1),1,0)
   if (myid.eq.0) then
     Call Writefield(imax*jmax,ff(:,:,1))
   endif
   ff = 0.
   Call Collect_Grid_Full(gtopoout,ff(:,:,1),1,0)
   if (myid.eq.0) then
     Call Writefield(imax*jmax,ff(:,:,1))
   endif
   ff = 0.
   Call Collect_Grid_Full(glnpsinp,ff(:,:,1),1,0)
   if (myid.eq.0) then
     Call Writefield(imax*jmax,ff(:,:,1))
   endif
   ff = 0.
   Call Collect_Grid_Full(glnpsout,ff(:,:,1),1,0)
   if (myid.eq.0) then
     Call Writefield(imax*jmax,ff(:,:,1))
   endif
 !
END SUBROUTINE dumpgra

SUBROUTINE WriteField (ndim, bfr)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ndim
  REAL(KIND=r8),  INTENT(IN) :: bfr(ndim)
  REAL(KIND=r4) :: bfr4(ndim)
  bfr4=REAL(bfr,r4)
  
  WRITE (UNIT=99) bfr4

END SUBROUTINE WriteField
 !
END MODULE Dumpgraph
