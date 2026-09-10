!
!  $Author: pkubota $
!  $Date: 2009/07/13 10:56:20 $
!  $Revision: 1.10 $
!
MODULE FileGrib

   USE Constants, ONLY : r4, r8, postclim,RegIntIn,kpds13,ENS,RunRecort,RecLat,RecLon,newlat0,newlat1,newlon0,newlon1
   USE Utils, ONLY : getpoint

   IMPLICIT NONE
   
   PUBLIC :: GDSPDSSETION, WriteGrbField

CONTAINS

SUBROUTINE GDSPDSSETION (kgds,kpds)
!
! only for gauss grid, regular in progress
! 
  USE Infgdspds, only : hh_fct, mm_fct, dd_fct, yy_fct  ,&
                        hh_anl, dd_anl, mm_anl, yy_anl
  USE Sizes,    ONLY : Imax, Jmax

  USE RegInterp, ONLY : Idim=>IdimOut, Jdim=>JdimOut,gLats,glond, glatd
 
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: kpds(200),kgds(200)
  INTEGER :: fct_hours,anl_hours
  integer :: idrt,dhfct
  CHARACTER(LEN=4) :: YYYY
  INTEGER :: YY1,YYA1
  INTEGER :: YY2,YYA2
  INTEGER :: CENTURY
  INTEGER :: YEAR_CENTURY

  IF(RunRecort)THEN
     newlat1=getpoint(glatd,RecLat(1))
     newlon0=getpoint(glond,RecLon(1))

     newlat0=getpoint(glatd,RecLat(2))
     newlon1=getpoint(glond,RecLon(2))
  END IF

  kpds=0
  kgds=0
!
! data representation type
!
  idrt=4                        ! gaussian
  if(RegIntIn)idrt=0            ! latlon
!
  CALL GetHours(hh_anl,dd_anl,mm_anl,yy_anl,&
                hh_fct,dd_fct,mm_fct,yy_fct,fct_hours,anl_hours)

!
!  GRID DESCRIPTION SECTION (GDS)
!
    kgds(1)=modulo(idrt,256)                    ! DATA REPRESENTATION TYPE
    IF(.NOT.RunRecort)THEN
       kgds(2)=Idim                                ! N(I) NR POINTS ON LATITUDE CIRCLE
       kgds(3)=Jdim                                ! N(J) NR POINTS ON LONGITUDE MERIDIAN        
       select case (idrt)
       case(0)
         kgds(4)=nint(-90.00*1000.)                ! LA(1) LATITUDE OF ORIGIN
       case(4)
         kgds(4)=nint(-90.00*1000.)                ! LA(1) LATITUDE OF ORIGIN
       end select
       kgds(5)=0                                   ! LO(1) LONGITUDE OF ORIGIN
       kgds(6)=128                                 ! RESOLUTION FLAG  (RIGHT ADJ COPY OF OCTET 17)
       kgds(7)=-kgds(4)                            ! LA(2) LATITUDE OF EXTREME POINT
       kgds(8)=nint(360000./kgds(2)*(kgds(2)-1))   ! LO(2) LONGITUDE OF EXTREME POINT
       kgds(9)=nint(360000./kgds(2))               ! DI LONGITUDINAL DIRECTION OF INCREMENT
       select case (idrt)
       case(0)
         kgds(10)=nint(180000./(kgds(3)-1))        ! DJ LATITUDINAL DIRECTION INCREMENT
       case(4)
         kgds(10)=jdim/2                           ! N - NR OF CIRCLES POLE TO EQUATOR
       end select                     
    ELSE
       kgds(2)=newlon1-newlon0+1                   ! N(I) NR POINTS ON LATITUDE CIRCLE
       kgds(3)=newlat1-newlat0+1                   ! N(J) NR POINTS ON LONGITUDE MERIDIAN        
       select case (idrt)
       case(0)
         kgds(4)=NINT(glatd(newlat0)*1000.0_r8)    ! LA(1) LATITUDE OF ORIGIN
       case(4)
         kgds(4)=NINT(gLats(newlat0)*1000.0_r8)    ! LA(1) LATITUDE OF ORIGIN
       end select
       kgds(5)=glond(newlon0)                      ! LO(1) LONGITUDE OF ORIGIN
       kgds(6)=128                                 ! RESOLUTION FLAG  (RIGHT ADJ COPY OF OCTET 17)
       kgds(7)=NINT(glatd(newlat1)*1000.0_r8)   ! LA(2) LATITUDE OF EXTREME POINT
       kgds(8)=NINT(glond(newlon1)*1000.0_r8)   ! LO(2) LONGITUDE OF EXTREME POINT
       kgds(9)=NINT((glond(newlon1) - glond(newlon0))*1000.0_r8/kgds(2))  ! DI LONGITUDINAL DIRECTION OF INCREMENT
       select case (idrt)
       case(0)
         kgds(10)=NINT((glatd(newlat1)-glatd(newlat0))*1000.0_r8/(kgds(3)))     ! DJ LATITUDINAL DIRECTION INCREMENT
       case(4)
         kgds(10)=jdim/2                           ! N - NR OF CIRCLES POLE TO EQUATOR
       end select                     
    END IF

    kgds(11:19)=0
    kgds(20)=255                                ! NEITHER ARE PRESENT
    kgds(21:200)=0

!
! PRODUCT DEFINITION SECTION (PDS)
!
  kpds(1)=46               ! ID of center (46--> CPTEC)
  kpds(2)=0                ! MODEL IDENTIFICATION
  kpds(3)=255              ! GRID IDENTIFICATION (RIGHT ADJ COPY OF OCTET 8)
  kpds(4)=128              ! GDS/BMS FLAGS
  kpds(5)=2                ! INDICATOR OF PARAMETER
  kpds(6)=1                ! TYPE OF LEVEL
  kpds(7)=102              ! HEIGHT/PRESSURE , ETC OF LEVEL
  IF(postclim)THEN
!     IF(yy_anl > 2000) THEN
!        kpds(8) =yy_fct-2000         ! YEAR INCLUDING (CENTURY-1) yy_anl 2001-2000= 1
!     ELSE
!        kpds(8) =yy_fct-1900         ! YEAR INCLUDING (CENTURY-1) yy_anl 1999-1900=99 
!     END IF
      WRITE(YYYY,'(I4.4)')yy_fct
      READ (YYYY(1:2),'(I2.2)')YY1
      READ (YYYY(3:4),'(I2.2)')YY2
      CENTURY=YY1
      YEAR_CENTURY=YY2
      If(YY2==0)YEAR_CENTURY=100
      If(YY2> 0)CENTURY=CENTURY+1
      kpds(8) =YEAR_CENTURY         ! YEAR INCLUDING (CENTURY-1) yy_anl 1999-1900=99 
      kpds(9) =mm_fct              ! MONTH OF YEAR   mm_anl,
      kpds(10)=dd_fct              ! DAY OF MONTH dd_anl,
      kpds(11)=hh_fct              ! HOUR OF DAY hh_anl,
  ELSE 
      !IF(yy_anl > 2000) THEN
      !   kpds(8) =yy_anl-2000         ! YEAR INCLUDING (CENTURY-1) yy_anl 2001-2000= 1
      !ELSE
      !   kpds(8) =yy_anl-1900         ! YEAR INCLUDING (CENTURY-1) yy_anl 1999-1900=99 
      !END IF
      WRITE(YYYY,'(I4.4)')yy_anl
      READ (YYYY(1:2),'(I2.2)')YY1
      READ (YYYY(3:4),'(I2.2)')YY2
      CENTURY=YY1
      YEAR_CENTURY=YY2
      IF(YY2==0)YEAR_CENTURY=100
      IF(YY2> 0)CENTURY=CENTURY+1
      kpds(8) =YEAR_CENTURY         ! YEAR INCLUDING (CENTURY-1) yy_anl 1999-1900=99 
      kpds(9) =mm_anl              ! MONTH OF YEAR   mm_anl,
      kpds(10)=dd_anl             ! DAY OF MONTH dd_anl,
      kpds(11)=hh_anl             ! HOUR OF DAY hh_anl,  
  END IF

  kpds(12)=0               ! MINUTE OF HOUR
!
! kpds(13) indicator of forecast time unit (warning should be automatic)
! 1 hour, 2 day, 3 month, 4 year, ... ON388 TABLE 4
! 10 3 hours, 11 6 hours, 12 12 hours, etc ...
! for CPTEC ensemble output is 6 hours
!
!Code Meaning table 4
!   0 Minute 
!   1 Hour 
!   2 Day 
!   3 Month 
!   4 Year 
!   5 Decade  (10 years) 
!   6 Normal  (30 years) 
!   7 Century (100 years) 
! 8-9 Reserved 
!  10 3 hours 
!  11 6 hours 
!  12 12 hours 
!  13-253 Reserved 
! 254 Second 

  kpds(13) = kpds13
  IF(kpds13==1)THEN
    dhfct=1               ! INDICATOR OF FORECAST TIME UNIT Hour 
  ELSE IF(kpds13==2)THEN
    dhfct=24              ! INDICATOR OF FORECAST TIME UNIT HOUR 
  ELSE IF(kpds13==3)THEN  
    dhfct=30*24           ! INDICATOR OF FORECAST TIME UNIT HOUR 
  ELSE IF(kpds13==11)THEN  
    dhfct= 6              ! INDICATOR OF FORECAST TIME UNIT HOUR 
  ELSE IF(kpds13==4)THEN  
    dhfct=365*24          ! INDICATOR OF FORECAST TIME UNIT HOUR 
  ELSE IF(kpds13==10)THEN  
    dhfct=3               ! INDICATOR OF FORECAST TIME UNIT HOUR 
  ELSE IF(kpds13==12)THEN  
    dhfct=12              ! INDICATOR OF FORECAST TIME UNIT HOUR 
  ELSE
    WRITE(6,*)'ERROR AT kpds13 INDICATOR OF FORECAST TIME UNIT'
    STOP
  END IF
  
  
  IF(postclim)THEN
  kpds(14)=0               ! TIME RANGE 1
  ELSE
  kpds(14)=(fct_hours-anl_hours)/dhfct  ! TIME RANGE 1
  END IF
  kpds(15)=0               ! TIME RANGE 2
  IF(postclim)THEN
    kpds(16)=1              ! TIME RANGE FLAG (1 Analysis)
  ELSE
    kpds(16)=0             ! TIME RANGE FLAG (0 Forecast product)
  ENDIF
  kpds(17)=0               ! NUMBER INCLUDED IN AVERAGE
  kpds(18)=1               ! VERSION NUMBER OF GRIB SPECIFICATION
  kpds(19)=254             ! VERSION NUMBER OF PARAMETER TABLE, key for use wgrib
  kpds(20)=0               ! NUMBER MISSING FROM AVERAGE/ACCUMULATION 
!  IF(yy_anl > 2000) THEN
!  kpds(21)=21              ! CENTURY OF REFERENCE TIME OF DATA
!  ELSE
!  kpds(21)=20              ! CENTURY OF REFERENCE TIME OF DATA
!  END IF
  kpds(21)=CENTURY         ! CENTURY OF REFERENCE TIME OF DATA
  kpds(22)=0               ! decimal scale factor, from table 1?
  kpds(23)=0               ! UNITS DECIMAL SCALE FACTOR
  kpds(24)=0               ! PDS BYTE 29, FOR NMC ENSEMBLE PRODUCTS
                           ! 128 IF FORECAST FIELD ERROR
                           ! 64 IF BIAS CORRECTED FCST FIELD
                           ! 32 IF SMOOTHED FIELD
                           ! WARNING: CAN BE COMBINATION OF MORE THAN 1
  kpds(25)=32              ! PDS BYTE 30, NOT USED
  kpds(26)=0
  kpds(27)=0
  kpds(28)=0
!
END SUBROUTINE GDSPDSSETION


SUBROUTINE WriteGrbField (var,ndim, kgds, kpds, bfr,l)

  USE tables, ONLY: table1,table2,table3,size_tb
  USE PrblSize, ONLY : Pmand
  USE Constants, ONLY : prefx,ENS

  IMPLICIT NONE
  CHARACTER(LEN=4), INTENT(IN   ) :: var
  INTEGER         , INTENT(IN   ) :: ndim
  integer         , INTENT(INOUT) :: kpds(200)
  integer         , INTENT(INOUT) :: kgds(200)
  integer                         :: kens(200)
  REAL (KIND=r8)  , INTENT(IN   ) :: bfr(:,:)
  INTEGER         , INTENT(IN   ) :: l

  integer :: ji,i,j,ijk,iret,ibs,nbits,kf, ifincr,ifhr,iMax
  real (kind=r4) :: var_min,var_max
  INTEGER :: itab1,itab2
     
  parameter(ji=1200*900) 
  logical bitmap(ndim)
  LOGICAL :: VarFound

  REAL (KIND=r8) :: bfrg(SIZE(bfr,1),SIZE(bfr,2))

  REAL (KIND=r4) :: bfr4(ndim)
  bfr4=0.0_r4
  bfrg=0.0_r8
  IF(RunRecort)THEN
      IF(RecLon(1) < 0.0_r8 .or. RecLon(2) < 0.0_r8)THEN
         bfrg=CSHIFT (bfr, SHIFT=SIZE(bfr,1)/2,DIM = 1) 
      ELSE
         bfrg=bfr
      END IF
      ijk=0
      DO j=newlat0,newlat1
         DO i=newlon0,newlon1
            ijk=ijk+1
            bfr4(ijk)=REAL(bfrg(i,j),r4)
         END DO
      END DO
  ELSE
      ijk=0
      DO j=1,SIZE(bfr,2)
         DO i=1,SIZE(bfr,1)
            ijk=ijk+1
            bfr4(ijk)=REAL(bfr(i,j),r4)
         END DO
      END DO
  END IF
  kf=ndim
  !
  ! Grib format
  !
  ! only for debug
  !
  var_min =  1.e8_r4
  var_max = -1.e8_r4

  ibs = 0
  nbits = 16     ! as used for Nilo
            
  var_min=MINVAL(bfr4(1:kf))
  var_max=MAXVAL(bfr4(1:kf))
    
  print*,'Parm = ',kpds(5),' Variable range:',minval(bfr4),' to ',maxval(bfr4)
  print*,'Level = ',Pmand(l)
! 
! only for purpuose debug
!
   !Search for var in table of vars
   VarFound=.FALSE.
   DO i=1,size_tb(1)
     IF(table1(i)%name==var) THEN
       VarFound=.TRUE.
       itab1=i
       EXIT
     END IF
   END DO
   IF(.NOT. VarFound) THEN
     PRINT *,'Variable '//var//' not found in table 1.'
     CALL flush(6)
     STOP
   END IF

   !Search for level in table of levels
   VarFound=.FALSE.
   DO i=1,size_tb(2)
     IF(trim(table2(i)%level_type)==trim(table1(itab1)%level)) THEN
       VarFound=.TRUE.
       itab2=i
       EXIT
     END IF
   END DO
   IF(.NOT. VarFound) THEN
     PRINT *,'Type of Level '//trim(table1(itab1)%level)//' not found in table 2.'     ! bug
     PRINT *,'May be one surface or plev variable, i.e. zonal wind, geopotential'
    kpds(5)=table1(itab1)%id    ! parameter identifier
    kpds(6)=100                 ! type of level	
    kpds(7)=Pmand(l)            ! level
!     CALL flush(6)
!     STOP
   ELSE
!
! surface variables
!
    kpds(5)=table1(itab1)%id
    kpds(6)=table2(itab2)%default
!    if(l.eq.1)kpds(7)=table2(itab2)%id           ! levels
    if(l.eq.1)kpds(7)=table2(itab2)%p2            ! levels
   END IF

   if(table1(itab1)%dec_scal_fact==-999)then      ! bug fixed
   ibs=0
   else
   ibs=table1(itab1)%dec_scal_fact
   endif

   nbits=table1(itab1)%precision
!
! only for output from old GCM: Precipitation and cloud cover
!
   do i=1,kf
   if(bfr4(i)<=1.e-15_r4.and.bfr4(i)>=-1.e-15_r4)bfr4(i)=0._r4  ! underflow if << 0., USE double precision
!    IF ( bfr4(i) <=1.e-15 ) then
!           bfr4(i) = 0.0
!    ENDIF
   enddo                          ! but w3lib should be in double precision
!
! correct accum variables (may be in GDSPDSSETION module )
! 
   IF ( var=='PREC' .or. var=='PRGE' .or. var=='PRCV' .or. var=='NEVE' ) THEN
!   if( kpds(13).eq.11 )IFINCR=6
!   if( kpds(13).eq.12 )IFINCR=12
   IFHR = kpds(14)
   kpds(14)= IFHR-1
   kpds(15)= IFHR
   kpds(16)= 4
!   ibs=-3.0		! to change table 1
   ENDIF 

   IF (ENS) THEN

!  write(6,*) 'Control and Perturbed forecast in coarse resolution'

!
!  ensemble information
!
   kens(1)=1            ! no change, 1=ensemble (identifies application)
!
!   kens(2)              type, 1=unperturbed crtl fcst, 2=ind neg,
!                        3=indv pos, 4=cluster,5=whole ensemble
!   kens(3)              ident number, if kens(2)=1 -> 1 high res, 2 low res
!                                      if kens(2)=2 or 3 -> will = 1, 2, 3, ...
!
!  control forecast
!
   if ( prefx=='AVN' ) then
   kens(2)=1
   kens(3)=1
   else
!
!  perturbed members
!
   if ( prefx=='01P' ) then
   kens(2)=3
   kens(3)=1
   elseif( prefx=='02P' ) then
   kens(2)=3
   kens(3)=2
   elseif( prefx=='03P' ) then
   kens(2)=3
   kens(3)=3
   elseif( prefx=='04P' ) then
   kens(2)=3
   kens(3)=4
   elseif( prefx=='05P' ) then
   kens(2)=3
   kens(3)=5
   elseif( prefx=='06P' ) then
   kens(2)=3
   kens(3)=6
   elseif( prefx=='07P' ) then
   kens(2)=3
   kens(3)=7
   elseif( prefx=='01N' ) then
   kens(2)=2
   kens(3)=1
   elseif( prefx=='02N' ) then
   kens(2)=2
   kens(3)=2
   elseif( prefx=='03N' ) then
   kens(2)=2
   kens(3)=3
   elseif( prefx=='04N' ) then
   kens(2)=2
   kens(3)=4
   elseif( prefx=='05N' ) then
   kens(2)=2
   kens(3)=5
   elseif( prefx=='06N' ) then
   kens(2)=2
   kens(3)=6
   elseif( prefx=='07N' ) then
   kens(2)=2
   kens(3)=7
   endif

   endif
!
   kens(4)=1            ! no change, 1=full field
   kens(5)=255          ! no change, 255 original resolution
!
! Write variable in grib format
!
   CALL putgben(51,kf,kpds,kgds,Kens,ibs,nbits,bitmap(1:kf),bfr4(1:kf),iret)
    
  if (iret.eq.0) then
!    write(6,*) 'PUTGBEN successful, iret=', iret
  else
    write(6,*) 'PUTGBEN failed!  iret=', iret
    stop 'WRITE_GRIB'
  endif
   ELSE

!   write(6,*) 'Determinist forecast in high resolution'
!
! Write variable in grib format
!
  CALL putgbn(51,kf,kpds,kgds,ibs,nbits,bitmap(1:kf),bfr4(1:kf),iret) ! control on # bytes
    
  if (iret.eq.0) then
    write(6,*) 'PUTGBN successful, iret=', iret
  else
    write(6,*) 'PUTGBN failed!  iret=', iret
    stop 'WRITE_GRIB'
  endif

   ENDIF

END SUBROUTINE WriteGrbField
SUBROUTINE GetHours(hh_anl,dd_anl,mm_anl,yy_anl,&
                 hh_fct,dd_fct,mm_fct,yy_fct,fct_hours,anl_hours)
  IMPLICIT NONE
  INTEGER, INTENT(in   ) :: hh_fct
  INTEGER, INTENT(in   ) :: dd_fct
  INTEGER, INTENT(in   ) :: mm_fct
  INTEGER, INTENT(in   ) :: yy_fct

  INTEGER, INTENT(in   ) :: hh_anl
  INTEGER, INTENT(in   ) :: dd_anl
  INTEGER, INTENT(in   ) :: mm_anl
  INTEGER, INTENT(in   ) :: yy_anl
  INTEGER, INTENT(OUT  ) :: fct_hours
  INTEGER, INTENT(OUT  ) :: anl_hours
  
  INTEGER :: iy,iday_anl,iday_fct,im
  INTEGER ::  MONLN(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  INTEGER ::  MONLB(12)=(/31,29,31,30,31,30,31,31,30,31,30,31/)
  iday_anl=0
  iday_fct=0
  DO iy=yy_anl,yy_fct
     IF(MOD(REAL(yy_fct),4.0) == 0.0) THEN
        IF(iy == yy_anl) THEN
           DO im=1,mm_anl-1
              iday_anl=iday_anl+MONLB(im)
           END DO
        END IF

        IF(iy == yy_fct) THEN
           DO im=1,mm_fct-1
              iday_fct=iday_fct+MONLB(im)
           END DO
        ELSE
           DO im=1,12
              iday_fct=iday_fct+MONLB(im)
           END DO
        END IF
     ELSE  
        IF(iy == yy_anl) THEN
           DO im=1,mm_anl-1
              iday_anl=iday_anl+MONLN(im)
           END DO
        END IF

        IF(iy == yy_fct) THEN
           DO im=1,mm_fct - 1
              iday_fct=iday_fct+MONLN(im)
           END DO
        ELSE
           DO im=1,12 
              iday_fct=iday_fct+MONLN(im)
           END DO
        END IF   
     END IF
  END DO
  
  iday_anl = iday_anl + (dd_anl-1) 
  iday_fct = iday_fct + (dd_fct-1) 
  anl_hours= 24*iday_anl + hh_anl
  fct_hours= 24*iday_fct + hh_fct
END SUBROUTINE GetHours
END MODULE FileGrib
