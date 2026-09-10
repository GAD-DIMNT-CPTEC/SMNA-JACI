!ctar
!ctar       Modified by Tarasova on December 2013- February 2014
!ctar
!ctar  -----------  analog of lwrad in Global model using Clirad-lw
!ctar
!ctar--------------- adding 2 artificial layers at the top: from 0 mb to
!ctar                 1 mb , and from 1 mb to the bottop of top model layer
!ctar-----------------------------------------------------------------------
!ctar         CALL cliradlw from RadiationDriver.f90
!ctar
!ctar         call cliradlw ( &
!ctar               ! Model Info and flags
!ctar              ncols, kmax  &    
!ctar               ! Atmospheric fields
!ctar        FlipPbot, FlipTe, FlipQe, FlipO3, gtg, &
!ctar                   ! LW Radiation fields 
!ctar               LwToaUpC, LwToaUp, asclr, asl, LwSfcNetC, LwSfcNet, &
!ctar               LwSfcDownC, LwSfcDown,    &
!ctar               ! Cloud field and Microphysics
!ctar               cld, clu )
!ctar
!ctar ------------------------------------------------------------------------    
!ctar              cliradlw call irrad
!ctar
!ctar              irrad call planck, plancd, h2oexps, conexps, co2exps, n20exps, 
!ctar              ch4exps,comexps, cfcexps, b10exps, tablup, h2okdis, co2kdis,
!ctar                   n2okdis, comkdis, cfckdis, b10kdis, cldovlp, sfcflux
!ctar
!ctar              sfcflux call planck, plancd
!ctar
!ctar------------------------------------------------------------------------   

MODULE Rad_Cliradlw

   USE Constants, ONLY :  r8

    USE Options, ONLY : co2val    ! co2val is wgne standard value in ppm


  IMPLICIT NONE
  PRIVATE

    !c-----parameters defining the size of the pre-computed tables for
    !c     transmittance using table look-up.

    !c     "nx" is the number of intervals in pressure
    !c     "no" is the number of intervals in o3 amount
    !c     "nc" is the number of intervals in co2 amount
    !c     "nh" is the number of intervals in h2o amount

    INTEGER, PARAMETER   ::  nx=26   
    INTEGER, PARAMETER   ::  no=21
    INTEGER, PARAMETER   ::  nc=30
    INTEGER, PARAMETER   ::  nh=31

    REAL(KIND=r8), DIMENSION(nx,nc) :: c1, c2, c3  

    REAL(KIND=r8), DIMENSION(nx,no) :: o1, o2, o3 

    REAL(KIND=r8), DIMENSION(nx,nh) :: h11, h12, h13, h21,h22,h23, h81,h82,h83    

    PUBLIC :: cliradlw,InitCliRadLW

CONTAINS
  SUBROUTINE InitCliRadLW()
    IMPLICIT NONE
    CALL read_table()
  END SUBROUTINE InitCliRadLW

  SUBROUTINE cliradlw ( &
                                ! Model Info and flags
       m         , &
       np        , &
                                ! Atmospheric fields
       pl20      , &
       Pint      , &
       tl        , &            !->   layer temperature (ta)                            k       m*np
       ql        , &            !->   layer specific humidity (wa)                      g/g     m*np
       o3l       , &            !->   layer ozone mixing ratio by mass (oa)             g/g     m*np
       tgg       , &            !->   land or ocean surface temperature (tg)            k       m*ns
                                ! LW Radiation fields 
       ulwclr    , &            !->   net downward flux, clear-sky (flc)             w/m**2  m*(np+1)
       ulwtop    , &            !->   net downward flux, all-sky   (flx)             w/m**2  m*(np+1)
       atlclr    , &            !->   Heating rate in clear case (K/s) compute cooling rate profile clear-sky (flc)   
       atl       , &            !->   Heating rate (K/s) compute cooling rate profile all-sky   (flx)
       rsclr     , &            !->   Net surface flux in clear case (W/m2 )specify output fluxes for lwrad clear-sky (flc)  
       rs        , &            !->   Net surface flux (W/m2 ) specify output fluxes for lwrad all-sky   (flx)
       dlwclr    , &            !->   Downward flux at surface in clear case (W/m2 )
       dlwbot    , &
                                ! Cloud field and Microphysics
       cld       , &
       clu         )     

    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar
    ! atlclr....Heating rate in clear case (K/s)
    ! atl.......Heating rate (K/s)
    ! rsclr.....Net surface flux in clear case (W/m2 )
    ! rs........Net surface flux
    ! dlwclr....Downward flux at surface in clear case (W/m2 )
    ! dlwbot....Downward flux at surface (W/m2 )

    ! Model Info and flags
    INTEGER, INTENT(in) :: m
    INTEGER, INTENT(in) :: np
    
    ! Atmospheric fields

    REAL(KIND=r8), INTENT(in) :: pl20(m,np)
    REAL(KIND=r8), INTENT(in) :: Pint(m,np+1)  
    REAL(KIND=r8), INTENT(in) :: tl  (m,np)
    REAL(KIND=r8), INTENT(in) :: ql  (m,np)
    REAL(KIND=r8), INTENT(in) :: o3l (m,np)

    REAL(KIND=r8), INTENT(in) :: tgg(m)

    ! LW Radiation fields 
    REAL(KIND=r8), INTENT(out) :: ulwclr (m)
    REAL(KIND=r8), INTENT(out) :: ulwtop (m)
    REAL(KIND=r8), INTENT(out) :: rsclr  (m)
    REAL(KIND=r8), INTENT(out) :: rs     (m)
    REAL(KIND=r8), INTENT(out) :: dlwclr (m)
    REAL(KIND=r8), INTENT(out) :: dlwbot (m)
    REAL(KIND=r8), INTENT(out) :: atlclr (m,np)
    REAL(KIND=r8), INTENT(out) :: atl    (m,np)

    ! Cloud field and Microphysics
    REAL(KIND=r8),    INTENT(in) :: cld  (m,np)
    REAL(KIND=r8),    INTENT(in) :: clu  (m,np)

    !c-----working parameters

    INTEGER  :: ict(m)
    INTEGER  :: icb(m)
    INTEGER  :: np1, i,j,k 

    INTEGER, PARAMETER  :: ns=1
    
    REAL(KIND=r8)   :: pl (m,np+2)      !c   level pressure (pl)                               mb      m*(np+1)
    REAL(KIND=r8)   :: ta (m,np+1)      !c   layer temperature (ta)                            k       m*np
    REAL(KIND=r8)   :: wa (m,np+1)     !c   layer specific humidity (wa)                      g/g     m*np
    REAL(KIND=r8)   :: oa (m,np+1)

    REAL(KIND=r8)   :: tb   (m)        
    REAL(KIND=r8)   :: fs   (m,ns)     
    REAL(KIND=r8)   :: tg   (m,ns)     
    REAL(KIND=r8)   :: eg   (m,ns,10)  
    REAL(KIND=r8)   :: taucl(m,np+1,3) 
    REAL(KIND=r8)   :: fcld (m,np+1)   
    REAL(KIND=r8)   :: flx  (m,np+2)   
    REAL(KIND=r8)   :: flc  (m,np+2)   
    REAL(KIND=r8)   :: sfcem(m)     
    REAL(KIND=r8)   :: dp   (m,np+1)   

    REAL(KIND=r8)   ::  co2

    !c ---- initialize output fluxes
    tb   =0.0_r8
    fs   =0.0_r8
    tg   =0.0_r8
    eg   =0.0_r8
    taucl=0.0_r8
    fcld =0.0_r8
    flx  =0.0_r8
    flc  =0.0_r8
    sfcem=0.0_r8
    dp   =0.0_r8

    ulwclr = 0.0_r8
    ulwtop = 0.0_r8
    atlclr = 0.0_r8
    atl    = 0.0_r8
    rsclr  = 0.0_r8
    rs     = 0.0_r8
    dlwclr = 0.0_r8 
    dlwbot = 0.0_r8

    !ctar     REAL(KIND=r8)      ::tck
    !ctar     REAL(KIND=r8)      ::tbk

    !c-----specify new ozone concentration

    !ctar       co2val=300.0_r8

    co2=co2val*1.0e-6_r8

 
    !c-----specify input errays in irrad

    !ctar      do k=1,np
    !ctar       do i=1,m
    !ctar         pl(i,k+1)=pl20(i,k)
    !ctar         ta(i,k)=tl(i,k)
    !ctar         wa(i,k)=ql(i,k)
    !ctar         oa(i,k)=o3l(i,k)
    !ctar       enddo
    !ctar      enddo
    
    DO k=1,np+1
       DO i=1,m          
          pl(i,k+1)=Pint(i,k)
       ENDDO
    ENDDO

    DO k=1,np
       DO i=1,m          
          ta(i,k+1)=tl(i,k)
          wa(i,k+1)=ql(i,k)
          oa(i,k+1)=o3l(i,k)
       ENDDO
    ENDDO

    !c adding of 2 layers at top 0 mb --1mb ---

    np1=np+1

    !c-----high pressure level

    DO i=1,m
       !tar         pl(i,1)=1.
       pl(i,1)=0.0_r8            !mb
       pl(i,2)=pl(i,3)*0.5_r8    !mb
       ta(i,1)=tl(i,1)-20.0_r8
       !tar                           
       wa(i,1)=1.0e-12_r8        
       IF (o3l(i,1).GE.2.5e-6_r8) THEN  
          oa(i,1)=o3l(i,1)
       ELSE 
          oa(i,1)=2.5e-6_r8        
       ENDIF
       !ctar                                    
    ENDDO

    !c-----specify near surface air temperature tb as mean of lowest layer and Tg

    DO i=1,m
       tb(i)=(ta(i,np+1)+tgg(i))*0.5_r8
       !c          tb(i)=tgg(i)
    ENDDO

    !c--------specify ground properties

    DO i=1,m
       fs(i,ns)=1.0_r8
       tg(i,ns)=tgg(i)
       DO j=1,10
          eg(i,ns,j)=1.0_r8 
       ENDDO
    ENDDO

    !c-----specify ict and icb levels

    !ctar          ict=26
    !ctar          icb=30

    ict=1
    icb=1

    !ctar      do i=1,m
    !ctar        do k=1,np1
    !ctar          if (pl(i,k)/pl(i,np1+1).le.0.4.and.pl(i,k+1)/pl(i,np1+1).gt.0.4) tck=k
    !ctar          if (pl(i,k)/pl(i,np1+1).le.0.7.and.pl(i,k+1)/pl(i,np1+1).gt.0.7) tbk=k
    !ctar        enddo
    !ctar          ict=ict+tck
    !ctar          icb=icb+tbk
    !ctar       enddo
    !ctar       ict=int(ict/m)
    !ctar       icb=int(icb/m)

    ! specify level indices separating high clouds from middle clouds
    ! (ict), and middle clouds from low clouds (icb).  this levels
    ! correspond to 400mb and 700 mb roughly.

    ! CPTEC-GCM works in sigma levels, hence in all columns the same
    ! layer will correspond to 0.4 and 0.7.Therefore, search is
    ! done only in the 1st column
    DO k=1,np
       DO i=1,m
          IF (pl(i,k)/pl(i,np1+1).LE.0.4_r8.AND.pl(i,k+1)/pl(i,np1+1).GT.0.4_r8) ict(i)=k
          IF (pl(i,k)/pl(i,np1+1).LE.0.7_r8.AND.pl(i,k+1)/pl(i,np1+1).GT.0.7_r8) icb(i)=k
       END DO
    ENDDO

    !c-----specify cloud optical depth

    DO k=1,np+1
       DO i=1,m
          dp(i,k)=pl(i,k+1)-pl(i,k)               
          taucl(i,k,1)=0.0_r8
          taucl(i,k,2)=0.0_r8              
          taucl(i,k,3)=0.0_r8 
          fcld(i,k)=0.0_r8                   
       ENDDO
    ENDDO

    DO k=9,np
       DO i=1,m
          IF (cld(i,k).GT.clu(i,k).AND.cld(i,k).GT.0.1_r8) THEN 
             fcld(i,k+1)=cld(i,k)

             IF (ta(i,k+1).GT.253.0_r8) THEN
                taucl(i,k+1,2)=0.05_r8*dp(i,k+1)*1.16_r8
             ELSE
                taucl(i,k+1,1)=0.025_r8*dp(i,k+1)*1.16_r8           
             ENDIF

          ENDIF
          IF (clu(i,k).GE.cld(i,k).AND.clu(i,k).GT.0.1_r8) THEN 
             fcld(i,k+1)=clu(i,k)
             taucl(i,k+1,2)=0.05_r8*dp(i,k+1)*1.16_r8
          ENDIF
       ENDDO
    ENDDO

    !c-----print out input data

    !       write (9,103)
    !       write (9,104)
    !       write (9,105)

    !      do i=1,m
    !       do k=1,np+1
    !         write (9,106) k,pl(i,k)
    !         write (9,107) ta(i,k),wa(i,k),oa(i,k), &
    !                      taucl(i,k,2),fcld(i,k)
    !       enddo
    !      enddo

    !         write (9,106) np+2,pl(m,np+2)

    !         write (9,*) 'ict=',ict,'icb=',icb 

    !         write (9,*) 'tb=', tb(1)
    !         write (9,*) 'tg=', tg(1,1)                  

    !c-----compute fluxes 

    CALL irrad ( &
         m      , & ! INTEGER      , INTENT(in   ) :: m
         np1    , & ! INTEGER      , INTENT(in   ) :: np1
         pl     , & ! REAL(KIND=r8), INTENT(in   ) :: pl    (m,np1+1) 
         ta     , & ! REAL(KIND=r8), INTENT(in   ) :: ta    (m,np1)
         wa     , & ! REAL(KIND=r8), INTENT(in   ) :: wa    (m,np1)
         oa     , & ! REAL(KIND=r8), INTENT(in   ) :: oa    (m,np1)
         tb     , & ! REAL(KIND=r8), INTENT(in   ) :: tb    (m)
         co2    , & ! REAL(KIND=r8), INTENT(in   ) :: co2
         taucl  , & ! REAL(KIND=r8), INTENT(in   ) :: taucl (m,np1,3) 
         fcld   , & ! REAL(KIND=r8), INTENT(in   ) :: fcld  (m,np1)
         ict    , & ! INTEGER      , INTENT(in   ) :: ict   (m)
         icb    , & ! INTEGER      , INTENT(in   ) :: icb   (m)
         ns     , & ! INTEGER      , INTENT(in   ) :: ns
         fs     , & ! REAL(KIND=r8), INTENT(in   ) :: fs    (m,ns)
         tg     , & ! REAL(KIND=r8), INTENT(in   ) :: tg    (m,ns)
         eg     , & ! REAL(KIND=r8), INTENT(in   ) :: eg    (m,ns,10)
         flx    , & ! REAL(KIND=r8), INTENT(out  ) :: flx   (m,np1+1)
         flc    , & ! REAL(KIND=r8), INTENT(out  ) :: flc   (m,np1+1)
         sfcem    ) ! REAL(KIND=r8), INTENT(out  ) :: sfcem (m)

    !c-----compute cooling rate profile
    !ctar
    DO k=1,np
       DO i=1,m
          !         atl(i,k)=-(flx(i,k+2)-flx(i,k+1))*8.441874/  &
          atl(i,k)=-(flx(i,k+2)-flx(i,k+1))*9.7706e-5_r8/  &                  
               dp(i,k+1)
          !         atlclr(i,k)=-(flc(i,k+2)-flc(i,k+1))*8.441874/  &
          atlclr(i,k)=-(flc(i,k+2)-flc(i,k+1))*9.7706e-5_r8/  &         
               dp(i,k+1)          
       ENDDO
    ENDDO

    DO k=1,np
       DO i=1,m
          IF(ABS(atlclr(i,k)).LT.1.e-22_r8) atlclr(i,k)=0.
          IF(ABS(atl(i,k)).LT.1.e-22_r8) atl(i,k) = 0.

          !        atlclr(i,k) = atlclr(i,k) * 1.1574e-5
          !        atl(i,k) = atl(i,k) * 1.1574e-5     
       ENDDO
    ENDDO

    !c-----specify output fluxes for lwrad
    !c   emission by the surface (sfcem)                w/m**2     m
    DO i=1,m         
       ulwclr(i)=-flc(i,1)
       ulwtop(i)=-flx(i,1)
       rsclr(i)=-flc(i,np+2)
       rs(i)=-flx(i,np+2)
       dlwclr(i)=flc(i,np+2)-sfcem(i)
       dlwbot(i)=flx(i,np+2)-sfcem(i)
    ENDDO

    !c-----print out results

    !       write (9,202)

    !       write (9,*) 'ulwclr=',ulwclr(1),'ulwtop=',ulwtop(1) 
    !       write (9,*) 'rsclr=',rsclr(1),'rs=',rs(1)
    !       write (9,*) 'dlwclr=',dlwclr(1),'dlwbot=',dlwbot(1)                     

    !       write (9,203)
    !       write (9,204)

    !      do i=1,m
    !       do k=1,np+2
    !           write (9,205) k,pl(i,k),flx(i,k),flc(i,k)
    !       enddo
    !      enddo
    !      
    !       write (9,207) sfcem(m)
    !       write (9,208) sfcem(m)-flx(m,np+2), sfcem(m)-flc(m,np+2)     

    !ctar      
    !      do i=1,m
    !       do k=1,np
    !           write (9,206) atl(i,k), atlclr(i,k)
    !       enddo
    !      enddo      
    !

!101 FORMAT (' high    = ',a6,/,' trace   = ',a6,/,' overcast= ',a6,/, &
!         ' cldwater= ',a6,/,' aerosol = ',a6)
!102 FORMAT (/,4x,'ts =',f7.2,' K',/,3x,'co2 = ',e9.3,' pppv',/,3x, &
!         'n2o = ',e9.3,' pppv',/,3x,'ch4 = ',e9.3,' pppv',/,3x, &
!         'f11 = ',e9.3,' pppv',/,3x,'f12 = ',e9.3, ' pppv',/,3x, &
!         'f22 = ',e9.3,' pppv')
!103 FORMAT (/,' ******  INPUT DATA  ******',/)
!104 FORMAT (10x,'p',7x,'T',7x,'q',8x,'o3',8x,'taucl',5x,'fcld', &
 !        4x,'taual(5)')
!105 FORMAT (9x,'(mb)',4x,'(k)',4x,'(g/g)',4x,'(g/g)',/)
!106 FORMAT (i4,f10.4)
!107 FORMAT (14x,f7.2,1p,2e10.3,2x,0p,f7.2,1x,f7.2)
!202 FORMAT (//,' ******  RESULTS  ******',/)
!203 FORMAT (10x,'p',8x,'flx',8x,'flc',5x,'dfdts',4x,'cooling rates')
!204 FORMAT (9x,'(mb)',4x,'(W/m^2)',4x,'(W/m^2)',2x,'(W/m^2)', &
!         6x,'(C/day)',/)
!205 FORMAT (i4,f10.4,f10.3,1x,f10.3)
!206 FORMAT (48x,f9.3,2x,f9.3)
!207 FORMAT (/,2x,'sfcem=',f8.2, ' W/m^2')
!208 FORMAT (/,2x,'flx_dn=',f8.2, ' W/m^2','flc_dn=',f8.2, ' W/m^2' )
!209 FORMAT (1x,f10.4,2(1x,f9.3))
!210 FORMAT (1x,f10.4,1x,f9.3)


  END SUBROUTINE cliradlw

  !c**********************   June 2003  *****************************

  SUBROUTINE irrad ( &
       m       , & ! INTEGER      , INTENT(in   ) :: m
       np      , & ! INTEGER      , INTENT(in   ) :: np
       pl      , & ! REAL(KIND=r8), INTENT(in   ) :: pl    (m,np+1) 
       ta      , & ! REAL(KIND=r8), INTENT(in   ) :: ta    (m,np)
       wa      , & ! REAL(KIND=r8), INTENT(in   ) :: wa    (m,np)
       oa      , & ! REAL(KIND=r8), INTENT(in   ) :: oa    (m,np)
       tb      , & ! REAL(KIND=r8), INTENT(in   ) :: tb    (m)
       co2     , & ! REAL(KIND=r8), INTENT(in   ) :: co2
       taucl   , & ! REAL(KIND=r8), INTENT(in   ) :: taucl (m,np,3) 
       fcld    , & ! REAL(KIND=r8), INTENT(in   ) :: fcld  (m,np)
       ict     , & ! INTEGER      , INTENT(in   ) :: ict   (m)
       icb     , & ! INTEGER      , INTENT(in   ) :: icb   (m)
       ns      , & ! INTEGER      , INTENT(in   ) :: ns
       fs      , & ! REAL(KIND=r8), INTENT(in   ) :: fs    (m,ns)
       tg      , & ! REAL(KIND=r8), INTENT(in   ) :: tg    (m,ns)
       eg      , & ! REAL(KIND=r8), INTENT(in   ) :: eg    (m,ns,10)
       flx     , & ! REAL(KIND=r8), INTENT(out  ) :: flx   (m,np+1)
       flc     , & ! REAL(KIND=r8), INTENT(out  ) :: flc   (m,np+1)
       sfcem     ) ! REAL(KIND=r8), INTENT(out  ) :: sfcem (m)
    !ctar      *                  high,trace,n2o,ch4,cfc11,cfc12,cfc22,
    !ctar      *                  vege,ns,fs,tg,eg,tv,ev,rv,
    !ctar      *                  overcast,cldwater,cwc,taucl,fcld,ict,icb,
    !ctar      *                  aerosol,na,taual,ssaal,asyal,
    !ctar      *                  flx,flc,dfdts,sfcem)

    !c*********************************************************************

    !c   THE EQUATION NUMBERS noted in this code follows the latest  
    !c    version (May 2003) of the NASA Tech. Memo. (2001), which can 
    !c    be accessed at ftp://climate.gsfc.nasa.gov/pub/chou/clirad_lw/

    !c*********************************************************************

    !c  CHANGE IN MAY 2003

    !c    The effective size of ice particles is replaced changed to the  
    !c    effective radius, and the definition of the effective radius 
    !c    follows that given in Chou, Lee and Yang (JGR, 2002).
    !c    The reff is no longer an input parameter.

    !c  CHANGE IN DECEMBER 2002

    !c    Do-loop 1500 is created inside the do-loop 1000 to compute 
    !c    the upward and downward emission of a layer

    !c   CHANGE IN JULY 2002

    !c    The effective Planck functions of a layer are separately
    !c    computed for the upward and downward emission (bu and bd).
    !c    For a optically thick cloud layer, the upward emission will be
    !c    at the cloud top temperature, and the downward emission will
    !c    at the cloud base temperature.

    !c   RECENT CHANGES:
    !c
    !c    Subroutines for planck functions
    !c    Subroutines for cloud overlapping
    !c    Eliminate "rflx" and "rflc". Fold the flux calculations in
    !c      Band 10 to that of the other bands.
    !c    Return the calculations when ibn=10 and trace=.false.
    !c    The number of aerosol types is allowed to be more than one.
    !c    Include sub-grid surface variability and vegetation canopy.
    !c    Include the CKD continuum absorption coefficient as an option.
    !c    
    !c********************************************************************

    !c Ice and liquid cloud particles are allowed to co-exist in each of the
    !c  np layers. 
    !c
    !c The maximum-random assumption is applied for cloud overlapping. 
    !c  Clouds are grouped into high, middle, and low clouds separated 
    !c  by the level indices ict and icb.  Within each of the three groups,
    !c  clouds are assumed maximally overlapped.  Clouds among the three 
    !c  groups are assumed randomly overlapped. The indices ict and icb 
    !c  correspond approximately to the 400 mb and 700 mb levels.
    !c
    !c Various types of aerosols are allowed to be in any of the np layers. 
    !c  Aerosol optical properties can be specified as functions of height  
    !c  and spectral band.
    !c
    !c The surface can be divided into a number of sub-regions either with or 
    !c  without vegetation cover. Reflectivity and emissivity can be 
    !c  specified for each sub-region.
    !c
    !c There are options for computing fluxes:
    !c
    !c   If high = .true., transmission functions in the co2, o3, and the
    !c   three water vapor bands with strong absorption are computed using
    !c   table look-up.  cooling rates are computed accurately from the
    !c   surface up to 0.01 mb.
    !c   If high = .false., transmission functions are computed using the
    !c   k-distribution method with linear pressure scaling for all spectral
    !c   bands except Band 5.  cooling rates are not accurately calculated 
    !c   for pressures less than 10 mb. the computation is faster with
    !c   high=.false. than with high=.true.
    !c
    !c   If trace = .true., absorption due to n2o, ch4, cfcs, and the 
    !c   two minor co2 bands in the window region is included.
    !c   Otherwise, absorption in those minor bands is neglected.
    !c
    !c   If vege=.true., a vegetation layer is added, and the emission and 
    !c   reflectivity are computed for the ground+vegetation surface.
    !c   Otherwise, only ground and ocean surfaces are considered.
    !c
    !c   If overcast=.true., the layer cloud cover is either 0 or 1.
    !c   If overcast=.false., the cloud cover can be anywhere between 0 and 1.
    !c   Computation is faster for the .true. option than the .false. option.
    !c
    !c   If cldwater=.true., taucl is computed from cwc and reff as a
    !c   function of height and spectral band. 
    !c   If cldwater=.false., taucl must be given as input to the radiation
    !c   routine. For this case, taucl is independent of spectral band.
    !c
    !c   If aerosol = .true., aerosols are included in calculating transmission
    !c   functions. Otherwise, aerosols are not included.
    !c   
    !c
    !c The IR spectrum is divided into nine bands:
    !c   
    !c   band     wavenumber (/cm)   absorber
    !c
    !c    1           0 - 340           h2o
    !c    2         340 - 540           h2o
    !c    3         540 - 800       h2o,cont,co2
    !c    4         800 - 980       h2o,cont
    !c                              co2,f11,f12,f22
    !c    5         980 - 1100      h2o,cont,o3
    !c                              co2,f11
    !c    6        1100 - 1215      h2o,cont
    !c                              n2o,ch4,f12,f22
    !c    7        1215 - 1380      h2o,cont
    !c                              n2o,ch4
    !c    8        1380 - 1900          h2o
    !c    9        1900 - 3000          h2o
    !c
    !c In addition, a narrow band in the 17 micrometer region (Band 10) is added
    !c    to compute flux reduction due to n2o
    !c
    !c    10        540 - 620       h2o,cont,co2,n2o
    !c
    !c Band 3 (540-800/cm) is further divided into 3 sub-bands :
    !c
    !c   subband   wavenumber (/cm)
    !c
    !c    3a        540 - 620
    !c    3b        620 - 720
    !c    3c        720 - 800
    !c
    !c---- Input parameters                               units    size
    !c
    !c   number of soundings (m)                            --      1
    !c   number of atmospheric layers (np)                  --      1
    !c   level pressure (pl)                               mb      m*(np+1)
    !c   layer temperature (ta)                            k       m*np
    !c   layer specific humidity (wa)                      g/g     m*np
    !c   layer ozone mixing ratio by mass (oa)             g/g     m*np
    !c   surface air temperature (tb)                      k        m
    !c   co2 mixing ratio by volume (co2)                  pppv     1
    !c   option (high) (see explanation above)              --      1
    !c   option (trace) (see explanation above)             --      1
    !c   n2o mixing ratio by volume (n2o)                  pppv     1
    !c   ch4 mixing ratio by volume (ch4)                  pppv     1
    !c   cfc11 mixing ratio by volume (cfc11)              pppv     1
    !c   cfc12 mixing ratio by volume (cfc12)              pppv     1
    !c   cfc22 mixing ratio by volume (cfc22)              pppv     1
    !c   option for including vegetation cover (vege)       --      1
    !c   number of sub-grid surface types (ns=2)              --    m
    !c   fractional cover of sub-grid regions (fs)       fraction  m*ns
    !c   land or ocean surface temperature (tg)            k       m*ns
    !c   land or ocean surface emissivity (eg)           fraction  m*ns*9
    !c   vegetation temperature (tv)                       k       m*ns
    !c   vegetation emissivity (ev)                      fraction  m*ns*9
    !c   vegetation reflectivity (rv)                    fraction  m*ns*9
    !c   option for cloud fractional cover                  --      1
    !c      (overcast)   (see explanation above)
    !c   option for cloud optical thickness                 --      1
    !c      (cldwater)   (see explanation above)
    !c   cloud water mixing ratio (cwc)                   gm/gm   m*np*3
    !c       index 1 for ice particles
    !c       index 2 for liquid drops
    !c       index 3 for rain drops
    !c   cloud optical thickness (taucl)                    --    m*np*3
    !c       index 1 for ice particles
    !c       index 2 for liquid drops
    !c       index 3 for rain drops
    !c   cloud amount (fcld)                             fraction  m*np
    !c   level index separating high and middle             --      1
    !c       clouds (ict)
    !c   level index separating middle and low              --      1
    !c       clouds (icb)
    !c   option for including aerosols (aerosol)            --      1
    !c   number of aerosol types (na)                       --      1
    !c   aerosol optical thickness (taual)                  --   m*np*10*na
    !c   aerosol single-scattering albedo (ssaal)           --   m*np*10*na
    !c   aerosol asymmetry factor (asyal)                   --   m*np*10*na
    !c
    !c---- output parameters
    !c
    !c   net downward flux, all-sky   (flx)             w/m**2  m*(np+1)
    !c   net downward flux, clear-sky (flc)             w/m**2  m*(np+1)
    !c   sensitivity of net downward flux  
    !c       to surface temperature (dfdts)            w/m**2/k m*(np+1)
    !c   emission by the surface (sfcem)                w/m**2     m
    !c
    !c Data used in table look-up for transmittance calculations:
    !c
    !c   c1 , c2, c3: for co2 (band 3)
    !c   o1 , o2, o3: for  o3 (band 5)
    !c   h11,h12,h13: for h2o (band 1)
    !c   h21,h22,h23: for h2o (band 2)
    !c   h81,h82,h83: for h2o (band 8)
    !c 
    !c Notes: 
    !c
    !c   (1) Scattering is parameterized for clouds and aerosols.
    !c   (2) Diffuse cloud and aerosol transmissions are computed
    !c       from exp(-1.66*tau).
    !c   (3) If there are no clouds, flx=flc.
    !c   (4) plevel(1) is the pressure at the top of the model atmosphere,
    !c        and plevel(np+1) is the surface pressure.
    !c   (5) Downward flux is positive and upward flux is negative.
    !c   (6) sfcem and dfdts are negative because upward flux is defined as negative.
    !c   (7) For questions and coding errors, please contact Ming-Dah Chou,
    !c       Code 913, NASA/Goddard Space Flight Center, Greenbelt, MD 20771.
    !c       Phone: 301-614-6192, Fax: 301-614-6307,
    !c       e-mail: chou@climate.gsfc.nasa.gov
    !c
    !c***************************************************************************

    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !ctar---- input parameters ------

    INTEGER      , INTENT(in   ) :: m
    INTEGER      , INTENT(in   ) :: np
    INTEGER      , INTENT(in   ) :: ict(m)
    INTEGER      , INTENT(in   ) :: icb(m)
    INTEGER      , INTENT(in   ) :: ns

    REAL(KIND=r8), INTENT(in) :: pl    (m,np+1) 
    REAL(KIND=r8), INTENT(in) :: ta    (m,np)
    REAL(KIND=r8), INTENT(in) :: wa    (m,np)
    REAL(KIND=r8), INTENT(in) :: oa    (m,np)
    REAL(KIND=r8), INTENT(in) :: tb    (m)
    REAL(KIND=r8), INTENT(in) :: co2
    REAL(KIND=r8), INTENT(in) :: fs    (m,ns)
    REAL(KIND=r8), INTENT(in) :: tg    (m,ns)
    REAL(KIND=r8), INTENT(in) :: eg    (m,ns,10)
    REAL(KIND=r8), INTENT(in) :: taucl (m,np,3) 
    REAL(KIND=r8), INTENT(in) :: fcld  (m,np)

    !ctar---- output parameters ------

    !ctar      real dfdts(m,np+1)

    REAL(KIND=r8), DIMENSION(m,np+1):: dfdts

    REAL(KIND=r8), INTENT(out) :: flx (m,np+1)
    REAL(KIND=r8), INTENT(out) :: flc (m,np+1)
    REAL(KIND=r8), INTENT(out) :: sfcem(m)


    !ctar---- input parameters defined inside the subroutine ------

    INTEGER, PARAMETER   :: na=1

    REAL(KIND=r8), PARAMETER :: n2o = 0.28e-6_r8
    REAL(KIND=r8), PARAMETER :: ch4 =1.75e-6_r8
    REAL(KIND=r8), PARAMETER :: cfc11 = 0.3e-9_r8
    REAL(KIND=r8), PARAMETER :: cfc12 = 0.5e-9_r8
    REAL(KIND=r8), PARAMETER :: cfc22 = 0.2e-9_r8

    !ctar ------------ not used at the moment-------------

    REAL(KIND=r8), DIMENSION(m,ns) :: tv
    REAL(KIND=r8), DIMENSION(m,ns,10) :: ev, rv   
    REAL(KIND=r8), DIMENSION(m,np,3) :: cwc
    REAL(KIND=r8), DIMENSION(m,np,10,na) :: taual, ssaal, asyal  

    !ctar---- input parameters defined inside the subroutine ------

    LOGICAL, PARAMETER  :: high = .FALSE. 
    LOGICAL, PARAMETER  :: vege = .FALSE.
    LOGICAL, PARAMETER  :: trace = .TRUE.
    LOGICAL, PARAMETER  :: overcast = .FALSE. 
    LOGICAL, PARAMETER  :: cldwater = .FALSE.
    LOGICAL, PARAMETER  :: aerosol= .FALSE.

    !c---- static data -----

    REAL(KIND=r8), PARAMETER, DIMENSION(9) :: xkw = (/ & 
         29.55_r8, 4.167e-1_r8, 1.328e-2_r8, 5.250e-4_r8, &
         5.25e-4_r8, 9.369e-3_r8, 4.719e-2_r8, 1.320e-0_r8, 5.250e-4_r8 /)
    REAL(KIND=r8), PARAMETER, DIMENSION(9) :: xke = (/ & 
         0.0_r8,    271.0_r8,    25.00_r8,   16.8_r8, &
         8.31_r8,   6.52_r8,    12.7_r8,    0.0_r8,  0.0_r8 /)
    REAL(KIND=r8), PARAMETER, DIMENSION(9) :: aw = (/ & 
         0.0021_r8, 0.0140_r8, 0.0167_r8, 0.0302_r8, &
         0.0307_r8, 0.0195_r8, 0.0152_r8, 0.0008_r8, 0.0096_r8 /)
    REAL(KIND=r8), PARAMETER, DIMENSION(9) :: bw = (/ & 
         -1.01e-5_r8, 5.57e-5_r8, 8.54e-5_r8, 2.96e-4_r8, &
         2.86e-4_r8, 1.108e-4_r8, 7.608e-5_r8, -3.52e-6_r8, 1.64e-5_r8 /)
    REAL(KIND=r8), PARAMETER, DIMENSION(9) :: pm = (/ & 
         1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, 0.77_r8, 0.5_r8, 1.0_r8, 1.0_r8 /)
    REAL(KIND=r8), PARAMETER, DIMENSION(6,9) :: fkw =  RESHAPE(  &
         SHAPE = (/ 6, 9 /), SOURCE = (/ & 
         0.2747_r8,0.2717_r8,0.2752_r8,0.1177_r8,0.0352_r8,0.0255_r8, &
         0.1521_r8,0.3974_r8,0.1778_r8,0.1826_r8,0.0374_r8,0.0527_r8, &
         1.00_r8, 1.00_r8,1.00_r8,1.00_r8,1.00_r8,1.00_r8, &
         0.4654_r8,0.2991_r8,0.1343_r8,0.0646_r8,0.0226_r8,0.0140_r8, &
         0.5543_r8,0.2723_r8,0.1131_r8,0.0443_r8,0.0160_r8,0.0000_r8, &
         0.5955_r8,0.2693_r8,0.0953_r8,0.0335_r8,0.0064_r8,0.0000_r8, &
         0.1958_r8,0.3469_r8,0.3147_r8,0.1013_r8,0.0365_r8,0.0048_r8, &
         0.0740_r8,0.1636_r8,0.4174_r8,0.1783_r8,0.1101_r8,0.0566_r8, &
         0.1437_r8,0.2197_r8,0.3185_r8,0.2351_r8,0.0647_r8,0.0183_r8 /) )
    REAL(KIND=r8), PARAMETER, DIMENSION(6,3) :: gkw = RESHAPE(  &
         SHAPE = (/ 6, 3 /), SOURCE = (/ & 
         0.1782_r8,0.0593_r8,0.0215_r8,0.0068_r8,0.0022_r8,0.0000_r8, &
         0.0923_r8,0.1675_r8,0.0923_r8,0.0187_r8,0.0178_r8,0.0000_r8, &
         0.0000_r8,0.1083_r8,0.1581_r8,0.0455_r8,0.0274_r8,0.0041_r8 /) )

    REAL(KIND=r8), PARAMETER, DIMENSION(3,10) :: aib  = RESHAPE(  &
         SHAPE = (/ 3, 10 /), SOURCE = (/ & 
         -0.44171_r8,    0.61222_r8,   0.06465_r8, &
         -0.13727_r8,    0.54102_r8,   0.28962_r8, &
         -0.01878_r8,    1.19270_r8,   0.79080_r8, &
         -0.01896_r8,    0.78955_r8,   0.69493_r8, &
         -0.04788_r8,    0.69729_r8,   0.54492_r8, &
         -0.02265_r8,    1.13370_r8,   0.76161_r8, &
         -0.01038_r8,    1.46940_r8,   0.89045_r8, &
         -0.00450_r8,    1.66240_r8,   0.95989_r8, &
         -0.00044_r8,    2.01500_r8,   1.03750_r8, &
         -0.02956_r8,    1.06430_r8,   0.71283_r8 /) ) 

    REAL(KIND=r8), PARAMETER, DIMENSION(4,10) :: awb = RESHAPE(  &
         SHAPE = (/ 4, 10 /), SOURCE = (/ & 
         0.08641_r8,    0.01769_r8,    -1.5572e-3_r8,   3.4896e-5_r8, &
         0.22027_r8,    0.00997_r8,    -1.8719e-3_r8,   5.3112e-5_r8, &
         0.38074_r8,   -0.03027_r8,     1.0154e-3_r8,  -1.1849e-5_r8, &
         0.15587_r8,    0.00371_r8,    -7.7705e-4_r8,   2.0547e-5_r8, &
         0.05518_r8,    0.04544_r8,    -4.2067e-3_r8,   1.0184e-4_r8, &
         0.12724_r8,    0.04751_r8,    -5.2037e-3_r8,   1.3711e-4_r8, &
         0.30390_r8,    0.01656_r8,    -3.5271e-3_r8,   1.0828e-4_r8, &
         0.63617_r8,   -0.06287_r8,     2.2350e-3_r8,  -2.3177e-5_r8, &
         1.15470_r8,   -0.19282_r8,     1.2084e-2_r8,  -2.5612e-4_r8, &
         0.34021_r8,   -0.02805_r8,     1.0654e-3_r8,  -1.5443e-5_r8 /) )
    REAL(KIND=r8), PARAMETER, DIMENSION(4,10)  :: aiw = RESHAPE(  &
         SHAPE = (/ 4, 10 /), SOURCE = (/ & 
         0.17201_r8,    1.8814e-2_r8,  -3.5117e-4_r8,   2.1127e-6_r8, &
         0.81470_r8,   -4.1989e-3_r8,   2.3152e-7_r8,   2.0992e-7_r8, &
         0.54859_r8,   -7.4266e-4_r8,   1.2865e-5_r8,  -5.7092e-8_r8, &
         0.39218_r8,    6.4180e-3_r8,  -1.1567e-4_r8,   6.9710e-7_r8, &
         0.71773_r8,   -5.1754e-3_r8,   4.6658e-5_r8,  -1.2085e-7_r8, &
         0.77345_r8,   -8.4966e-3_r8,   1.1451e-4_r8,  -5.5170e-7_r8, &
         0.74975_r8,   -8.7083e-3_r8,   1.3367e-4_r8,  -7.1603e-7_r8, &
         0.69011_r8,   -6.9766e-3_r8,   1.1674e-4_r8,  -6.6472e-7_r8, &
         0.83963_r8,   -1.0347e-2_r8,   1.4651e-4_r8,  -7.5965e-7_r8, &
         0.64860_r8,   -4.4142e-3_r8,   6.5458e-5_r8,  -3.2655e-7_r8 /) )
    REAL(KIND=r8), PARAMETER, DIMENSION(4,10)  :: aww = RESHAPE(  &
         SHAPE = (/ 4, 10 /), SOURCE = (/ & 
         -7.8566e-2_r8,  8.0875e-2_r8,  -4.3403e-3_r8,   8.1341e-5_r8, &
         -1.3384e-2_r8,  9.3134e-2_r8,  -6.0491e-3_r8,   1.3059e-4_r8, &
         3.7096e-2_r8,  7.3211e-2_r8,  -4.4211e-3_r8,   9.2448e-5_r8, &
         -3.7600e-3_r8,  9.3344e-2_r8,  -5.6561e-3_r8,   1.1387e-4_r8, &
         0.40212_r8,    7.8083e-2_r8,  -5.9583e-3_r8,   1.2883e-4_r8, &
         0.57928_r8,    5.9094e-2_r8,  -5.4425e-3_r8,   1.2725e-4_r8, &
         0.68974_r8,    4.2334e-2_r8,  -4.9469e-3_r8,   1.2863e-4_r8, &
         0.80122_r8,    9.4578e-3_r8,  -2.8508e-3_r8,   9.0078e-5_r8, &
         1.02340_r8,   -2.6204e-2_r8,   4.2552e-4_r8,   3.2160e-6_r8, &
         0.05092_r8,    7.5409e-2_r8,  -4.7305e-3_r8,   1.0121e-4_r8 /) )
    REAL(KIND=r8), PARAMETER, DIMENSION(4,10)  :: aig = RESHAPE(  &
         SHAPE = (/ 4, 10 /), SOURCE = (/ & 
         0.57867_r8,    1.5592e-2_r8,  -2.6372e-4_r8,   1.5125e-6_r8, &
         0.72259_r8,    4.7922e-3_r8,  -4.7164e-5_r8,   2.0400e-7_r8, &
         0.76109_r8,    6.9922e-3_r8,  -1.0935e-4_r8,   5.9885e-7_r8, &
         0.86934_r8,    4.2268e-3_r8,  -7.4085e-5_r8,   4.3547e-7_r8, &
         0.89103_r8,    2.8482e-3_r8,  -3.9174e-5_r8,   2.0098e-7_r8, &
         0.86325_r8,    3.2935e-3_r8,  -3.9872e-5_r8,   1.8015e-7_r8, &
         0.85064_r8,    3.8505e-3_r8,  -4.9259e-5_r8,   2.3096e-7_r8, &
         0.86945_r8,    3.7869e-3_r8,  -5.6525e-5_r8,   3.0016e-7_r8, &
         0.80122_r8,    4.9086e-3_r8,  -5.8831e-5_r8,   2.6367e-7_r8, &
         0.73290_r8,    7.3898e-3_r8,  -1.0515e-4_r8,   5.4034e-7_r8 /) )
    REAL(KIND=r8), PARAMETER, DIMENSION(4,10)  :: awg = RESHAPE(  &
         SHAPE = (/ 4, 10 /), SOURCE = (/ & 
         -0.51930_r8,    0.20290_r8,    -1.1747e-2_r8,   2.3868e-4_r8, &
         -0.22151_r8,    0.19708_r8,    -1.2462e-2_r8,   2.6646e-4_r8, &
         0.14157_r8,    0.14705_r8,    -9.5802e-3_r8,   2.0819e-4_r8, &
         0.41590_r8,    0.10482_r8,    -6.9118e-3_r8,   1.5115e-4_r8, &
         0.55338_r8,    7.7016e-2_r8,  -5.2218e-3_r8,   1.1587e-4_r8, &
         0.61384_r8,    6.4402e-2_r8,  -4.6241e-3_r8,   1.0746e-4_r8, &
         0.67891_r8,    4.8698e-2_r8,  -3.7021e-3_r8,   9.1966e-5_r8, &
         0.78169_r8,    2.0803e-2_r8,  -1.4749e-3_r8,   3.9362e-5_r8, &
         0.93218_r8,   -3.3425e-2_r8,   2.9632e-3_r8,  -6.9362e-5_r8, &
         0.01649_r8,    0.16561_r8,    -1.0723e-2_r8,   2.3220e-4_r8 /)  )
    INTEGER, PARAMETER, DIMENSION(9)  :: mw = (/ & 
         6,6,8,6,6,8,9,6,16 /)

!   !c-----parameters defining the size of the pre-computed tables for
!   !c     transmittance using table look-up.

!   !c     "nx" is the number of intervals in pressure
!   !c     "no" is the number of intervals in o3 amount
!   !c     "nc" is the number of intervals in co2 amount
!   !c     "nh" is the number of intervals in h2o amount

!   INTEGER, PARAMETER   ::  nx=26   
!   INTEGER, PARAMETER   ::  no=21
!   INTEGER, PARAMETER   ::  nc=30
!   INTEGER, PARAMETER   ::  nh=31

!   REAL(KIND=r8), DIMENSION(nx,nc) :: c1, c2, c3  

!   REAL(KIND=r8), DIMENSION(nx,no) :: o1, o2, o3 

!   REAL(KIND=r8), DIMENSION(nx,nh) :: h11, h12, h13, h21,h22,h23, h81,h82,h83    


    !c---- temporary arrays -----

    REAL(KIND=r8), DIMENSION(m,np) :: pa, dt
    REAL(KIND=r8), DIMENSION(m) :: tx, xlayer
    REAL(KIND=r8), DIMENSION(m,np,3) :: reff

    REAL(KIND=r8), DIMENSION(m) :: x1, x2,  x3

    REAL(KIND=r8), DIMENSION(m,np) :: dh2o, dcont, dco2, do3 

    REAL(KIND=r8), DIMENSION(m,np) :: dn2o, dch4

    REAL(KIND=r8), DIMENSION(m,np) :: df11, df12, df22

    REAL(KIND=r8), DIMENSION(m,6) :: th2o
    REAL(KIND=r8), DIMENSION(m,3) :: tcon
    REAL(KIND=r8), DIMENSION(m,6,2) :: tco2

    REAL(KIND=r8), DIMENSION(m,4) :: tn2o   
    REAL(KIND=r8), DIMENSION(m,4) :: tch4
    REAL(KIND=r8), DIMENSION(m,6) :: tcom

    REAL(KIND=r8), DIMENSION(m) :: tf11, tf12, tf22 

    REAL(KIND=r8), DIMENSION(m,np,6) :: h2oexp
    REAL(KIND=r8), DIMENSION(m,np,3) :: conexp   
    REAL(KIND=r8), DIMENSION(m,np,6,2) :: co2exp

    REAL(KIND=r8), DIMENSION(m,4) :: n2oexp(m,np,4)
    REAL(KIND=r8), DIMENSION(m,4) :: ch4exp(m,np,4)
    REAL(KIND=r8), DIMENSION(m,6) :: comexp(m,np,6)

    REAL(KIND=r8), DIMENSION(m,np) :: f11exp, f12exp, f22exp

    REAL(KIND=r8), DIMENSION(m,0:np+1) :: blayer
    REAL(KIND=r8), DIMENSION(m,np+1) :: blevel

    REAL(KIND=r8), DIMENSION(m,0:np+1) :: bd
    REAL(KIND=r8), DIMENSION(m,0:np+1) :: bu

    REAL(KIND=r8), DIMENSION(m) :: bs, dbs, rflxs

    REAL(KIND=r8), DIMENSION(m,np) ::  dp   
    REAL(KIND=r8), DIMENSION(m,np,3) :: cwp

    REAL(KIND=r8), DIMENSION(m,np+1) :: trant
    REAL(KIND=r8), DIMENSION(m) :: tranal(m)
    REAL(KIND=r8), DIMENSION(m,np+1) :: transfc
    REAL(KIND=r8), DIMENSION(m,np+1) :: trantcr

    REAL(KIND=r8), DIMENSION(m,np+1) :: flxu, flxd, flcu, flcd 

    REAL(KIND=r8), DIMENSION(m,np) :: taua, ssaa, asya, taerlyr   

    INTEGER, DIMENSION(m)  :: it, im, ib

    INTEGER, DIMENSION(m,np)  :: itx, imx, ibx 

    REAL(KIND=r8), DIMENSION(m) :: cldhi, cldmd, cldlw 

    REAL(KIND=r8), DIMENSION(m,np) :: tcldlyr
    REAL(KIND=r8), DIMENSION(m,np+1) :: fclr

    INTEGER  :: i,j,k,ibn,ik,iq,isb,k1,k2,ne

    REAL(KIND=r8) :: x, xx, p1, dwe, dpe, a1, b1, fk1, a2, b2, fk2 

    REAL(KIND=r8) :: yy

    REAL(KIND=r8) :: w1, w2, w3, g1, g2, g3, ww, gg, ff, tauc 

    LOGICAL :: oznbnd, co2bnd, h2otbl, conbnd, n2obnd 

    LOGICAL :: ch4bnd, combnd, f11bnd , f12bnd, f22bnd, b10bnd 


    !c-----coefficients for computing effective particle size following
    !c     McFarquhar (QJRMS, 2000)

    REAL(KIND=r8), PARAMETER, DIMENSION(5)  :: ai = (/ &
         2.076_r8, 2.054_r8, 2.035_r8, 2.019_r8, 2.003_r8 /)
    REAL(KIND=r8), PARAMETER, DIMENSION(5)  :: bi = (/ &
         0.148_r8, 0.130_r8, 0.119_r8,    0.111_r8,    0.102_r8 /)
    REAL(KIND=r8), PARAMETER, DIMENSION(5)  :: ci = (/ &
         -0.0453_r8, -0.0491_r8, -0.0507_r8, -0.0517_r8, -0.0532_r8 /)
    REAL(KIND=r8), PARAMETER, DIMENSION(5)  :: di = (/ &
         -0.00686_r8, -0.00711_r8, -0.00716_r8, -0.00717_r8, -0.00725_r8 /)


    !c-----xkw is the absorption coefficient for the first k-distribution
    !c     interval due to water vapor line absorption (Table 4)
    !c     Units are cm**2/g    

    !ctar      data xkw / 29.55_r8  , 4.167e-1_r8, 1.328e-2_r8, 5.250e-4_r8, &
    !ctar                5.25e-4_r8, 9.369e-3_r8, 4.719e-2_r8, 1.320e-0_r8, 5.250e-4_r8/

    !c-----xke is the absorption coefficient for the first k-distribution
    !c     function due to water vapor continuum absorption (Table 9).
    !c     Units are cm**2/g

    !c-----Roberts et al's continuum k data

    !c     data xke /  0.00,   339.00,  27.40,   15.8, &
    !c                9.40,   7.75,    7.70,    0.0,   0.0/

    !c-----CKD (Version 2.3) continuum k data

    !ctar      data xke /  0.0_r8,    271.0_r8,    25.00_r8,   16.8_r8, &
    !ctar                 8.31_r8,   6.52_r8,    12.7_r8,    0.0_r8,  0.0_r8/

    !c-----mw is the ratio between neighboring absorption coefficients
    !c     for water vapor line absorption (Table 4).

    !ctar       data mw /6,6,8,6,6,8,9,6,16/

    !c-----aw and bw (Table 3) are the coefficients for temperature scaling
    !c     for water vapor in Eq. (4.2).

    !ctar      data aw/ 0.0021_r8, 0.0140_r8, 0.0167_r8, 0.0302_r8, &
    !ctar              0.0307_r8, 0.0195_r8, 0.0152_r8, 0.0008_r8, 0.0096_r8/
    !ctar      data bw/ -1.01e-5_r8, 5.57e-5_r8, 8.54e-5_r8, 2.96e-4_r8, &
    !ctar               2.86e-4_r8, 1.108e-4_r8, 7.608e-5_r8, -3.52e-6_r8, 1.64e-5_r8/

    !c-----pm is the pressure-scaling parameter for water vapor absorption
    !c     Eq. (4.1) and Table 3.

    !ctar      data pm/ 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, 0.77_r8, 0.5_r8, 1.0_r8, 1.0_r8/

    !c-----fkw is the planck-weighted k-distribution function due to h2o
    !c     line absorption (Table 4).
    !c     The k-distribution function for Band 3, fkw(*,3), 
    !c     is not used (see the parameter gkw below).

    !ctar      data fkw / 0.2747_r8,0.2717_r8,0.2752_r8,0.1177_r8,0.0352_r8,0.0255_r8, &
    !ctar                0.1521_r8,0.3974_r8,0.1778_r8,0.1826_r8,0.0374_r8,0.0527_r8, &
    !ctar                 6*1.00_r8, &
    !ctar                 0.4654_r8,0.2991_r8,0.1343_r8,0.0646_r8,0.0226_r8,0.0140_r8, &
    !ctar                 0.5543_r8,0.2723_r8,0.1131_r8,0.0443_r8,0.0160_r8,0.0000_r8, &
    !ctar                 0.5955_r8,0.2693_r8,0.0953_r8,0.0335_r8,0.0064_r8,0.0000_r8, &
    !ctar                 0.1958_r8,0.3469_r8,0.3147_r8,0.1013_r8,0.0365_r8,0.0048_r8, &
    !ctar                 0.0740_r8,0.1636_r8,0.4174_r8,0.1783_r8,0.1101_r8,0.0566_r8, &
    !ctar                 0.1437_r8,0.2197_r8,0.3185_r8,0.2351_r8,0.0647_r8,0.0183_r8/

    !c-----gkw is the planck-weighted k-distribution function due to h2o
    !c     line absorption in the 3 subbands (800-720,620-720,540-620 /cm)
    !c     of band 3 (Table 10).  Note that the order of the sub-bands
    !c     is reversed.

    !ctar       data gkw/  0.1782_r8,0.0593_r8,0.0215_r8,0.0068_r8,0.0022_r8,0.0000_r8, &
    !ctar                0.0923_r8,0.1675_r8,0.0923_r8,0.0187_r8,0.0178_r8,0.0000_r8, &
    !ctar                 0.0000_r8,0.1083_r8,0.1581_r8,0.0455_r8,0.0274_r8,0.0041_r8/


    !c-----Coefficients for computing the extinction coefficient
    !c     for cloud ice particles (Table 11a, Eq. 6.4a).
    !c
    !ctar       data aib /  -0.44171_r8,    0.61222_r8,   0.06465_r8, &
    !ctar                  -0.13727_r8,    0.54102_r8,   0.28962_r8, &
    !ctar                  -0.01878_r8,    1.19270_r8,   0.79080_r8, &
    !ctar                  -0.01896_r8,    0.78955_r8,   0.69493_r8, &
    !ctar                 -0.04788_r8,    0.69729_r8,   0.54492_r8, &
    !ctar                  -0.02265_r8,    1.13370_r8,   0.76161_r8, &
    !ctar                  -0.01038_r8,    1.46940_r8,   0.89045_r8, &
    !ctar                  -0.00450_r8,    1.66240_r8,   0.95989_r8, &
    !ctar                  -0.00044_r8,    2.01500_r8,   1.03750_r8, &
    !ctar                  -0.02956_r8,    1.06430_r8,   0.71283_r8/
    !c
    !c-----coefficients for computing the extinction coefficient
    !c     for cloud liquid drops. (Table 11b, Eq. 6.4b)
    !c
    !ctar      data awb /   0.08641_r8,    0.01769_r8,    -1.5572e-3_r8,   3.4896e-5_r8, &
    !ctar                  0.22027_r8,    0.00997_r8,    -1.8719e-3_r8,   5.3112e-5_r8, &
    !ctar                  0.38074_r8,   -0.03027_r8,     1.0154e-3_r8,  -1.1849e-5_r8, &
    !ctar                  0.15587_r8,    0.00371_r8,    -7.7705e-4_r8,   2.0547e-5_r8, &
    !ctar                  0.05518_r8,    0.04544_r8,    -4.2067e-3_r8,   1.0184e-4_r8, &
    !ctar                  0.12724_r8,    0.04751_r8,    -5.2037e-3_r8,   1.3711e-4_r8, &
    !ctar                  0.30390_r8,    0.01656_r8,    -3.5271e-3_r8,   1.0828e-4_r8, &
    !ctar                  0.63617_r8,   -0.06287_r8,     2.2350e-3_r8,  -2.3177e-5_r8, &
    !ctar                  1.15470_r8,   -0.19282_r8,     1.2084e-2_r8,  -2.5612e-4_r8, &
    !ctar                 0.34021_r8,   -0.02805_r8,     1.0654e-3_r8,  -1.5443e-5_r8/
    !c
    !c-----coefficients for computing the single-scattering albedo
    !c     for cloud ice particles. (Table 12a, Eq. 6.5)
    !c
    !ctar      data aiw/    0.17201_r8,    1.8814e-2_r8,  -3.5117e-4_r8,   2.1127e-6_r8, &
    !ctar                  0.81470_r8,   -4.1989e-3_r8,   2.3152e-7_r8,   2.0992e-7_r8, &
    !ctar                  0.54859_r8,   -7.4266e-4_r8,   1.2865e-5_r8,  -5.7092e-8_r8, &
    !ctar                  0.39218_r8,    6.4180e-3_r8,  -1.1567e-4_r8,   6.9710e-7_r8, &
    !ctar                  0.71773_r8,   -5.1754e-3_r8,   4.6658e-5_r8,  -1.2085e-7_r8, &
    !ctar                  0.77345_r8,   -8.4966e-3_r8,   1.1451e-4_r8,  -5.5170e-7_r8, &
    !ctar                  0.74975_r8,   -8.7083e-3_r8,   1.3367e-4_r8,  -7.1603e-7_r8, &
    !ctar                  0.69011_r8,   -6.9766e-3_r8,   1.1674e-4_r8,  -6.6472e-7_r8, &
    !ctar                  0.83963_r8,   -1.0347e-2_r8,   1.4651e-4_r8,  -7.5965e-7_r8, &
    !ctar                  0.64860_r8,   -4.4142e-3_r8,   6.5458e-5_r8,  -3.2655e-7_r8/

    !c-----coefficients for computing the single-scattering albedo
    !c     for cloud liquid drops. (Table 12b, Eq. 6.5)
    !c
    !ctar      data aww/   -7.8566e-2_r8,  8.0875e-2_r8,  -4.3403e-3_r8,   8.1341e-5_r8, &
    !ctar                 -1.3384e-2_r8,  9.3134e-2_r8,  -6.0491e-3_r8,   1.3059e-4_r8, &
    !ctar                  3.7096e-2_r8,  7.3211e-2_r8,  -4.4211e-3_r8,   9.2448e-5_r8, &
    !ctar                 -3.7600e-3_r8,  9.3344e-2_r8,  -5.6561e-3_r8,   1.1387e-4_r8, &
    !ctar                  0.40212_r8,    7.8083e-2_r8,  -5.9583e-3_r8,   1.2883e-4_r8, &
    !ctar                  0.57928_r8,    5.9094e-2_r8,  -5.4425e-3_r8,   1.2725e-4_r8, &
    !ctar                  0.68974_r8,    4.2334e-2_r8,  -4.9469e-3_r8,   1.2863e-4_r8, &
    !ctar                  0.80122_r8,    9.4578e-3_r8,  -2.8508e-3_r8,   9.0078e-5_r8, &
    !ctar                  1.02340_r8,   -2.6204e-2_r8,   4.2552e-4_r8,   3.2160e-6_r8, &
    !ctar                  0.05092_r8,    7.5409e-2_r8,  -4.7305e-3_r8,   1.0121e-4_r8/ 
    !c
    !c-----coefficients for computing the asymmetry factor for cloud ice 
    !c     particles. (Table 13a, Eq. 6.6)
    !c
    !ctar      data aig /   0.57867_r8,    1.5592e-2_r8,  -2.6372e-4_r8,   1.5125e-6_r8, &
    !ctar                  0.72259_r8,    4.7922e-3_r8,  -4.7164e-5_r8,   2.0400e-7_r8, &
    !ctar                  0.76109_r8,    6.9922e-3_r8,  -1.0935e-4_r8,   5.9885e-7_r8, &
    !ctar                  0.86934_r8,    4.2268e-3_r8,  -7.4085e-5_r8,   4.3547e-7_r8, &
    !ctar                  0.89103_r8,    2.8482e-3_r8,  -3.9174e-5_r8,   2.0098e-7_r8, &
    !ctar                  0.86325_r8,    3.2935e-3_r8,  -3.9872e-5_r8,   1.8015e-7_r8, &
    !ctar                  0.85064_r8,    3.8505e-3_r8,  -4.9259e-5_r8,   2.3096e-7_r8, &
    !ctar                  0.86945_r8,    3.7869e-3_r8,  -5.6525e-5_r8,   3.0016e-7_r8, &
    !ctar                  0.80122_r8,    4.9086e-3_r8,  -5.8831e-5_r8,   2.6367e-7_r8, &
    !ctar                  0.73290_r8,    7.3898e-3_r8,  -1.0515e-4_r8,   5.4034e-7_r8/
    !c
    !c-----coefficients for computing the asymmetry factor for cloud liquid 
    !c     drops. (Table 13b, Eq. 6.6)
    !c
    !ctar       data awg /  -0.51930_r8,    0.20290_r8,    -1.1747e-2_r8,   2.3868e-4_r8, &
    !ctar                  -0.22151_r8,    0.19708_r8,    -1.2462e-2_r8,   2.6646e-4_r8, &
    !ctar                  0.14157_r8,    0.14705_r8,    -9.5802e-3_r8,   2.0819e-4_r8, &
    !ctar                   0.41590_r8,    0.10482_r8,    -6.9118e-3_r8,   1.5115e-4_r8, &
    !ctar                   0.55338_r8,    7.7016e-2_r8,  -5.2218e-3_r8,   1.1587e-4_r8, &
    !ctar                   0.61384_r8,    6.4402e-2_r8,  -4.6241e-3_r8,   1.0746e-4_r8, &
    !ctar                   0.67891_r8,    4.8698e-2_r8,  -3.7021e-3_r8,   9.1966e-5_r8, &
    !ctar                   0.78169_r8,    2.0803e-2_r8,  -1.4749e-3_r8,   3.9362e-5_r8, &
    !ctar                   0.93218_r8,   -3.3425e-2_r8,   2.9632e-3_r8,  -6.9362e-5_r8, &
    !ctar                  0.01649_r8,    0.16561_r8,    -1.0723e-2_r8,   2.3220e-4_r8/ 
    !c
    !c-----include tables used in the table look-up for co2 (band 3), 
    !c     o3 (band 5), and h2o (bands 1, 2, and 8) transmission functions.
    !c     "co2.tran4" is the co2 transmission table applicable to a large
    !c     range of co2 amount (up to 100 times of the present-time value).

    !INCLUDE "h2o.tran3_90"
    !INCLUDE "co2.tran4_90"
    !INCLUDE "o3.tran3_90"
    !ctar
    !ctar   Definition of some input parameters:     
    !ctar 
    !       high = .false. 
    !       trace = .true.
    !ctar       trace = .false.          
    !      n2o = 0.28e-6_r8
    !       ch4 = 1.75e-6_r8
    !       cfc11 = 0.3e-9_r8
    !       cfc12 = 0.5e-9_r8
    !       cfc22 = 0.2e-9_r8
    !       vege = .false.
    !       overcast = .false. 
    !      cldwater = .false. 
    !       aerosol = .false.             
    !ctar
    !----------------initialize variables currently not used

    tv=0.0_r8
    ev=0.0_r8
    rv=0.0_r8
    cwc=0.0_r8
    taual=0.0_r8
    ssaal=0.0_r8
    asyal=0.0_r8

    !ctar
    !c-----compute layer pressure (pa) and layer temperature minus 250K (dt)

    DO k=1,np
       DO i=1,m
          pa(i,k)=0.5_r8*(pl(i,k)+pl(i,k+1))
          dt(i,k)=ta(i,k)-250.0_r8
       ENDDO
    ENDDO

    !c-----compute layer absorber amount

    !c     dh2o : water vapor amount (g/cm**2)
    !c     dcont: scaled water vapor amount for continuum absorption
    !c            (g/cm**2)
    !c     dco2 : co2 amount (cm-atm)stp
    !c     do3  : o3 amount (cm-atm)stp
    !c     dn2o : n2o amount (cm-atm)stp
    !c     dch4 : ch4 amount (cm-atm)stp
    !c     df11 : cfc11 amount (cm-atm)stp
    !c     df12 : cfc12 amount (cm-atm)stp
    !c     df22 : cfc22 amount (cm-atm)stp
    !c     the factor 1.02 is equal to 1000/980
    !c     factors 789 and 476 are for unit conversion
    !c     the factor 0.001618 is equal to 1.02/(.622*1013.25) 
    !c     the factor 6.081 is equal to 1800/296

    DO k=1,np
       DO i=1,m

          dp   (i,k) = pl(i,k+1)-pl(i,k)

          dh2o (i,k) = 1.02_r8*wa(i,k)*dp(i,k)
          dh2o (i,k) = MAX(dh2o (i,k),1.e-10_r8)
          do3  (i,k) = 476.0_r8*oa(i,k)*dp(i,k)
          do3 (i,k) = MAX(do3 (i,k),1.e-6_r8)
          dco2 (i,k) = 789.0_r8*co2*dp(i,k)
          dco2 (i,k) = MAX(dco2 (i,k),1.e-4_r8)

          dch4 (i,k) = 789.0_r8*ch4*dp(i,k)
          dn2o (i,k) = 789.0_r8*n2o*dp(i,k)
          df11 (i,k) = 789.0_r8*cfc11*dp(i,k)
          df12 (i,k) = 789.0_r8*cfc12*dp(i,k)
          df22 (i,k) = 789.0_r8*cfc22*dp(i,k)

          !c-----compute scaled water vapor amount for h2o continuum absorption
          !c     following eq. (4.21).

          xx=pa(i,k)*0.001618_r8*wa(i,k)*wa(i,k)*dp(i,k)
          dcont(i,k) = xx*EXP(1800.0_r8/ta(i,k)-6.081_r8)

       ENDDO
    ENDDO

    !c-----Set default values for reff.
    !c     Index is 1 for ice, 2 for waterdrops and 3 for raindrops.

    DO k=1,np
       DO i=1,m
          reff(i,k,1)=40.0_r8
          reff(i,k,2)=10.0_r8
       ENDDO
    ENDDO

    !c-----compute layer cloud water amount (gm/m**2)

    IF (cldwater) THEN
       DO k=1,np
          DO i=1,m
             xx=1.02_r8*10000.0_r8*(pl(i,k+1)-pl(i,k))
             cwp(i,k,1)=xx*cwc(i,k,1)
             cwp(i,k,2)=xx*cwc(i,k,2)
             cwp(i,k,3)=xx*cwc(i,k,3)
          ENDDO
       ENDDO

       DO k=1,np
          DO i=1,m

             IF (cwp(i,k,1) .GT. 0.000001_r8) THEN

                !c-----Compute effective radius of ice cloud particles following Equation (6.9)

                j=INT((ta(i,k)-193.0_r8)*0.1_r8)
                IF (j.LT.1) j=1
                IF (j.GT.5) j=5

                !c-----Conversion of the unit of cwc in g/g to the unit of x in g/m^3.
                !c     The constant 348.43 is equal to (100/0.287), 
                !c     where the constant 0.287 is related to the gas constant of dry air.

                x=cwc(i,k,1)*348.43_r8*pa(i,k)/ta(i,k)
                x=LOG10(x)
                reff(i,k,1)=0.65_r8*10.0_r8**(ai(j)+bi(j)*x+ci(j)*x*x+di(j)*x*x*x)
                reff(i,k,1)=MAX(reff(i,k,1),10.0_r8)
                reff(i,k,1)=MIN(reff(i,k,1),70.0_r8)

             ENDIF

             IF(cwp(i,k,2) .GT. 0.000001_r8) THEN

                !c-----Effective radius of water cloud particles following Equation (6.13).

                x=cwc(i,k,2)*348.43_r8*pa(i,k)/ta(i,k)
                reff(i,k,2)=14.3_r8*x**0.1667_r8
                reff(i,k,2)=MAX(reff(i,k,2),4.0_r8)
                reff(i,k,2)=MIN(reff(i,k,2),20.0_r8)

             ENDIF

          ENDDO
       ENDDO

    ENDIF

    !c-----the surface (np+1) is treated as a layer filled with black clouds.
    !c     transfc is the transmittance between the surface and a pressure level.
    !c     trantcr is the clear-sky transmittance between the surface and a
    !c     pressure level.

    DO i=1,m
       sfcem(i)       =0.0_r8
       transfc(i,np+1)=1.0_r8
       trantcr(i,np+1)=1.0_r8
    ENDDO

    !c-----initialize fluxes

    DO k=1,np+1
       DO i=1,m
          flx(i,k)  = 0.0_r8
          flc(i,k)  = 0.0_r8
          dfdts(i,k)= 0.0_r8
       ENDDO
    ENDDO

    !c-----integration over spectral bands

    DO  ibn=1,10 !do 1000 ibn=1,10

       IF (ibn.EQ.10 .AND. .NOT.trace) RETURN

       !c-----if h2otbl, compute h2o (line) transmittance using table look-up.
       !c     if conbnd, compute h2o (continuum) transmittance in bands 2-7.
       !c     if co2bnd, compute co2 transmittance in band 3.
       !c     if oznbnd, compute  o3 transmittance in band 5.
       !c     if n2obnd, compute n2o transmittance in bands 6 and 7.
       !c     if ch4bnd, compute ch4 transmittance in bands 6 and 7.
       !c     if combnd, compute co2-minor transmittance in bands 4 and 5.
       !c     if f11bnd, compute cfc11 transmittance in bands 4 and 5.
       !c     if f12bnd, compute cfc12 transmittance in bands 4 and 6.
       !c     if f22bnd, compute cfc22 transmittance in bands 4 and 6.
       !c     if b10bnd, compute flux reduction due to n2o in band 10.

       h2otbl=high.AND.(ibn.EQ.1.OR.ibn.EQ.2.OR.ibn.EQ.8)
       conbnd=ibn.GE.2.AND.ibn.LE.7
       co2bnd=ibn.EQ.3
       oznbnd=ibn.EQ.5
       n2obnd=ibn.EQ.6.OR.ibn.EQ.7
       ch4bnd=ibn.EQ.6.OR.ibn.EQ.7
       combnd=ibn.EQ.4.OR.ibn.EQ.5
       f11bnd=ibn.EQ.4.OR.ibn.EQ.5
       f12bnd=ibn.EQ.4.OR.ibn.EQ.6
       f22bnd=ibn.EQ.4.OR.ibn.EQ.6
       b10bnd=ibn.EQ.10

       !c-----blayer is the spectrally integrated planck flux of the mean layer
       !c     temperature derived from eq. (3.11)
       !c     The fitting for the planck flux is valid for the range 160-345 K.

       DO k=1,np

          DO i=1,m
             tx(i)=ta(i,k)
          ENDDO
          CALL planck(ibn,m,tx,xlayer)

          DO i=1,m
             blayer(i,k)=xlayer(i)
          ENDDO

       ENDDO

       !c-----Index "0" is the layer above the top of the atmosphere.

       DO i=1,m
          blayer(i,0)=0.0_r8
       ENDDO

       !c-----Surface emission and reflectivity. See Section 9.
       !c     bs and dbs include the effect of surface emissivity.

       CALL sfcflux (ibn,m,ns,fs,tg,eg,tv,ev,rv,vege,bs,dbs,rflxs) 

       DO i=1,m
          blayer(i,np+1)=bs(i)
       ENDDO

       !c------interpolate Planck function at model levels (linear in p)

       DO k=2,np
          DO i=1,m
             blevel(i,k)=(blayer(i,k-1)*dp(i,k)+blayer(i,k)*dp(i,k-1))/ &
                  (dp(i,k-1)+dp(i,k))
          ENDDO
       ENDDO

       !c-----Extrapolate blevel(i,1) from blayer(i,2) and blayer(i,1)

       DO i=1,m
          blevel(i,1)=blayer(i,1)+(blayer(i,1)-blayer(i,2))*dp(i,1)/ &
               (dp(i,1)+dp(i,2))
       ENDDO

       !c-----If the surface air temperature tb is known, compute blevel(i,np+1)

       CALL planck(ibn,m,tb,xlayer)
       DO i=1,m
          blevel(i,np+1)=xlayer(i)
       ENDDO

       !c-----Otherwise, extrapolate blevel(np+1) from blayer(np-1) and blayer(np)

       !c      do i=1,m
       !c        blevel(i,np+1)=blayer(i,np)+(blayer(i,np)-blayer(i,np-1)) &
       !c                     *dp(i,np)/(dp(i,np)+dp(i,np-1))
       !c      enddo

       !c-----Compute cloud optical thickness following Eqs. (6.4a,b) and (6.7)
       !c     Rain optical thickness is set to 0.00307 /(gm/m**2).
       !c     It is for a specific drop size distribution provided by Q. Fu.

       !ctar      if (cldwater) then
       !ctar        do k=1,np
       !ctar        do i=1,m
       !ctar           taucl(i,k,1)=cwp(i,k,1)*(aib(1,ibn)+aib(2,ibn)/ &
       !ctar            reff(i,k,1)**aib(3,ibn))
       !ctar           taucl(i,k,2)=cwp(i,k,2)*(awb(1,ibn)+(awb(2,ibn)+ &
       !ctar            (awb(3,ibn)+awb(4,ibn)*reff(i,k,2))*reff(i,k,2)) &
       !ctar            *reff(i,k,2))
       !ctar           taucl(i,k,3)=0.00307_r8*cwp(i,k,3)
       !ctar         enddo
       !ctar        enddo
       !ctar       endif

       !c-----Compute cloud single-scattering albedo and asymmetry factor for
       !c     a mixture of ice particles and liquid drops following 
       !c     Eqs. (6.5), (6.6), (6.15) and (6.16).
       !c     Single-scattering albedo and asymmetry factor of rain are set
       !c     to 0.54 and 0.95, respectively, based on the information provided
       !c     by Prof. Qiang Fu.

       DO k=1,np
          DO i=1,m

             tcldlyr(i,k) = 1.0_r8
             tauc=taucl(i,k,1)+taucl(i,k,2)+taucl(i,k,3)

             IF (tauc.GT.0.02 .AND. fcld(i,k).GT.0.01) THEN

                w1=taucl(i,k,1)*(aiw(1,ibn)+(aiw(2,ibn)+(aiw(3,ibn) &
                     +aiw(4,ibn)*reff(i,k,1))*reff(i,k,1))*reff(i,k,1))
                w2=taucl(i,k,2)*(aww(1,ibn)+(aww(2,ibn)+(aww(3,ibn) &
                     +aww(4,ibn)*reff(i,k,2))*reff(i,k,2))*reff(i,k,2))
                w3=taucl(i,k,3)*0.54_r8
                ww=(w1+w2+w3)/tauc

                g1=w1*(aig(1,ibn)+(aig(2,ibn)+(aig(3,ibn) &
                     +aig(4,ibn)*reff(i,k,1))*reff(i,k,1))*reff(i,k,1))
                g2=w2*(awg(1,ibn)+(awg(2,ibn)+(awg(3,ibn) &
                     +awg(4,ibn)*reff(i,k,2))*reff(i,k,2))*reff(i,k,2))
                g3=w3*0.95_r8

                gg=(g1+g2+g3)/(w1+w2+w3)

                !c-----Parameterization of LW scattering following Eqs. (6.11) and (6.12). 

                ff=0.5_r8+(0.3739_r8+(0.0076_r8+0.1185_r8*gg)*gg)*gg
                tauc=(1.0_r8-ww*ff)*tauc

                !c-----compute cloud diffuse transmittance. It is approximated by using 
                !c     a diffusivity factor of 1.66.

                tcldlyr(i,k)=EXP(-1.66_r8*tauc)

             ENDIF

          ENDDO
       ENDDO

       !c-----Compute optical thickness, single-scattering albedo and asymmetry
       !c     factor for a mixture of "na" aerosol types. Eqs. (7.1)-(7.3)

       IF (aerosol) THEN

          DO k=1,np

             DO i=1,m
                taua(i,k)=0.0_r8
                ssaa(i,k)=0.0_r8
                asya(i,k)=0.0_r8
             ENDDO

             DO j=1,na
                DO i=1,m
                   taua(i,k)=taua(i,k)+taual(i,k,ibn,j)
                   w1=ssaal(i,k,ibn,j)*taual(i,k,ibn,j)
                   ssaa(i,k)=ssaa(i,k)+w1
                   asya(i,k)=asya(i,k)+asyal(i,k,ibn,j)*w1
                ENDDO
             ENDDO

             !c-----taerlyr is the aerosol diffuse transmittance

             DO i=1,m
                taerlyr(i,k)=1.0_r8

                IF (taua(i,k) .GT. 0.001_r8) THEN 
                   IF (ssaa(i,k) .GT. 0.001_r8) THEN
                      asya(i,k)=asya(i,k)/ssaa(i,k)
                      ssaa(i,k)=ssaa(i,k)/taua(i,k)

                      !c-----Parameterization of aerosol scattering following Eqs. (6.11) and (6.12). 

                      ff=0.5_r8+(0.3739_r8+(0.0076_r8+0.1185_r8*asya(i,k))*asya(i,k))*asya(i,k)
                      taua(i,k)=taua(i,k)*(1.-ssaa(i,k)*ff)

                   ENDIF
                   taerlyr(i,k)=EXP(-1.66_r8*taua(i,k))
                ENDIF

             ENDDO
          ENDDO

       ENDIF

       !c-----Compute the exponential terms (Eq. 8.21) at each layer due to
       !c     water vapor line absorption when k-distribution is used

       IF (.NOT.h2otbl .AND. .NOT.b10bnd) THEN
          CALL h2oexps(ibn,m,np,dh2o,pa,dt,xkw,aw,bw,pm,mw,h2oexp)
       ENDIF

       !c-----compute the exponential terms (Eq. 4.24) at each layer due to
       !c     water vapor continuum absorption.
       !c     ne is the number of terms used in each band to compute water 
       !c     vapor continuum transmittance (Table 9).

       ne=0
       IF (conbnd) THEN

          ne=1
          IF (ibn.EQ.3) ne=3

          CALL conexps(ibn,m,np,dcont,xke,conexp)

       ENDIF

       !c-----compute the exponential terms (Eq. 8.21) at each layer due to
       !c     co2 absorption

       IF (.NOT.high .AND. co2bnd) THEN
          CALL co2exps(m,np,dco2,pa,dt,co2exp)
       ENDIF

       !c***** for trace gases *****

       IF (trace) THEN

          !c-----compute the exponential terms at each layer due to n2o absorption

          IF (n2obnd) THEN
             CALL n2oexps(ibn,m,np,dn2o,pa,dt,n2oexp)
          ENDIF

          !c-----compute the exponential terms at each layer due to ch4 absorption

          IF (ch4bnd) THEN
             CALL ch4exps(ibn,m,np,dch4,pa,dt,ch4exp)
          ENDIF

          !c-----Compute the exponential terms due to co2 minor absorption

          IF (combnd) THEN
             CALL comexps(ibn,m,np,dco2,dt,comexp)
          ENDIF

          !c-----Compute the exponential terms due to cfc11 absorption.
          !c     The values of the parameters are given in Table 7.

          IF (f11bnd) THEN
             a1  = 1.26610e-3_r8
             b1  = 3.55940e-6_r8
             fk1 = 1.89736e+1_r8
             a2  = 8.19370e-4_r8
             b2  = 4.67810e-6_r8
             fk2 = 1.01487e+1_r8
             CALL cfcexps(ibn,m,np,a1,b1,fk1,a2,b2,fk2,df11,dt,f11exp)
          ENDIF

          !c-----Compute the exponential terms due to cfc12 absorption.

          IF (f12bnd) THEN
             a1  = 8.77370e-4_r8
             b1  =-5.88440e-6_r8
             fk1 = 1.58104e+1_r8
             a2  = 8.62000e-4_r8
             b2  =-4.22500e-6_r8
             fk2 = 3.70107e+1_r8
             CALL cfcexps(ibn,m,np,a1,b1,fk1,a2,b2,fk2,df12,dt,f12exp)
          ENDIF

          !c-----Compute the exponential terms due to cfc22 absorption.

          IF (f22bnd) THEN
             a1  = 9.65130e-4_r8
             b1  = 1.31280e-5_r8
             fk1 = 6.18536e+0_r8
             a2  =-3.00010e-5_r8 
             b2  = 5.25010e-7_r8
             fk2 = 3.27912e+1_r8
             CALL cfcexps(ibn,m,np,a1,b1,fk1,a2,b2,fk2,df22,dt,f22exp)
          ENDIF

          !c-----Compute the exponential terms at each layer in band 10 due to
          !c     h2o line and continuum, co2, and n2o absorption

          IF (b10bnd) THEN
             CALL b10exps(m,np,dh2o,dcont,dco2,dn2o,pa,dt &
                  ,h2oexp,conexp,co2exp,n2oexp)
          ENDIF

       ENDIF

       !c-----blayer(i,np+1) includes the effect of surface emissivity.

       DO i=1,m
          bd(i,0)=0.0_r8
          bu(i,np+1)=blayer(i,np+1)
       ENDDO

       !c-----do-loop 1500 is for computing upward (bu) and downward (bd)
       !c     emission of a layer following Eqs. (8.17), (8.18), (8.19).
       !c     Here, trant(i,k2) is the transmittance of the layer k2-1.

       DO  k2=2,np+1 !do 1500 k2=2,np+1

          !c-----for h2o line transmission

          IF (.NOT. h2otbl) THEN
             DO ik=1,6
                DO i=1,m
                   th2o(i,ik)=1.0_r8
                ENDDO
             ENDDO
          ENDIF

          !c-----for h2o continuum transmission

          DO iq=1,3
             DO i=1,m
                tcon(i,iq)=1.0_r8
             ENDDO
          ENDDO

          !c-----for co2 transmission using k-distribution method.
          !c     band 3 is divided into 3 sub-bands, but sub-bands 3a and 3c
          !c     are combined in computing the co2 transmittance.

          IF (.NOT.high .AND. co2bnd) THEN
             DO isb=1,2
                DO ik=1,6
                   DO i=1,m
                      tco2(i,ik,isb)=1.0_r8
                   ENDDO
                ENDDO
             ENDDO
          ENDIF


          DO i=1,m
             x1(i)=0.0_r8
             x2(i)=0.0_r8
             x3(i)=0.0_r8
             trant(i,k2)=1.0_r8
          ENDDO

          IF (h2otbl) THEN

             !c-----Compute water vapor transmittance using table look-up.
             !c     The following values are taken from Table 8.

             w1=-8.0_r8
             p1=-2.0_r8
             dwe=0.3_r8
             dpe=0.2_r8

             IF (ibn.EQ.1) THEN
                CALL tablup(k2,m,np,nx,nh,dh2o,pa,dt,x1,x2,x3, &
                     w1,p1,dwe,dpe,h11,h12,h13,trant)

             ENDIF
             IF (ibn.EQ.2) THEN
                CALL tablup(k2,m,np,nx,nh,dh2o,pa,dt,x1,x2,x3, &
                     w1,p1,dwe,dpe,h21,h22,h23,trant)

             ENDIF
             IF (ibn.EQ.8) THEN
                CALL tablup(k2,m,np,nx,nh,dh2o,pa,dt,x1,x2,x3, &
                     w1,p1,dwe,dpe,h81,h82,h83,trant)
             ENDIF

             !c-----for water vapor continuum absorption

             IF (conbnd) THEN
                DO i=1,m
                   tcon(i,1)=tcon(i,1)*conexp(i,k2-1,1)
                   trant(i,k2)=trant(i,k2)*tcon(i,1)
                ENDDO
             ENDIF

          ELSE

             !c-----compute water vapor transmittance using k-distribution

             IF (.NOT.b10bnd) THEN
                CALL h2okdis(ibn,m,np,k2-1,fkw,gkw,ne,h2oexp,conexp, &
                     th2o,tcon,trant)

             ENDIF

          ENDIF

          IF (co2bnd) THEN

             IF (high) THEN

                !c-----Compute co2 transmittance using table look-up method.
                !c     The following values are taken from Table 8.

                w1=-4.0_r8
                p1=-2.0_r8
                dwe=0.3_r8
                dpe=0.2_r8
                CALL tablup(k2,m,np,nx,nc,dco2,pa,dt,x1,x2,x3, &
                     w1,p1,dwe,dpe,c1,c2,c3,trant)
             ELSE

                !c-----compute co2 transmittance using k-distribution method
                CALL co2kdis(m,np,k2-1,co2exp,tco2,trant)

             ENDIF

          ENDIF

          !c-----Always use table look-up to compute o3 transmittance.
          !c     The following values are taken from Table 8.

          IF (oznbnd) THEN
             w1=-6.0_r8
             p1=-2.0_r8
             dwe=0.3_r8
             dpe=0.2_r8
             CALL tablup(k2,m,np,nx,no,do3,pa,dt,x1,x2,x3, &
                  w1,p1,dwe,dpe,o1,o2,o3,trant)
          ENDIF

          !c-----include aerosol effect

          IF (aerosol) THEN
             DO i=1,m
                trant(i,k2)=trant(i,k2)*taerlyr(i,k2-1)
             ENDDO
          ENDIF

          !c-----Compute upward and downward emission of the layer k2-1

          DO i=1,m

             xx=(blayer(i,k2-1)-blevel(i,k2-1))* &
                  (blayer(i,k2-1)-blevel(i,k2))

             IF (xx.GT.0.0_r8) THEN

                !c-----If xx>0, there is a local temperature minimum or maximum.
                !c     Computations of bd and bu follow Eq. (8.20).

                bd(i,k2-1)=.5_r8*blayer(i,k2-1)+.25_r8*(blevel(i,k2-1)+blevel(i,k2))
                bu(i,k2-1)=bd(i,k2-1)

             ELSE

                !c-----Computations of bd and bu following Eqs.(8.17) and (8.18).
                !c     The effect of clouds on the transmission of a layer is taken
                !c     into account, following Eq. (8.19).

                xx=(fcld(i,k2-1)*tcldlyr(i,k2-1)+(1.-fcld(i,k2-1))) &
                     *trant(i,k2)

                yy=MIN(0.9999_r8,xx)
                yy=MAX(0.00001_r8,yy)
                xx=(blevel(i,k2-1)-blevel(i,k2))/log(yy)
                bd(i,k2-1)=(blevel(i,k2)-blevel(i,k2-1)*yy)/(1.0_r8-yy)-xx
                bu(i,k2-1)=(blevel(i,k2-1)+blevel(i,k2))-bd(i,k2-1)

             ENDIF

          ENDDO

       END DO !1500 continue

       !c-----initialize fluxes

       DO k=1,np+1
          DO i=1,m
             flxu(i,k) = 0.0_r8
             flxd(i,k) = 0.0_r8
             flcu(i,k) = 0.0_r8
             flcd(i,k) = 0.0_r8
          ENDDO
       ENDDO

       !c-----

       DO  k1=1,np ! do 2000 k1=1,np

          !c-----initialization
          !c
          !c     it, im, and ib are the numbers of cloudy layers in the high,
          !c     middle, and low cloud groups between levels k1 and k2.
          !c     cldlw, cldmd, and cldhi are the equivalent black-cloud fractions
          !c     of low, middle, and high troposphere.
          !c     tranal is the aerosol transmission function

          DO i=1,m
             it(i) = 0
             im(i) = 0
             ib(i) = 0
             cldlw(i) = 0.0_r8
             cldmd(i) = 0.0_r8
             cldhi(i) = 0.0_r8
             tranal(i)= 1.0_r8
          ENDDO

          !c-----for h2o line transmission

          IF (.NOT. h2otbl) THEN
             DO ik=1,6
                DO i=1,m
                   th2o(i,ik)=1.0_r8
                ENDDO
             ENDDO
          ENDIF

          !c-----for h2o continuum transmission

          DO iq=1,3
             DO i=1,m
                tcon(i,iq)=1.0_r8
             ENDDO
          ENDDO

          !c-----for co2 transmission using k-distribution method.
          !c     band 3 is divided into 3 sub-bands, but sub-bands 3a and 3c
          !c     are combined in computing the co2 transmittance.

          IF (.NOT.high .AND. co2bnd) THEN
             DO isb=1,2
                DO ik=1,6
                   DO i=1,m
                      tco2(i,ik,isb)=1.0_r8
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          !c***** for trace gases *****

          IF (trace) THEN

             !c-----for n2o transmission using k-distribution method.

             IF (n2obnd) THEN
                DO ik=1,4
                   DO i=1,m
                      tn2o(i,ik)=1.0_r8
                   ENDDO
                ENDDO
             ENDIF

             !c-----for ch4 transmission using k-distribution method.

             IF (ch4bnd) THEN
                DO ik=1,4
                   DO i=1,m
                      tch4(i,ik)=1.0_r8
                   ENDDO
                ENDDO
             ENDIF

             !c-----for co2-minor transmission using k-distribution method.

             IF (combnd) THEN
                DO ik=1,6
                   DO i=1,m
                      tcom(i,ik)=1.0_r8
                   ENDDO
                ENDDO
             ENDIF

             !c-----for cfc-11 transmission using k-distribution method.

             IF (f11bnd) THEN
                DO i=1,m
                   tf11(i)=1.0_r8
                ENDDO
             ENDIF

             !c-----for cfc-12 transmission using k-distribution method.

             IF (f12bnd) THEN
                DO i=1,m
                   tf12(i)=1.0_r8
                ENDDO
             ENDIF

             !c-----for cfc-22 transmission when using k-distribution method.

             IF (f22bnd) THEN
                DO i=1,m
                   tf22(i)=1.0_r8
                ENDDO
             ENDIF

             !c-----for the transmission in band 10 using k-distribution method.

             IF (b10bnd) THEN
                DO ik=1,5
                   DO i=1,m
                      th2o(i,ik)=1.0_r8
                   ENDDO
                ENDDO

                DO ik=1,6
                   DO i=1,m
                      tco2(i,ik,1)=1.0_r8
                   ENDDO
                ENDDO

                DO i=1,m
                   tcon(i,1)=1.0_r8
                ENDDO

                DO ik=1,2
                   DO i=1,m
                      tn2o(i,ik)=1.0_r8
                   ENDDO
                ENDDO
             ENDIF

          ENDIF

          !c***** end trace gases *****

          DO i=1,m
             x1(i)=0.0_r8
             x2(i)=0.0_r8
             x3(i)=0.0_r8
          ENDDO

          !c-----trant is the total transmittance between levels k1 and k2.
          !c     fclr is the clear line-of-sight  between levels k1 and k2.

          DO k=1,np+1
             DO i=1,m
                trant(i,k)=1.0_r8
                fclr(i,k) =1.0_r8
             ENDDO
          ENDDO

          !c-----do-loop 3000 are for computing (a) transmittance, trant(i,k2),
          !c     and (b) clear line-of-sight, fclr(i,k2), between levels k1 and k2.

          DO  k2=k1+1,np+1 !do 3000 k2=k1+1,np+1

             IF (h2otbl) THEN

                !c-----Compute water vapor transmittance using table look-up.
                !c     The following values are taken from Table 8.

                w1=-8.0_r8
                p1=-2.0_r8
                dwe=0.3_r8
                dpe=0.2_r8

                IF (ibn.EQ.1) THEN
                   CALL tablup(k2,m,np,nx,nh,dh2o,pa,dt,x1,x2,x3, &
                        w1,p1,dwe,dpe,h11,h12,h13,trant)

                ENDIF
                IF (ibn.EQ.2) THEN
                   CALL tablup(k2,m,np,nx,nh,dh2o,pa,dt,x1,x2,x3, &
                        w1,p1,dwe,dpe,h21,h22,h23,trant)

                ENDIF
                IF (ibn.EQ.8) THEN
                   CALL tablup(k2,m,np,nx,nh,dh2o,pa,dt,x1,x2,x3, &
                        w1,p1,dwe,dpe,h81,h82,h83,trant)
                ENDIF

                IF (conbnd) THEN
                   DO i=1,m
                      tcon(i,1)=tcon(i,1)*conexp(i,k2-1,1)
                      trant(i,k2)=trant(i,k2)*tcon(i,1)
                   ENDDO
                ENDIF

             ELSE

                !c-----compute water vapor transmittance using k-distribution

                IF (.NOT.b10bnd) THEN
                   CALL h2okdis(ibn,m,np,k2-1,fkw,gkw,ne,h2oexp,conexp, &
                        th2o,tcon,trant)
                ENDIF

             ENDIF

             IF (co2bnd) THEN

                IF (high) THEN

                   !c-----Compute co2 transmittance using table look-up method.
                   !c     The following values are taken from Table 8.

                   w1=-4.0_r8
                   p1=-2.0_r8
                   dwe=0.3_r8
                   dpe=0.2_r8
                   CALL tablup(k2,m,np,nx,nc,dco2,pa,dt,x1,x2,x3, &
                        w1,p1,dwe,dpe,c1,c2,c3,trant)
                ELSE

                   !c-----compute co2 transmittance using k-distribution method
                   CALL co2kdis(m,np,k2-1,co2exp,tco2,trant)

                ENDIF

             ENDIF

             !c-----Always use table look-up to compute o3 transmittance.
             !c     The following values are taken from Table 8.

             IF (oznbnd) THEN
                w1=-6.0_r8
                p1=-2.0_r8
                dwe=0.3_r8
                dpe=0.2_r8
                CALL tablup(k2,m,np,nx,no,do3,pa,dt,x1,x2,x3, &
                     w1,p1,dwe,dpe,o1,o2,o3,trant)
             ENDIF

             !c***** for trace gases *****

             IF (trace) THEN

                !c-----compute n2o transmittance using k-distribution method

                IF (n2obnd) THEN
                   CALL n2okdis(ibn,m,np,k2-1,n2oexp,tn2o,trant)
                ENDIF

                !c-----compute ch4 transmittance using k-distribution method

                IF (ch4bnd) THEN
                   CALL ch4kdis(ibn,m,np,k2-1,ch4exp,tch4,trant)
                ENDIF

                !c-----compute co2-minor transmittance using k-distribution method

                IF (combnd) THEN
                   CALL comkdis(ibn,m,np,k2-1,comexp,tcom,trant)
                ENDIF

                !c-----compute cfc11 transmittance using k-distribution method

                IF (f11bnd) THEN
                   CALL cfckdis(m,np,k2-1,f11exp,tf11,trant)
                ENDIF

                !c-----compute cfc12 transmittance using k-distribution method

                IF (f12bnd) THEN
                   CALL cfckdis(m,np,k2-1,f12exp,tf12,trant)
                ENDIF

                !c-----compute cfc22 transmittance using k-distribution method

                IF (f22bnd) THEN
                   CALL cfckdis(m,np,k2-1,f22exp,tf22,trant)
                ENDIF

                !c-----Compute transmittance in band 10 using k-distribution method.
                !c     For band 10, trant is the change in transmittance due to n2o 
                !c     absorption.

                IF (b10bnd) THEN
                   CALL b10kdis(m,np,k2-1,h2oexp,conexp,co2exp,n2oexp &
                        ,th2o,tcon,tco2,tn2o,trant)

                ENDIF

             ENDIF

             !c*****   end trace gases  *****

             !c-----include aerosol effect

             IF (aerosol) THEN
                DO i=1,m
                   tranal(i)=tranal(i)*taerlyr(i,k2-1)
                   trant(i,k2)=trant(i,k2) *tranal(i)
                ENDDO
             ENDIF

             !c***** cloud overlapping *****

             IF (.NOT. overcast) THEN
                CALL cldovlp (m,np,k2,ict,icb,it,im,ib,itx,imx,ibx, &
                     cldhi,cldmd,cldlw,fcld,tcldlyr,fclr)

             ELSE

                DO i=1,m
                   fclr(i,k2)=fclr(i,k2)*tcldlyr(i,k2-1)
                ENDDO

             ENDIF

          END DO !3000 continue

          !c-----do-loop 4000 is for computing upward and downward fluxes
          !c     for each spectral band
          !c     flcu, flcd: clear-sky upward and downward fluxes
          !c     flxu, flxd: all-sky   upward and downward fluxes

          DO  k2=k1+1,np+1 ! do 4000 k2=k1+1,np+1

             IF (k2.EQ.k1+1 .AND. ibn .NE. 10) THEN

                !c-----The first terms on the rhs of Eqs. (8.15) and (8.16)

                DO i=1,m
                   flcu(i,k1)=flcu(i,k1)-bu(i,k1)
                   flcd(i,k2)=flcd(i,k2)+bd(i,k1)
                   flxu(i,k1)=flxu(i,k1)-bu(i,k1)
                   flxd(i,k2)=flxd(i,k2)+bd(i,k1)
                ENDDO

             ENDIF

             !c-----The summation terms on the rhs of Eqs. (8.15) and (8.16).
             !c     Also see Eqs. (5.4) and (5.5) for Band 10.

             DO i=1,m
                xx=trant(i,k2)*(bu(i,k2-1)-bu(i,k2))
                flcu(i,k1) =flcu(i,k1)+xx
                flxu(i,k1) =flxu(i,k1)+xx*fclr(i,k2)
                xx=trant(i,k2)*(bd(i,k1-1)-bd(i,k1))
                flcd(i,k2) =flcd(i,k2)+xx
                flxd(i,k2) =flxd(i,k2)+xx*fclr(i,k2)
             ENDDO

          END DO! 4000 continue

          !c-----Here, fclr and trant are, respectively, the clear line-of-sight 
          !c     and the transmittance between k1 and the surface.

          DO i=1,m
             trantcr(i,k1) =trant(i,np+1)
             transfc(i,k1) =trant(i,np+1)*fclr(i,np+1)
          ENDDO

          !c-----compute the partial derivative of fluxes with respect to
          !c     surface temperature (Eq. 3.12). 
          !c     Note: upward flux is negative, and so is dfdts.

          DO i=1,m
             dfdts(i,k1) =dfdts(i,k1)-dbs(i)*transfc(i,k1)
          ENDDO

       END DO! 2000 continue

       IF (.NOT. b10bnd) THEN

          !c-----For surface emission.
          !c     Note: blayer(i,np+1) and dbs include the surface emissivity effect.
          !c     Both dfdts and sfcem are negative quantities.

          DO i=1,m
             flcu(i,np+1)=-blayer(i,np+1)
             flxu(i,np+1)=-blayer(i,np+1)
             sfcem(i)=sfcem(i)-blayer(i,np+1)
             dfdts(i,np+1)=dfdts(i,np+1)-dbs(i)
          ENDDO


          !c-----Add the flux reflected by the surface. (Second term on the
          !c     rhs of Eq. 8.16)

          DO k=1,np+1
             DO i=1,m
                flcu(i,k)=flcu(i,k)- &
                     flcd(i,np+1)*trantcr(i,k)*rflxs(i)
                flxu(i,k)=flxu(i,k)-   &                      
                     flxd(i,np+1)*transfc(i,k)*rflxs(i)   
             ENDDO
          ENDDO

       ENDIF

       !c-----Summation of fluxes over spectral bands

       DO k=1,np+1
          DO i=1,m
             flc(i,k)=flc(i,k)+flcd(i,k)+flcu(i,k)
             flx(i,k)=flx(i,k)+flxd(i,k)+flxu(i,k)
          ENDDO
       ENDDO

    END DO! 1000 continue


  END SUBROUTINE irrad

  !c***********************************************************************
  SUBROUTINE planck(ibn,m,t,xlayer)
    !c***********************************************************************
    !c
    !c-----Compute spectrally integrated Planck flux
    !c
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    INTEGER, INTENT(IN) :: ibn,m

    REAL(KIND=r8), INTENT(IN), DIMENSION(m)  :: t

    REAL(KIND=r8), INTENT(OUT), DIMENSION(m)  ::  xlayer

    REAL(KIND=r8), PARAMETER, DIMENSION(6,10) :: cb = RESHAPE(  &
         SHAPE = (/ 6, 10 /), SOURCE = (/ &
         5.3443e+0_r8,  -2.0617e-1_r8,   2.5333e-3_r8, &
         -6.8633e-6_r8,   1.0115e-8_r8,  -6.2672e-12_r8, &
         2.7148e+1_r8,  -5.4038e-1_r8,   2.9501e-3_r8, &
         2.7228e-7_r8,  -9.3384e-9_r8,   9.9677e-12_r8, &
         -3.4860e+1_r8,   1.1132e+0_r8,  -1.3006e-2_r8, &
         6.4955e-5_r8,  -1.1815e-7_r8,   8.0424e-11_r8, &
         -6.0513e+1_r8,   1.4087e+0_r8,  -1.2077e-2_r8, &
         4.4050e-5_r8,  -5.6735e-8_r8,   2.5660e-11_r8, &
         -2.6689e+1_r8,   5.2828e-1_r8,  -3.4453e-3_r8, &
         6.0715e-6_r8,   1.2523e-8_r8,  -2.1550e-11_r8, &
         -6.7274e+0_r8,   4.2256e-2_r8,   1.0441e-3_r8, &
         -1.2917e-5_r8,   4.7396e-8_r8,  -4.4855e-11_r8, &
         1.8786e+1_r8,  -5.8359e-1_r8,   6.9674e-3_r8, &
         -3.9391e-5_r8,   1.0120e-7_r8,  -8.2301e-11_r8, &
         1.0344e+2_r8,  -2.5134e+0_r8,   2.3748e-2_r8, &
         -1.0692e-4_r8,   2.1841e-7_r8,  -1.3704e-10_r8, &
         -1.0482e+1_r8,   3.8213e-1_r8,  -5.2267e-3_r8, &
         3.4412e-5_r8,  -1.1075e-7_r8,   1.4092e-10_r8, &
         1.6769e+0_r8,   6.5397e-2_r8,  -1.8125e-3_r8, &
         1.2912e-5_r8,  -2.6715e-8_r8,   1.9792e-11_r8 /) )
    INTEGER  :: i

    !ctar      integer ibn                   ! spectral band index
    !ctar      integer m                     ! no of points
    !ctar      real t(m)                     ! temperature (K)
    !ctar      real xlayer(m)                ! planck flux (w/m2)
    !ctar      real cb(6,10)
    !ctar      integer i

    !c-----the following coefficients are given in Table 2 for computing  
    !c     spectrally integrated planck fluxes using Eq. (3.11)

    !ctar       data cb/ &
    !ctar           5.3443e+0,  -2.0617e-1,   2.5333e-3, &
    !ctar          -6.8633e-6,   1.0115e-8,  -6.2672e-12, &
    !ctar           2.7148e+1,  -5.4038e-1,   2.9501e-3, &
    !ctar           2.7228e-7,  -9.3384e-9,   9.9677e-12, &
    !ctar          -3.4860e+1,   1.1132e+0,  -1.3006e-2, &
    !ctar           6.4955e-5,  -1.1815e-7,   8.0424e-11, &
    !ctar          -6.0513e+1,   1.4087e+0,  -1.2077e-2, &
    !ctar           4.4050e-5,  -5.6735e-8,   2.5660e-11, &
    !ctar         -2.6689e+1,   5.2828e-1,  -3.4453e-3, &
    !ctar           6.0715e-6,   1.2523e-8,  -2.1550e-11, &
    !ctar          -6.7274e+0,   4.2256e-2,   1.0441e-3, &
    !ctar         -1.2917e-5,   4.7396e-8,  -4.4855e-11, &
    !ctar           1.8786e+1,  -5.8359e-1,   6.9674e-3, &
    !ctar          -3.9391e-5,   1.0120e-7,  -8.2301e-11, &
    !ctar           1.0344e+2,  -2.5134e+0,   2.3748e-2, &
    !ctar          -1.0692e-4,   2.1841e-7,  -1.3704e-10, &
    !ctar          -1.0482e+1,   3.8213e-1,  -5.2267e-3, &
    !ctar           3.4412e-5,  -1.1075e-7,   1.4092e-10, &
    !ctar           1.6769e+0,   6.5397e-2,  -1.8125e-3, &
    !ctar          1.2912e-5,  -2.6715e-8,   1.9792e-11/
    !c
    DO i=1,m
       xlayer(i)=t(i)*(t(i)*(t(i)*(t(i)*(t(i)*cb(6,ibn)+cb(5,ibn)) &
            +cb(4,ibn))+cb(3,ibn))+cb(2,ibn))+cb(1,ibn)
    ENDDO


  END SUBROUTINE planck

  !c***********************************************************************
  SUBROUTINE plancd(ibn,m,t,dbdt) 
    !c***********************************************************************
    !c
    !c-----Compute the derivative of Planck flux wrt temperature
    !c
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    INTEGER, INTENT(IN) :: ibn,m

    REAL(KIND=r8), INTENT(IN), DIMENSION(m)  :: t
    REAL(KIND=r8), INTENT(OUT), DIMENSION(m)  :: dbdt 

    REAL(KIND=r8), PARAMETER, DIMENSION(5,10) :: dcb =  RESHAPE(  &
         SHAPE = (/ 5, 10 /), SOURCE = (/ &
         -2.0617E-01_r8, 5.0666E-03_r8,-2.0590E-05_r8, 4.0460E-08_r8,-3.1336E-11_r8, &
         -5.4038E-01_r8, 5.9002E-03_r8, 8.1684E-07_r8,-3.7354E-08_r8, 4.9839E-11_r8, &
         1.1132E+00_r8,-2.6012E-02_r8, 1.9486E-04_r8,-4.7260E-07_r8, 4.0212E-10_r8, &
         1.4087E+00_r8,-2.4154E-02_r8, 1.3215E-04_r8,-2.2694E-07_r8, 1.2830E-10_r8, &
         5.2828E-01_r8,-6.8906E-03_r8, 1.8215E-05_r8, 5.0092E-08_r8,-1.0775E-10_r8, &
         4.2256E-02_r8, 2.0882E-03_r8,-3.8751E-05_r8, 1.8958E-07_r8,-2.2428E-10_r8, &
         -5.8359E-01_r8, 1.3935E-02_r8,-1.1817E-04_r8, 4.0480E-07_r8,-4.1150E-10_r8, &
         -2.5134E+00_r8, 4.7496E-02_r8,-3.2076E-04_r8, 8.7364E-07_r8,-6.8520E-10_r8, &
         3.8213E-01_r8,-1.0453E-02_r8, 1.0324E-04_r8,-4.4300E-07_r8, 7.0460E-10_r8, &
         6.5397E-02_r8,-3.6250E-03_r8, 3.8736E-05_r8,-1.0686E-07_r8, 9.8960E-11_r8 /) )
    INTEGER :: i

    !ctar      integer ibn               ! spectral band index
    !ctar      integer m                 ! no of points
    !ctar      real t(m)                 ! temperature (K)
    !ctar      real dbdt(m)              ! derivative of Planck flux wrt temperature
    !ctar      real dcb(5,10)
    !ctar      integer i

    !c-----Coefficients for computing the derivative of Planck function
    !c     with respect to temperature (Eq. 3.12).
    !c     dcb(1)=1*cb(2), dcb(2)=2*cb(3), dcb(3)=3*cb(4) ...  etc

    !ctar       data dcb/ &
    !ctar       -2.0617E-01, 5.0666E-03,-2.0590E-05, 4.0460E-08,-3.1336E-11, &
    !ctar      -5.4038E-01, 5.9002E-03, 8.1684E-07,-3.7354E-08, 4.9839E-11, &
    !ctar        1.1132E+00,-2.6012E-02, 1.9486E-04,-4.7260E-07, 4.0212E-10, &
    !ctar        1.4087E+00,-2.4154E-02, 1.3215E-04,-2.2694E-07, 1.2830E-10, &
    !ctar        5.2828E-01,-6.8906E-03, 1.8215E-05, 5.0092E-08,-1.0775E-10, &
    !ctar        4.2256E-02, 2.0882E-03,-3.8751E-05, 1.8958E-07,-2.2428E-10, &
    !ctar       -5.8359E-01, 1.3935E-02,-1.1817E-04, 4.0480E-07,-4.1150E-10, &
    !ctar       -2.5134E+00, 4.7496E-02,-3.2076E-04, 8.7364E-07,-6.8520E-10, &
    !ctar        3.8213E-01,-1.0453E-02, 1.0324E-04,-4.4300E-07, 7.0460E-10, &
    !ctar       6.5397E-02,-3.6250E-03, 3.8736E-05,-1.0686E-07, 9.8960E-11/
    !c
    DO i=1,m
       dbdt(i)=t(i)*(t(i)*(t(i)*(t(i)*dcb(5,ibn)+dcb(4,ibn)) &
            +dcb(3,ibn))+dcb(2,ibn))+dcb(1,ibn)
    ENDDO


  END SUBROUTINE plancd

  !c**********************************************************************
  SUBROUTINE h2oexps(ib,m,np,dh2o,pa,dt,xkw,aw,bw,pm,mw,h2oexp)
    !c**********************************************************************
    !c   Compute exponentials for water vapor line absorption
    !c   in individual layers using Eqs. (8.21) and (8.22).
    !c
    !c---- input parameters
    !c  spectral band (ib)
    !c  number of grid intervals (m)
    !c  number of layers (np)
    !c  layer water vapor amount for line absorption (dh2o) 
    !c  layer pressure (pa)
    !c  layer temperature minus 250K (dt)
    !c  absorption coefficients for the first k-distribution
    !c     function due to h2o line absorption (xkw)
    !c  coefficients for the temperature and pressure scaling (aw,bw,pm)
    !c  ratios between neighboring absorption coefficients for
    !c     h2o line absorption (mw)
    !c
    !c---- output parameters
    !c  6 exponentials for each layer  (h2oexp)
    !c**********************************************************************
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !c---- input parameters ------

    INTEGER, INTENT(IN) :: ib,m,np
    INTEGER :: i,k,ik

    !ctar      integer ib,m,np,i,k,ik


    !ctar       real dh2o(m,np),pa(m,np),dt(m,np)

    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np):: dh2o, pa, dt

    !c---- output parameters -----

    !ctar       real h2oexp(m,np,6)

    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np,6):: h2oexp

    !c---- static data -----

    !ctar       integer mw(9)
    !ctar       real xkw(9),aw(9),bw(9),pm(9)

    INTEGER, INTENT(IN), DIMENSION(9):: mw
    REAL(KIND=r8), INTENT(IN) , DIMENSION(9):: xkw,aw,bw,pm


    !c---- temporary arrays -----

    !ctar       real xh

    REAL(KIND=r8) ::  xh

    !c**********************************************************************
    !c    note that the 3 sub-bands in band 3 use the same set of xkw, aw,
    !c    and bw,  therefore, h2oexp for these sub-bands are identical.
    !c**********************************************************************

    DO k=1,np
       DO i=1,m

          !c-----xh is the scaled water vapor amount for line absorption
          !c     computed from Eq. (4.4).

          xh = dh2o(i,k)*(pa(i,k)/500.0_r8)**pm(ib) &
               * ( 1.0_r8+(aw(ib)+bw(ib)* dt(i,k))*dt(i,k) )

          !c-----h2oexp is the water vapor transmittance of the layer k
          !c     due to line absorption

          h2oexp(i,k,1) = EXP(-xh*xkw(ib))

       ENDDO
    ENDDO

    !c-----compute transmittances from Eq. (8.22)

    DO ik=2,6

       IF (mw(ib).EQ.6) THEN

          DO k=1,np
             DO i=1,m
                xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
                h2oexp(i,k,ik) = xh*xh*xh
             ENDDO
          ENDDO

       ELSEIF (mw(ib).EQ.8) THEN

          DO k=1,np
             DO i=1,m
                xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
                xh = xh*xh
                h2oexp(i,k,ik) = xh*xh
             ENDDO
          ENDDO

       ELSEIF (mw(ib).EQ.9) THEN

          DO k=1,np
             DO i=1,m
                xh=h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
                h2oexp(i,k,ik) = xh*xh*xh
             ENDDO
          ENDDO

       ELSE

          DO k=1,np
             DO i=1,m
                xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
                xh = xh*xh
                xh = xh*xh
                h2oexp(i,k,ik) = xh*xh
             ENDDO
          ENDDO

       ENDIF
    ENDDO

  END SUBROUTINE h2oexps

  !c**********************************************************************
  SUBROUTINE conexps(ib,m,np,dcont,xke,conexp)
    !c**********************************************************************
    !c   compute exponentials for continuum absorption in individual layers.
    !c
    !c---- input parameters
    !c  spectral band (ib)
    !c  number of grid intervals (m)
    !c  number of layers (np)
    !c  layer scaled water vapor amount for continuum absorption (dcont) 
    !c  absorption coefficients for the first k-distribution function
    !c     due to water vapor continuum absorption (xke)
    !c
    !c---- output parameters
    !c  1 or 3 exponentials for each layer (conexp)
    !c**********************************************************************
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !ctar      integer ib,m,np,i,k

    INTEGER, INTENT(IN) :: ib,m,np
    INTEGER :: i,k

    !c---- input parameters ------

    !ctar      real dcont(m,np)

    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np):: dcont

    !c---- updated parameters -----

    !ctar       real conexp(m,np,3)

    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np,3):: conexp

    !c---- static data -----

    !ctar       real xke(9)

    REAL(KIND=r8), INTENT(IN), DIMENSION(9):: xke

    !c****************************************************************

    DO k=1,np
       DO i=1,m
          conexp(i,k,1) = EXP(-dcont(i,k)*xke(ib))
       ENDDO
    ENDDO

    IF (ib .EQ. 3) THEN

       !c-----The absorption coefficients for sub-bands 3b and 3a are, respectively,
       !c     two and four times the absorption coefficient for sub-band 3c (Table 9).
       !c     Note that conexp(i,k,3) is for sub-band 3a. 

       DO k=1,np
          DO i=1,m
             conexp(i,k,2) = conexp(i,k,1) *conexp(i,k,1)
             conexp(i,k,3) = conexp(i,k,2) *conexp(i,k,2)
          ENDDO
       ENDDO

    ENDIF

  END SUBROUTINE conexps

  !c**********************************************************************
  SUBROUTINE co2exps(m,np,dco2,pa,dt,co2exp)
    !c**********************************************************************
    !c   Compute co2 exponentials for individual layers.
    !c
    !c---- input parameters
    !c  number of grid intervals (m)
    !c  number of layers (np)
    !c  layer co2 amount (dco2)
    !c  layer pressure (pa)
    !c  layer temperature minus 250K (dt)
    !c
    !c---- output parameters
    !c  6 exponentials for each layer (co2exp)
    !c**********************************************************************
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !ctar      integer m,np,i,k

    INTEGER, INTENT(IN) :: m,np
    INTEGER :: i,k

    !c---- input parameters -----

    !ctar      real dco2(m,np),pa(m,np),dt(m,np)

    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np):: dco2,pa,dt

    !c---- output parameters -----

    !ctar      real co2exp(m,np,6,2)

    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np,6,2):: co2exp

    !c---- temporary arrays -----

    !ctar      real xc

    REAL(KIND=r8):: xc


    !c**********************************************************************

    DO k=1,np
       DO i=1,m

          !c-----The scaling parameters are given in Table 3, and values of
          !c     the absorption coefficient are given in Table 10.

          !c     Scaled co2 amount for band-wings (sub-bands 3a and 3c)

          xc = dco2(i,k)*(pa(i,k)/300.0_r8)**0.5_r8 &
               *(1.0_r8+(0.0182_r8+1.07e-4_r8*dt(i,k))*dt(i,k))

          !c-----six exponentials by powers of 8 (See Eqs. 8.21, 8.22 and Table 10).

          co2exp(i,k,1,1)=EXP(-xc*2.656e-5_r8)

          xc=co2exp(i,k,1,1)*co2exp(i,k,1,1)
          xc=xc*xc
          co2exp(i,k,2,1)=xc*xc

          xc=co2exp(i,k,2,1)*co2exp(i,k,2,1)
          xc=xc*xc
          co2exp(i,k,3,1)=xc*xc

          xc=co2exp(i,k,3,1)*co2exp(i,k,3,1)
          xc=xc*xc
          co2exp(i,k,4,1)=xc*xc

          xc=co2exp(i,k,4,1)*co2exp(i,k,4,1)
          xc=xc*xc
          co2exp(i,k,5,1)=xc*xc

          xc=co2exp(i,k,5,1)*co2exp(i,k,5,1)
          xc=xc*xc
          co2exp(i,k,6,1)=xc*xc

          !c-----For band-center region (sub-band 3b)

          xc = dco2(i,k)*(pa(i,k)/30.0_r8)**0.85_r8 &
               *(1.0_r8+(0.0042_r8+2.00e-5_r8*dt(i,k))*dt(i,k))

          co2exp(i,k,1,2)=EXP(-xc*2.656e-3_r8)

          xc=co2exp(i,k,1,2)*co2exp(i,k,1,2)
          xc=xc*xc
          co2exp(i,k,2,2)=xc*xc

          xc=co2exp(i,k,2,2)*co2exp(i,k,2,2)
          xc=xc*xc
          co2exp(i,k,3,2)=xc*xc

          xc=co2exp(i,k,3,2)*co2exp(i,k,3,2)
          xc=xc*xc
          co2exp(i,k,4,2)=xc*xc

          xc=co2exp(i,k,4,2)*co2exp(i,k,4,2)
          xc=xc*xc
          co2exp(i,k,5,2)=xc*xc

          xc=co2exp(i,k,5,2)*co2exp(i,k,5,2)
          xc=xc*xc
          co2exp(i,k,6,2)=xc*xc

       ENDDO
    ENDDO

  END SUBROUTINE co2exps

  !c**********************************************************************
  SUBROUTINE n2oexps(ib,m,np,dn2o,pa,dt,n2oexp)
    !c**********************************************************************
    !c   Compute n2o exponentials for individual layers 
    !c
    !c---- input parameters
    !c  spectral band (ib)
    !c  number of grid intervals (m)
    !c  number of layers (np)
    !c  layer n2o amount (dn2o)
    !c  layer pressure (pa)
    !c  layer temperature minus 250K (dt)
    !c
    !c---- output parameters
    !c  2 or 4 exponentials for each layer (n2oexp)
    !c**********************************************************************
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !ctar       integer ib,m,np,i,k

    INTEGER, INTENT(IN) :: ib,m,np
    INTEGER :: i,k

    !c---- input parameters -----

    !ctar       real dn2o(m,np),pa(m,np),dt(m,np)

    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np):: dn2o,pa,dt

    !c---- output parameters -----

    !ctar       real n2oexp(m,np,4)

    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np,4):: n2oexp

    !c---- temporary arrays -----

    !ctar       real xc,xc1,xc2

    REAL(KIND=r8):: xc,xc1,xc2

    !c-----Scaling and absorption data are given in Table 5.
    !c     Transmittances are computed using Eqs. (8.21) and (8.22).

    DO k=1,np
       DO i=1,m

          !c-----four exponential by powers of 21 for band 6.

          IF (ib.EQ.6) THEN

             xc=dn2o(i,k)*(1.0_r8+(1.9297e-3_r8+4.3750e-6_r8*dt(i,k))*dt(i,k))
             n2oexp(i,k,1)=EXP(-xc*6.31582e-2_r8)

             xc=n2oexp(i,k,1)*n2oexp(i,k,1)*n2oexp(i,k,1)
             xc1=xc*xc
             xc2=xc1*xc1
             n2oexp(i,k,2)=xc*xc1*xc2

             !c-----four exponential by powers of 8 for band 7

          ELSE

             xc=dn2o(i,k)*(pa(i,k)/500.0_r8)**0.48_r8  &
                  *(1.0_r8+(1.3804e-3_r8+7.4838e-6_r8*dt(i,k))*dt(i,k))
             n2oexp(i,k,1)=EXP(-xc*5.35779e-2_r8)

             xc=n2oexp(i,k,1)*n2oexp(i,k,1)
             xc=xc*xc
             n2oexp(i,k,2)=xc*xc
             xc=n2oexp(i,k,2)*n2oexp(i,k,2)
             xc=xc*xc
             n2oexp(i,k,3)=xc*xc
             xc=n2oexp(i,k,3)*n2oexp(i,k,3)
             xc=xc*xc
             n2oexp(i,k,4)=xc*xc

          ENDIF

       ENDDO
    ENDDO

  END SUBROUTINE n2oexps

  !c**********************************************************************
  SUBROUTINE ch4exps(ib,m,np,dch4,pa,dt,ch4exp)
    !c**********************************************************************
    !c   Compute ch4 exponentials for individual layers
    !c
    !c---- input parameters
    !c  spectral band (ib)
    !c  number of grid intervals (m)
    !c  number of layers (np)
    !c  layer ch4 amount (dch4)
    !c  layer pressure (pa)
    !c  layer temperature minus 250K (dt)
    !c
    !c---- output parameters
    !c  1 or 4 exponentials for each layer (ch4exp)
    !c**********************************************************************
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !ctar      integer ib,m,np,i,k

    INTEGER, INTENT(IN) :: ib,m,np
    INTEGER :: i,k

    !c---- input parameters -----

    !ctar      real dch4(m,np),pa(m,np),dt(m,np)

    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np):: dch4,pa,dt

    !c---- output parameters -----

    !ctar      real ch4exp(m,np,4)

    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np,4):: ch4exp

    !c---- temporary arrays -----

    !ctar      real xc

    REAL(KIND=r8):: xc

    !c*****  Scaling and absorption data are given in Table 5  *****

    DO k=1,np
       DO i=1,m

          !c-----four exponentials for band 6

          IF (ib.EQ.6) THEN

             xc=dch4(i,k)*(1.0_r8+(1.7007e-2_r8+1.5826e-4_r8*dt(i,k))*dt(i,k))
             ch4exp(i,k,1)=EXP(-xc*5.80708e-3_r8)

             !c-----four exponentials by powers of 12 for band 7

          ELSE

             xc=dch4(i,k)*(pa(i,k)/500.0_r8)**0.65_r8  &
                  *(1.0_r8+(5.9590e-4_r8-2.2931e-6_r8*dt(i,k))*dt(i,k))
             ch4exp(i,k,1)=EXP(-xc*6.29247e-2_r8)

             xc=ch4exp(i,k,1)*ch4exp(i,k,1)*ch4exp(i,k,1)
             xc=xc*xc
             ch4exp(i,k,2)=xc*xc

             xc=ch4exp(i,k,2)*ch4exp(i,k,2)*ch4exp(i,k,2)
             xc=xc*xc
             ch4exp(i,k,3)=xc*xc

             xc=ch4exp(i,k,3)*ch4exp(i,k,3)*ch4exp(i,k,3)
             xc=xc*xc
             ch4exp(i,k,4)=xc*xc

          ENDIF

       ENDDO
    ENDDO

  END SUBROUTINE ch4exps

  !c**********************************************************************
  SUBROUTINE comexps(ib,m,np,dcom,dt,comexp)
    !c**********************************************************************
    !c   Compute co2-minor exponentials for individual layers using 
    !c   Eqs. (8.21) and (8.22).
    !c
    !c---- input parameters
    !c  spectral band (ib)
    !c  number of grid intervals (m)
    !c  number of layers (np)
    !c  layer co2 amount (dcom)
    !c  layer temperature minus 250K (dt)
    !c
    !c---- output parameters
    !c  6 exponentials for each layer (comexp)
    !c**********************************************************************
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !ctar      integer ib,m,np,i,k,ik

    INTEGER, INTENT(IN) :: ib,m,np
    INTEGER :: i,k,ik

    !c---- input parameters -----

    !ctar      real dcom(m,np),dt(m,np)

    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np):: dcom,dt

    !c---- output parameters -----

    !ctar      real comexp(m,np,6)

    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np,6):: comexp

    !c---- temporary arrays -----

    !ctar      real xc

    REAL(KIND=r8):: xc

    !c*****  Scaling and absorpton data are given in Table 6  *****

    DO k=1,np
       DO i=1,m

          IF (ib.EQ.4) THEN
             xc=dcom(i,k)*(1.0_r8+(3.5775e-2_r8+4.0447e-4_r8*dt(i,k))*dt(i,k))
          ENDIF

          IF (ib.EQ.5) THEN
             xc=dcom(i,k)*(1.0_r8+(3.4268e-2_r8+3.7401e-4_r8*dt(i,k))*dt(i,k))
          ENDIF

          comexp(i,k,1)=EXP(-xc*1.922e-7_r8)

          DO ik=2,6
             xc=comexp(i,k,ik-1)*comexp(i,k,ik-1)
             xc=xc*xc
             comexp(i,k,ik)=xc*comexp(i,k,ik-1)
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE comexps

  !c**********************************************************************
  SUBROUTINE cfcexps(ib,m,np,a1,b1,fk1,a2,b2,fk2,dcfc,dt,cfcexp)
    !c**********************************************************************
    !c   compute cfc(-11, -12, -22) exponentials for individual layers.
    !c
    !c---- input parameters
    !c  spectral band (ib)
    !c  number of grid intervals (m)
    !c  number of layers (np)
    !c  parameters for computing the scaled cfc amounts
    !c             for temperature scaling (a1,b1,a2,b2)
    !c  the absorption coefficients for the
    !c     first k-distribution function due to cfcs (fk1,fk2)
    !c  layer cfc amounts (dcfc)
    !c  layer temperature minus 250K (dt)
    !c
    !c---- output parameters
    !c  1 exponential for each layer (cfcexp)
    !c**********************************************************************
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !ctar      integer ib,m,np,i,k

    INTEGER, INTENT(IN) :: ib,m,np
    INTEGER :: i,k

    !c---- input parameters -----

    !ctar      real dcfc(m,np),dt(m,np)

    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np):: dcfc,dt

    !c---- output parameters -----

    !ctar      real cfcexp(m,np)

    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np):: cfcexp

    !c---- static data -----

    !ctar      real a1,b1,fk1,a2,b2,fk2

    REAL(KIND=r8), INTENT(IN):: a1,b1,fk1,a2,b2,fk2

    !c---- temporary arrays -----

    !ctar      real xf

    REAL(KIND=r8):: xf

    !c**********************************************************************

    DO k=1,np
       DO i=1,m

          !c-----compute the scaled cfc amount (xf) and exponential (cfcexp)

          IF (ib.EQ.4) THEN
             xf=dcfc(i,k)*(1.0_r8+(a1+b1*dt(i,k))*dt(i,k))
             cfcexp(i,k)=EXP(-xf*fk1)
          ELSE
             xf=dcfc(i,k)*(1.0_r8+(a2+b2*dt(i,k))*dt(i,k))
             cfcexp(i,k)=EXP(-xf*fk2)
          ENDIF

       ENDDO
    ENDDO

  END SUBROUTINE cfcexps

  !c**********************************************************************
  SUBROUTINE b10exps(m,np,dh2o,dcont,dco2,dn2o,pa,dt  &
       ,h2oexp,conexp,co2exp,n2oexp)
    !c**********************************************************************
    !c   Compute band3a exponentials for individual layers
    !c
    !c---- input parameters
    !c  number of grid intervals (m)
    !c  number of layers (np)
    !c  layer h2o amount for line absorption (dh2o)
    !c  layer h2o amount for continuum absorption (dcont)
    !c  layer co2 amount (dco2)
    !c  layer n2o amount (dn2o)
    !c  layer pressure (pa)
    !c  layer temperature minus 250K (dt)
    !c
    !c---- output parameters
    !c
    !c  exponentials for each layer (h2oexp,conexp,co2exp,n2oexp)
    !c**********************************************************************
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !ctar      integer m,np,i,k

    INTEGER, INTENT(IN) :: m,np
    INTEGER :: i,k

    !c---- input parameters -----

    !ctar      real dh2o(m,np),dcont(m,np),dn2o(m,np)
    !ctar      real dco2(m,np),pa(m,np),dt(m,np)

    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np):: dh2o,dcont,dn2o,dco2,pa,dt

    !c---- output parameters -----

    !ctar      real h2oexp(m,np,6),conexp(m,np,3),co2exp(m,np,6,2)  &
    !ctar         ,n2oexp(m,np,4)

    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np,6):: h2oexp
    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np,3):: conexp
    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np,6,2):: co2exp
    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np,4):: n2oexp

    !c---- temporary arrays -----

    !ctar      real xx,xx1,xx2,xx3

    REAL(KIND=r8):: xx,xx1,xx2,xx3

    !c**********************************************************************

    DO k=1,np
       DO i=1,m

          !c-----Compute scaled h2o-line amount for Band 10 (Eq. 4.4 and Table 3).

          xx=dh2o(i,k)*(pa(i,k)/500.0_r8) &
               *(1.0_r8+(0.0149_r8+6.20e-5_r8*dt(i,k))*dt(i,k))

          !c-----six exponentials by powers of 8

          h2oexp(i,k,1)=EXP(-xx*0.10624_r8)

          xx=h2oexp(i,k,1)*h2oexp(i,k,1)
          xx=xx*xx
          h2oexp(i,k,2)=xx*xx

          xx=h2oexp(i,k,2)*h2oexp(i,k,2)
          xx=xx*xx
          h2oexp(i,k,3)=xx*xx

          xx=h2oexp(i,k,3)*h2oexp(i,k,3)
          xx=xx*xx
          h2oexp(i,k,4)=xx*xx

          xx=h2oexp(i,k,4)*h2oexp(i,k,4)
          xx=xx*xx
          h2oexp(i,k,5)=xx*xx

          !c-----one exponential of h2o continuum for sub-band 3a (Table 9).

          conexp(i,k,1)=EXP(-dcont(i,k)*109.0_r8)

          !c-----Scaled co2 amount for the Band 10 (Eq. 4.4, Tables 3 and 6).

          xx=dco2(i,k)*(pa(i,k)/300.0_r8)**0.5_r8  &
               *(1.0_r8+(0.0179_r8+1.02e-4_r8*dt(i,k))*dt(i,k))

          !c-----six exponentials by powers of 8

          co2exp(i,k,1,1)=EXP(-xx*2.656e-5_r8)

          xx=co2exp(i,k,1,1)*co2exp(i,k,1,1)
          xx=xx*xx
          co2exp(i,k,2,1)=xx*xx

          xx=co2exp(i,k,2,1)*co2exp(i,k,2,1)
          xx=xx*xx
          co2exp(i,k,3,1)=xx*xx

          xx=co2exp(i,k,3,1)*co2exp(i,k,3,1)
          xx=xx*xx
          co2exp(i,k,4,1)=xx*xx

          xx=co2exp(i,k,4,1)*co2exp(i,k,4,1)
          xx=xx*xx
          co2exp(i,k,5,1)=xx*xx

          xx=co2exp(i,k,5,1)*co2exp(i,k,5,1)
          xx=xx*xx
          co2exp(i,k,6,1)=xx*xx

          !c-----Compute the scaled n2o amount for Band 10 (Table 5).

          xx=dn2o(i,k)*(1.0_r8+(1.4476e-3_r8+3.6656e-6_r8*dt(i,k))*dt(i,k))

          !c-----Two exponentials by powers of 58

          n2oexp(i,k,1)=EXP(-xx*0.25238_r8)

          xx=n2oexp(i,k,1)*n2oexp(i,k,1)
          xx1=xx*xx
          xx1=xx1*xx1
          xx2=xx1*xx1
          xx3=xx2*xx2
          n2oexp(i,k,2)=xx*xx1*xx2*xx3

       ENDDO
    ENDDO

  END SUBROUTINE b10exps

  !c**********************************************************************
  SUBROUTINE tablup(k2,m,np,nx,nh,dw,p,dt,s1,s2,s3,w1,p1, &
       dwe,dpe,coef1,coef2,coef3,tran)
    !c**********************************************************************
    !c   Compute water vapor, co2 and o3 transmittances between level
    !c   k1 and and level k2 for m soundings, using table look-up.
    !c
    !c   Calculations follow Eq. (4.16).
    !c
    !c---- input ---------------------
    !c
    !c  index for level (k2)
    !c  number of grid intervals (m)
    !c  number of atmospheric layers (np)
    !c  number of pressure intervals in the table (nx)
    !c  number of absorber amount intervals in the table (nh)
    !c  layer absorber amount (dw)
    !c  layer pressure in mb (p)
    !c  deviation of layer temperature from 250K (dt)
    !c  first value of absorber amount (log10) in the table (w1) 
    !c  first value of pressure (log10) in the table (p1) 
    !c  size of the interval of absorber amount (log10) in the table (dwe)
    !c  size of the interval of pressure (log10) in the table (dpe)
    !c  pre-computed coefficients (coef1, coef2, and coef3)
    !c
    !c---- updated ---------------------
    !c
    !c  column integrated absorber amount (s1)
    !c  absorber-weighted column pressure (s2)
    !c  absorber-weighted column temperature (s3)
    !c  transmittance (tran)
    !c
    !c  Note: Units of s1 are g/cm**2 for water vapor and
    !c       (cm-atm)stp for co2 and o3.
    !c   
    !c**********************************************************************
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !ctar       integer k2,m,np,nx,nh,i

    INTEGER, INTENT(IN) :: k2,m,np,nx,nh
    INTEGER :: i

    !c---- input parameters -----

    !ctar       real w1,p1,dwe,dpe
    !ctar       real dw(m,np),p(m,np),dt(m,np)
    !ctar       real coef1(nx,nh),coef2(nx,nh),coef3(nx,nh)

    REAL(KIND=r8), INTENT(IN) :: w1,p1,dwe,dpe
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np) :: dw,p,dt
    REAL(KIND=r8), INTENT(IN), DIMENSION(nx,nh) :: coef1,coef2,coef3

    !c---- update parameter -----

    !ctar       real s1(m),s2(m),s3(m),tran(m,np+1)

    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m)::s1,s2,s3
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1):: tran

    !c---- temporary variables -----

    !ctar       real we,pe,fw,fp,pa,pb,pc,ax,ba,bb,t1,ca,cb,t2,x1,x2,x3
    !ctar       integer iw,ip

    REAL(KIND=r8):: we,pe,fw,fp,pa,pb,pc,ax,ba,bb,t1,ca,cb,t2,x1,x2,x3
    INTEGER :: iw,ip


    !c-----Compute effective pressure (x2) and temperature (x3) following 
    !c     Eqs. (8.28) and (8.29)

    DO i=1,m

       s1(i)=s1(i)+dw(i,k2-1)
       s2(i)=s2(i)+p(i,k2-1)*dw(i,k2-1)
       s3(i)=s3(i)+dt(i,k2-1)*dw(i,k2-1)

       x1=s1(i)
       x2=s2(i)/s1(i)
       x3=s3(i)/s1(i)

       !c-----normalize we and pe

       we=(LOG10(x1)-w1)/dwe
       pe=(LOG10(x2)-p1)/dpe

       !c-----restrict the magnitudes of the normalized we and pe.

       we=MIN(we,float(nh-1))
       pe=MIN(pe,float(nx-1))

       !c-----assign iw and ip and compute the distance of we and pe 
       !c     from iw and ip.

       iw=INT(we+1.0_r8)
       iw=MIN(iw,nh-1)
       iw=MAX(iw, 2)
       fw=we-float(iw-1)

       ip=INT(pe+1.0_r8)
       ip=MIN(ip,nx-1)
       ip=MAX(ip, 1)
       fp=pe-float(ip-1)

       !c-----linear interpolation in pressure

       pa = coef1(ip,iw-1)*(1.0_r8-fp)+coef1(ip+1,iw-1)*fp
       pb = coef1(ip,  iw)*(1.0_r8-fp)+coef1(ip+1,  iw)*fp
       pc = coef1(ip,iw+1)*(1.0_r8-fp)+coef1(ip+1,iw+1)*fp

       !c-----quadratic interpolation in absorber amount for coef1

       ax = (-pa*(1.0_r8-fw)+pc*(1.0_r8+fw)) *fw*0.5_r8 + pb*(1.0_r8-fw*fw)

       !c-----linear interpolation in absorber amount for coef2 and coef3

       ba = coef2(ip,  iw)*(1.0_r8-fp)+coef2(ip+1,  iw)*fp
       bb = coef2(ip,iw+1)*(1.0_r8-fp)+coef2(ip+1,iw+1)*fp
       t1 = ba*(1.0_r8-fw) + bb*fw

       ca = coef3(ip,  iw)*(1.0_r8-fp)+coef3(ip+1,  iw)*fp
       cb = coef3(ip,iw+1)*(1.0_r8-fp)+coef3(ip+1,iw+1)*fp
       t2 = ca*(1.0_r8-fw) + cb*fw

       !c-----update the total transmittance between levels k1 and k2

       tran(i,k2)= (ax + (t1+t2*x3) * x3)*tran(i,k2)
       tran(i,k2)=MIN(tran(i,k2),0.9999999_r8)
       tran(i,k2)=MAX(tran(i,k2),0.0000001_r8)

    ENDDO

  END SUBROUTINE tablup

  !c**********************************************************************
  SUBROUTINE h2okdis(ib,m,np,k,fkw,gkw,ne,h2oexp,conexp,  &
       th2o,tcon,tran)
    !c**********************************************************************
    !c   compute water vapor transmittance between levels k1 and k2 for
    !c   m soundings, using the k-distribution method.
    !c
    !c---- input parameters
    !c  spectral band (ib)
    !c  number of grid intervals (m)
    !c  number of levels (np)
    !c  current level (k)
    !c  planck-weighted k-distribution function due to
    !c    h2o line absorption (fkw)
    !c  planck-weighted k-distribution function due to
    !c    h2o continuum absorption (gkw)
    !c  number of terms used in each band to compute water vapor
    !c     continuum transmittance (ne)
    !c  exponentials for line absorption (h2oexp) 
    !c  exponentials for continuum absorption (conexp) 
    !c
    !c---- updated parameters
    !c  transmittance between levels k1 and k2 due to
    !c    water vapor line absorption (th2o)
    !c  transmittance between levels k1 and k2 due to
    !c    water vapor continuum absorption (tcon)
    !c  total transmittance (tran)
    !c
    !c**********************************************************************
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !c---- input parameters ------

    !ctar      integer ib,m,np,k,ne
    !ctar      real conexp(m,np,3),h2oexp(m,np,6)
    !ctar      real  fkw(6,9),gkw(6,3)

    INTEGER, INTENT(IN) :: ib,m,np,k,ne
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,3):: conexp
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,6):: h2oexp
    REAL(KIND=r8), INTENT(IN), DIMENSION(6,9):: fkw
    REAL(KIND=r8), INTENT(IN), DIMENSION(6,3):: gkw

    !c---- updated parameters -----

    !ctar      real th2o(m,6),tcon(m,3),tran(m,np+1)

    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,6):: th2o
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,3):: tcon
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1):: tran

    !c---- temporary arrays -----

    !ctar      real trnth2o
    !ctar      integer i

    REAL(KIND=r8):: trnth2o
    INTEGER :: i


    !c-----tco2 are the six exp factors between levels k1 and k2 
    !c     tran is the updated total transmittance between levels k1 and k2

    !c-----th2o is the 6 exp factors between levels k1 and k2 due to
    !c     h2o line absorption. 

    !c-----tcon is the 3 exp factors between levels k1 and k2 due to
    !c     h2o continuum absorption.

    !c-----trnth2o is the total transmittance between levels k1 and k2 due
    !c     to both line and continuum absorption.



    !c-----Compute th2o following Eq. (8.23).

    DO i=1,m
       th2o(i,1) = th2o(i,1)*h2oexp(i,k,1)
       th2o(i,2) = th2o(i,2)*h2oexp(i,k,2)
       th2o(i,3) = th2o(i,3)*h2oexp(i,k,3)
       th2o(i,4) = th2o(i,4)*h2oexp(i,k,4)
       th2o(i,5) = th2o(i,5)*h2oexp(i,k,5)
       th2o(i,6) = th2o(i,6)*h2oexp(i,k,6)
    ENDDO


    IF (ne.EQ.0) THEN

       !c-----Compute trnh2o following Eq. (8.25). fkw is given in Table 4.

       DO i=1,m

          trnth2o      =(fkw(1,ib)*th2o(i,1) &
               + fkw(2,ib)*th2o(i,2) &
               + fkw(3,ib)*th2o(i,3) &
               + fkw(4,ib)*th2o(i,4) &
               + fkw(5,ib)*th2o(i,5) &
               + fkw(6,ib)*th2o(i,6))

          tran(i,k+1)=tran(i,k+1)*trnth2o

       ENDDO

    ELSEIF (ne.EQ.1) THEN

       !c-----Compute trnh2o following Eqs. (8.25) and (4.27).

       DO i=1,m

          tcon(i,1)= tcon(i,1)*conexp(i,k,1)

          trnth2o      =(fkw(1,ib)*th2o(i,1)  &
               + fkw(2,ib)*th2o(i,2)  &
               + fkw(3,ib)*th2o(i,3)  &
               + fkw(4,ib)*th2o(i,4)  &
               + fkw(5,ib)*th2o(i,5)  &
               + fkw(6,ib)*th2o(i,6))*tcon(i,1)

          tran(i,k+1)=tran(i,k+1)*trnth2o

       ENDDO

    ELSE

       !c-----For band 3. This band is divided into 3 subbands.

       DO i=1,m

          tcon(i,1)= tcon(i,1)*conexp(i,k,1)
          tcon(i,2)= tcon(i,2)*conexp(i,k,2)
          tcon(i,3)= tcon(i,3)*conexp(i,k,3)

          !c-----Compute trnh2o following Eqs. (4.29) and (8.25).

          trnth2o      = (  gkw(1,1)*th2o(i,1) &
               + gkw(2,1)*th2o(i,2) &
               + gkw(3,1)*th2o(i,3) &
               + gkw(4,1)*th2o(i,4) &
               + gkw(5,1)*th2o(i,5) &
               + gkw(6,1)*th2o(i,6) ) * tcon(i,1) &
               + (  gkw(1,2)*th2o(i,1) &
               + gkw(2,2)*th2o(i,2) &
               + gkw(3,2)*th2o(i,3) &
               + gkw(4,2)*th2o(i,4) &
               + gkw(5,2)*th2o(i,5) &
               + gkw(6,2)*th2o(i,6) ) * tcon(i,2) &
               + (  gkw(1,3)*th2o(i,1) &
               + gkw(2,3)*th2o(i,2) &
               + gkw(3,3)*th2o(i,3) &
               + gkw(4,3)*th2o(i,4) &
               + gkw(5,3)*th2o(i,5) &
               + gkw(6,3)*th2o(i,6) ) * tcon(i,3)

          tran(i,k+1)=tran(i,k+1)*trnth2o

       ENDDO

    ENDIF

  END SUBROUTINE h2okdis

  !c**********************************************************************
  SUBROUTINE co2kdis(m,np,k,co2exp,tco2,tran)
    !c**********************************************************************
    !c   compute co2 transmittances between levels k1 and k2 for
    !c    m soundings, using the k-distribution method with linear
    !c    pressure scaling.
    !c
    !c---- input parameters
    !c   number of grid intervals (m)
    !c   number of levels (np)
    !c   current level (k)
    !c   exponentials for co2 absorption (co2exp)
    !c
    !c---- updated parameters
    !c   transmittance between levels k1 and k2 due to co2 absorption
    !c     for the various values of the absorption coefficient (tco2)
    !c   total transmittance (tran)
    !c
    !c**********************************************************************
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !ctar       integer m,np,i,k

    INTEGER, INTENT(IN) :: m,np,k
    INTEGER :: i


    !c---- input parameters -----

    !ctar       real co2exp(m,np,6,2)

    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,6,2):: co2exp

    !c---- updated parameters -----

    !ctar       real tco2(m,6,2),tran(m,np+1)

    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,6,2):: tco2
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1):: tran


    !c---- temporary arrays -----

    !ctar       real xc

    REAL(KIND=r8):: xc


    !c-----tco2 is the 6 exp factors between levels k1 and k2 computed
    !c     from Eqs. (8.23) and (8.25). Also see Eq. (4.30).
    !c     The k-distribution functions are given in Table 10.

    DO i=1,m

       !c-----band-wings

       tco2(i,1,1)=tco2(i,1,1)*co2exp(i,k,1,1)
       xc=   0.1395_r8 *tco2(i,1,1)

       tco2(i,2,1)=tco2(i,2,1)*co2exp(i,k,2,1)
       xc=xc+0.1407_r8 *tco2(i,2,1)

       tco2(i,3,1)=tco2(i,3,1)*co2exp(i,k,3,1)
       xc=xc+0.1549_r8 *tco2(i,3,1)

       tco2(i,4,1)=tco2(i,4,1)*co2exp(i,k,4,1)
       xc=xc+0.1357_r8 *tco2(i,4,1)

       tco2(i,5,1)=tco2(i,5,1)*co2exp(i,k,5,1)
       xc=xc+0.0182_r8 *tco2(i,5,1)

       tco2(i,6,1)=tco2(i,6,1)*co2exp(i,k,6,1)
       xc=xc+0.0220_r8 *tco2(i,6,1)

       !c-----band-center region

       tco2(i,1,2)=tco2(i,1,2)*co2exp(i,k,1,2)
       xc=xc+0.0766_r8 *tco2(i,1,2)

       tco2(i,2,2)=tco2(i,2,2)*co2exp(i,k,2,2)
       xc=xc+0.1372_r8 *tco2(i,2,2)

       tco2(i,3,2)=tco2(i,3,2)*co2exp(i,k,3,2)
       xc=xc+0.1189_r8 *tco2(i,3,2)

       tco2(i,4,2)=tco2(i,4,2)*co2exp(i,k,4,2)
       xc=xc+0.0335_r8 *tco2(i,4,2)

       tco2(i,5,2)=tco2(i,5,2)*co2exp(i,k,5,2)
       xc=xc+0.0169_r8 *tco2(i,5,2)

       tco2(i,6,2)=tco2(i,6,2)*co2exp(i,k,6,2)
       xc=xc+0.0059_r8 *tco2(i,6,2)

       tran(i,k+1)=tran(i,k+1)*xc

    ENDDO

  END SUBROUTINE co2kdis

  !c**********************************************************************
  SUBROUTINE n2okdis(ib,m,np,k,n2oexp,tn2o,tran)
    !c**********************************************************************
    !c   compute n2o transmittances between levels k1 and k2 for
    !c    m soundings, using the k-distribution method with linear
    !c    pressure scaling.
    !c
    !c---- input parameters
    !c   spectral band (ib)
    !c   number of grid intervals (m)
    !c   number of levels (np)
    !c   current level (k)
    !c   exponentials for n2o absorption (n2oexp)
    !c
    !c---- updated parameters
    !c   transmittance between levels k1 and k2 due to n2o absorption
    !c     for the various values of the absorption coefficient (tn2o)
    !c   total transmittance (tran)
    !c
    !c**********************************************************************
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !ctar       integer ib,m,np,i,k

    INTEGER, INTENT(IN) :: ib,m,np,k
    INTEGER :: i

    !c---- input parameters -----

    !ctar       real n2oexp(m,np,4)

    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,4):: n2oexp

    !c---- updated parameters -----

    !ctar       real tn2o(m,4),tran(m,np+1)

    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,4):: tn2o
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1):: tran

    !c---- temporary arrays -----

    !ctar       real xc

    REAL(KIND=r8):: xc


    !c-----tn2o is computed from Eq. (8.23). 
    !c     xc is the total n2o transmittance computed from (8.25)
    !c     The k-distribution functions are given in Table 5.

    DO i=1,m

       !c-----band 6

       IF (ib.EQ.6) THEN

          tn2o(i,1)=tn2o(i,1)*n2oexp(i,k,1)
          xc=   0.940414_r8*tn2o(i,1)

          tn2o(i,2)=tn2o(i,2)*n2oexp(i,k,2)
          xc=xc+0.059586_r8*tn2o(i,2)

          !c-----band 7

       ELSE

          tn2o(i,1)=tn2o(i,1)*n2oexp(i,k,1)
          xc=   0.561961_r8*tn2o(i,1)

          tn2o(i,2)=tn2o(i,2)*n2oexp(i,k,2)
          xc=xc+0.138707_r8*tn2o(i,2)

          tn2o(i,3)=tn2o(i,3)*n2oexp(i,k,3)
          xc=xc+0.240670_r8*tn2o(i,3)

          tn2o(i,4)=tn2o(i,4)*n2oexp(i,k,4)
          xc=xc+0.058662_r8*tn2o(i,4)

       ENDIF

       tran(i,k+1)=tran(i,k+1)*xc

    ENDDO

  END SUBROUTINE n2okdis

  !c**********************************************************************
  SUBROUTINE ch4kdis(ib,m,np,k,ch4exp,tch4,tran)
    !c**********************************************************************
    !c   compute ch4 transmittances between levels k1 and k2 for
    !c    m soundings, using the k-distribution method with
    !c    linear pressure scaling.
    !c
    !c---- input parameters
    !c   spectral band (ib)
    !c   number of grid intervals (m)
    !c   number of levels (np)
    !c   current level (k)
    !c   exponentials for ch4 absorption (ch4exp)
    !c
    !c---- updated parameters
    !c   transmittance between levels k1 and k2 due to ch4 absorption
    !c     for the various values of the absorption coefficient (tch4)
    !c   total transmittance (tran)
    !c
    !c**********************************************************************
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !ctar      integer ib,m,np,i,k

    INTEGER, INTENT(IN) :: ib,m,np,k
    INTEGER :: i

    !c---- input parameters -----

    !ctar      real ch4exp(m,np,4)

    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,4):: ch4exp

    !c---- updated parameters -----

    !ctar      real tch4(m,4),tran(m,np+1)

    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,4):: tch4
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1):: tran

    !c---- temporary arrays -----

    !ctar      real xc

    REAL(KIND=r8):: xc

    !c-----tch4 is computed from Eq. (8.23). 
    !c     xc is the total ch4 transmittance computed from (8.25)
    !c     The k-distribution functions are given in Table 5.

    DO i=1,m

       !c-----band 6

       IF (ib.EQ.6) THEN

          tch4(i,1)=tch4(i,1)*ch4exp(i,k,1)
          xc= tch4(i,1)

          !c-----band 7

       ELSE

          tch4(i,1)=tch4(i,1)*ch4exp(i,k,1)
          xc=   0.610650_r8*tch4(i,1)

          tch4(i,2)=tch4(i,2)*ch4exp(i,k,2)
          xc=xc+0.280212_r8*tch4(i,2)

          tch4(i,3)=tch4(i,3)*ch4exp(i,k,3)
          xc=xc+0.107349_r8*tch4(i,3)

          tch4(i,4)=tch4(i,4)*ch4exp(i,k,4)
          xc=xc+0.001789_r8*tch4(i,4)

       ENDIF

       tran(i,k+1)=tran(i,k+1)*xc

    ENDDO

  END SUBROUTINE ch4kdis

  !c**********************************************************************
  SUBROUTINE comkdis(ib,m,np,k,comexp,tcom,tran)
    !c**********************************************************************
    !c  compute co2-minor transmittances between levels k1 and k2
    !c   for m soundings, using the k-distribution method
    !c   with linear pressure scaling.
    !c
    !c---- input parameters
    !c   spectral band (ib)
    !c   number of grid intervals (m)
    !c   number of levels (np)
    !c   current level (k)
    !c   exponentials for co2-minor absorption (comexp)
    !c
    !c---- updated parameters
    !c   transmittance between levels k1 and k2 due to co2-minor absorption
    !c     for the various values of the absorption coefficient (tcom)
    !c   total transmittance (tran)
    !c
    !c**********************************************************************
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !ctar      integer ib,m,np,i,k

    INTEGER, INTENT(IN) :: ib,m,np,k
    INTEGER :: i

    !c---- input parameters -----

    !ctar      real comexp(m,np,6)

    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,6):: comexp

    !c---- updated parameters -----

    !ctar      real tcom(m,6),tran(m,np+1)

    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,6):: tcom
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1):: tran

    !c---- temporary arrays -----

    !ctar      real xc

    REAL(KIND=r8):: xc

    !c-----tcom is computed from Eq. (8.23). 
    !c     xc is the total co2 transmittance computed from (8.25)
    !c     The k-distribution functions are given in Table 6.

    DO i=1,m

       !c-----band 4

       IF (ib.EQ.4) THEN

          tcom(i,1)=tcom(i,1)*comexp(i,k,1)
          xc=   0.12159_r8*tcom(i,1)
          tcom(i,2)=tcom(i,2)*comexp(i,k,2)
          xc=xc+0.24359_r8*tcom(i,2)
          tcom(i,3)=tcom(i,3)*comexp(i,k,3)
          xc=xc+0.24981_r8*tcom(i,3)
          tcom(i,4)=tcom(i,4)*comexp(i,k,4)
          xc=xc+0.26427_r8*tcom(i,4)
          tcom(i,5)=tcom(i,5)*comexp(i,k,5)
          xc=xc+0.07807_r8*tcom(i,5)
          tcom(i,6)=tcom(i,6)*comexp(i,k,6)
          xc=xc+0.04267_r8*tcom(i,6)

          !c-----band 5

       ELSE

          tcom(i,1)=tcom(i,1)*comexp(i,k,1)
          xc=   0.06869_r8*tcom(i,1)
          tcom(i,2)=tcom(i,2)*comexp(i,k,2)
          xc=xc+0.14795_r8*tcom(i,2)
          tcom(i,3)=tcom(i,3)*comexp(i,k,3)
          xc=xc+   0.19512_r8*tcom(i,3)
          tcom(i,4)=tcom(i,4)*comexp(i,k,4)
          xc=xc+   0.33446_r8*tcom(i,4)
          tcom(i,5)=tcom(i,5)*comexp(i,k,5)
          xc=xc+   0.17199_r8*tcom(i,5)
          tcom(i,6)=tcom(i,6)*comexp(i,k,6)
          xc=xc+   0.08179_r8*tcom(i,6)
       ENDIF

       tran(i,k+1)=tran(i,k+1)*xc

    ENDDO

  END SUBROUTINE comkdis

  !c**********************************************************************
  SUBROUTINE cfckdis(m,np,k,cfcexp,tcfc,tran)
    !c**********************************************************************
    !c  compute cfc-(11,12,22) transmittances between levels k1 and k2
    !c   for m soundings, using the k-distribution method with
    !c   linear pressure scaling.
    !c
    !c---- input parameters
    !c   number of grid intervals (m)
    !c   number of levels (np)
    !c   current level (k)
    !c   exponentials for cfc absorption (cfcexp)
    !c
    !c---- updated parameters
    !c   transmittance between levels k1 and k2 due to cfc absorption
    !c     for the various values of the absorption coefficient (tcfc)
    !c   total transmittance (tran)
    !c
    !c**********************************************************************
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !ctar      integer m,np,i,k

    INTEGER, INTENT(IN) :: m,np,k
    INTEGER :: i

    !c---- input parameters -----

    !ctar      real cfcexp(m,np)

    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np):: cfcexp

    !c---- updated parameters -----

    !ctar      real tcfc(m),tran(m,np+1)

    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m):: tcfc
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1):: tran

    !c-----tcfc is the exp factors between levels k1 and k2. 

    DO i=1,m

       tcfc(i)=tcfc(i)*cfcexp(i,k)
       tran(i,k+1)=tran(i,k+1)*tcfc(i)

    ENDDO

  END SUBROUTINE cfckdis

  !c**********************************************************************
  SUBROUTINE b10kdis(m,np,k,h2oexp,conexp,co2exp,n2oexp &
       ,th2o,tcon,tco2,tn2o,tran)
    !c**********************************************************************
    !c
    !c   compute h2o (line and continuum),co2,n2o transmittances between
    !c   levels k1 and k2 for m soundings, using the k-distribution
    !c   method with linear pressure scaling.
    !c
    !c---- input parameters
    !c   number of grid intervals (m)
    !c   number of levels (np)
    !c   current level (k)
    !c   exponentials for h2o line absorption (h2oexp)
    !c   exponentials for h2o continuum absorption (conexp)
    !c   exponentials for co2 absorption (co2exp)
    !c   exponentials for n2o absorption (n2oexp)
    !c
    !c---- updated parameters
    !c   transmittance between levels k1 and k2 due to h2o line absorption
    !c     for the various values of the absorption coefficient (th2o)
    !c   transmittance between levels k1 and k2 due to h2o continuum
    !c     absorption for the various values of the absorption
    !c     coefficient (tcon)
    !c   transmittance between levels k1 and k2 due to co2 absorption
    !c     for the various values of the absorption coefficient (tco2)
    !c   transmittance between levels k1 and k2 due to n2o absorption
    !c     for the various values of the absorption coefficient (tn2o)
    !c   total transmittance (tran)
    !c
    !c**********************************************************************
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !ctar      integer m,np,i,k

    INTEGER, INTENT(IN) :: m,np,k
    INTEGER :: i


    !c---- input parameters -----

    !ctar      real h2oexp(m,np,6),conexp(m,np,3),co2exp(m,np,6,2) &
    !ctar         ,n2oexp(m,np,4)

    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,6):: h2oexp
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,3):: conexp
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,6,2):: co2exp
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,4):: n2oexp

    !c---- updated parameters -----

    !ctar      real th2o(m,6),tcon(m,3),tco2(m,6,2),tn2o(m,4)  &
    !ctar         ,tran(m,np+1)

    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,6):: th2o
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,3):: tcon
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,6,2):: tco2
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,4):: tn2o
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1):: tran

    !c---- temporary arrays -----

    !ctar      real xx

    REAL(KIND=r8):: xx

    !c-----For h2o line. The k-distribution functions are given in Table 4.

    DO i=1,m

       th2o(i,1)=th2o(i,1)*h2oexp(i,k,1)
       xx=   0.3153_r8*th2o(i,1)

       th2o(i,2)=th2o(i,2)*h2oexp(i,k,2)
       xx=xx+0.4604_r8*th2o(i,2)

       th2o(i,3)=th2o(i,3)*h2oexp(i,k,3)
       xx=xx+0.1326_r8*th2o(i,3)

       th2o(i,4)=th2o(i,4)*h2oexp(i,k,4)
       xx=xx+0.0798_r8*th2o(i,4)

       th2o(i,5)=th2o(i,5)*h2oexp(i,k,5)
       xx=xx+0.0119_r8*th2o(i,5)

       tran(i,k+1)=xx

    ENDDO

    !c-----For h2o continuum. Note that conexp(i,k,3) is for subband 3a.

    DO i=1,m

       tcon(i,1)=tcon(i,1)*conexp(i,k,1)
       tran(i,k+1)=tran(i,k+1)*tcon(i,1)

    ENDDO

    !c-----For co2 (Table 6)

    DO i=1,m

       tco2(i,1,1)=tco2(i,1,1)*co2exp(i,k,1,1)
       xx=    0.2673_r8*tco2(i,1,1)

       tco2(i,2,1)=tco2(i,2,1)*co2exp(i,k,2,1)
       xx=xx+ 0.2201_r8*tco2(i,2,1)

       tco2(i,3,1)=tco2(i,3,1)*co2exp(i,k,3,1)
       xx=xx+ 0.2106_r8*tco2(i,3,1)

       tco2(i,4,1)=tco2(i,4,1)*co2exp(i,k,4,1)
       xx=xx+ 0.2409_r8*tco2(i,4,1)

       tco2(i,5,1)=tco2(i,5,1)*co2exp(i,k,5,1)
       xx=xx+ 0.0196_r8*tco2(i,5,1)

       tco2(i,6,1)=tco2(i,6,1)*co2exp(i,k,6,1)
       xx=xx+ 0.0415_r8*tco2(i,6,1)

       tran(i,k+1)=tran(i,k+1)*xx

    ENDDO

    !c-----For n2o (Table 5)

    DO i=1,m

       tn2o(i,1)=tn2o(i,1)*n2oexp(i,k,1)
       xx=   0.970831_r8*tn2o(i,1)

       tn2o(i,2)=tn2o(i,2)*n2oexp(i,k,2)
       xx=xx+0.029169_r8*tn2o(i,2)
       tran(i,k+1)=tran(i,k+1)*(xx-1.0_r8)

    ENDDO

  END SUBROUTINE b10kdis

  !c***********************************************************************
  SUBROUTINE cldovlp (m,np,k2,ict,icb,it,im,ib,itx,imx,ibx,  &
       cldhi,cldmd,cldlw,fcld,tcldlyr,fclr)

    !c***********************************************************************
    !c     compute the fractional clear line-of-sight between levels k1
    !c     and k2 following Eqs.(6.18)-(6.21).
    !c
    !c input parameters
    !c
    !c  m:       number of soundings
    !c  np:      number of layers
    !c  k2:      index for the level
    !c  ict:     the level separating high and middle clouds
    !c  icb:     the level separating middle and low clouds
    !c  it:      number of cloudy layers in the high-cloud group
    !c  im:      number of cloudy layers in the middle-cloud group
    !c  ib:      number of cloudy layers in the low-cloud group
    !c  fcld:    fractional cloud cover of a layer
    !c  tcldlyr: transmittance of a cloud layer
    !c  
    !c output parameter
    !c
    !c  fclr:    clear line-of-sight between levels k1 and k2
    !c***********************************************************************

    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !ctar      integer m,np,k2,ict,icb

    INTEGER, INTENT(IN) :: m,np,k2
    INTEGER, INTENT(IN) :: ict(m)
    INTEGER, INTENT(IN) :: icb(m)

    !ctar      integer i,j,k,ii,it(m),im(m),ib(m),itx(m,np),imx(m,np),ibx(m,np)

    INTEGER :: i,j,k,ii
    INTEGER, INTENT(INOUT), DIMENSION(m):: it,im,ib
    INTEGER, INTENT(INOUT), DIMENSION(m,np)::  itx,imx,ibx

    !ctar      real cldhi(m),cldmd(m),cldlw(m)

    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m):: cldhi,cldmd,cldlw

    !ctar      real fcld(m,np),tcldlyr(m,np),fclr(m,np+1)

    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np):: fcld,tcldlyr

    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np+1):: fclr 

    !c***********************************************************************
    DO i=1,m

       !c-----For high clouds
       !c     "it" is the number of high-cloud layers

       IF (k2.LE.ict(i)) THEN
          IF(fcld(i,k2-1).GT.0.001) THEN

             it(i)=it(i)+1
             ii=it(i)
             itx(i,ii)=k2-1

             IF (ii .EQ. 1) go to 11

             !c-----Rearrange the order of cloud layers with increasing cloud amount

             DO k=1,ii-1
                j=itx(i,k)
                IF(fcld(i,j).GT.fcld(i,k2-1)) THEN
                   DO j=ii-1,k,-1
                      itx(i,j+1)=itx(i,j)
                   ENDDO
                   itx(i,k)=k2-1
                   go to 11
                ENDIF
             ENDDO

11           CONTINUE

             !c-----compute equivalent black-body high cloud amount

             cldhi(i)=0.0_r8
             DO k=1,ii
                j=itx(i,k)
                cldhi(i)=fcld(i,j)-tcldlyr(i,j)*(fcld(i,j)-cldhi(i))
             ENDDO

          ENDIF
       ENDIF

       !c-----For middle clouds
       !c     "im" is the number of middle-cloud layers

       IF (k2.GT.ict(i) .AND. k2.LE.icb(i)) THEN
          IF(fcld(i,k2-1).GT.0.001_r8) THEN

             im(i)=im(i)+1
             ii=im(i)
             imx(i,ii)=k2-1

             IF (ii .EQ. 1) go to 21

             !c-----Rearrange the order of cloud layers with increasing cloud amount

             DO k=1,ii-1
                j=imx(i,k)
                IF(fcld(i,j).GT.fcld(i,k2-1)) THEN
                   DO j=ii-1,k,-1
                      imx(i,j+1)=imx(i,j)
                   ENDDO
                   imx(i,k)=k2-1
                   go to 21
                ENDIF
             ENDDO

21           CONTINUE

             !c-----compute equivalent black-body middle cloud amount

             cldmd(i)=0.0_r8
             DO k=1,ii
                j=imx(i,k)
                cldmd(i)=fcld(i,j)-tcldlyr(i,j)*(fcld(i,j)-cldmd(i))
             ENDDO

          ENDIF
       ENDIF

       !c-----For low clouds
       !c     "ib" is the number of low-cloud layers

       IF (k2.GT.icb(i)) THEN
          IF(fcld(i,k2-1).GT.0.001_r8) THEN

             ib(i)=ib(i)+1
             ii=ib(i)
             ibx(i,ii)=k2-1

             IF (ii .EQ. 1) go to 31

             !c-----Rearrange the order of cloud layers with increasing cloud amount

             DO k=1,ii-1
                j=ibx(i,k)
                IF(fcld(i,j).GT.fcld(i,k2-1)) THEN
                   DO j=ii-1,k,-1
                      ibx(i,j+1)=ibx(i,j)
                   ENDDO
                   ibx(i,k)=k2-1
                   go to 31
                ENDIF
             ENDDO

31           CONTINUE

             !c-----compute equivalent black-body low cloud amount

             cldlw(i)=0.0_r8
             DO k=1,ii
                j=ibx(i,k)
                cldlw(i)=fcld(i,j)-tcldlyr(i,j)*(fcld(i,j)-cldlw(i))
             ENDDO

          ENDIF
       ENDIF

       !c-----fclr is the equivalent clear fraction between levels k1 and k2
       !c     assuming the three cloud groups are randomly overlapped.
       !c     It follows Eqs. (6.20) and (6.21).

       fclr(i,k2)=(1.0_r8-cldhi(i))*(1.0_r8-cldmd(i))*(1.0_r8-cldlw(i))   

    ENDDO

  END SUBROUTINE cldovlp

  !c***********************************************************************
  SUBROUTINE sfcflux (ibn,m,ns,fs,tg,eg,tv,ev,rv,vege,  &
       bs,dbs,rflxs)
    !c***********************************************************************
    !c Compute emission and reflection by an homogeneous/inhomogeneous 
    !c  surface with vegetation cover.
    !c
    !c-----Input parameters
    !c  index for the spectral band (ibn)
    !c  number of grid box (m)
    !c  number of sub-grid box (ns)
    !c  fractional cover of sub-grid box (fs)
    !c  sub-grid ground temperature (tg)
    !c  sub-grid ground emissivity (eg)
    !c  sub-grid vegetation temperature (tv)
    !c  sub-grid vegetation emissivity (ev)
    !c  sub-grid vegetation reflectivity (rv)
    !c  if there is vegetation cover, vege=.true.
    !c
    !c-----Output parameters
    !c  Emission by the surface (ground+vegetation) (bs)
    !c  Derivative of bs rwt temperature (dbs)
    !c  Reflection by the surface (rflxs)

    !c**********************************************************************
    IMPLICIT NONE

    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar

    !c---- input parameters -----

    !ctar       integer ibn,m,ns

    INTEGER, INTENT(IN) :: ibn,m,ns

    !ctar       real fs(m,ns),tg(m,ns),eg(m,ns,10)

    REAL(KIND=r8), INTENT(IN), DIMENSION(m,ns):: fs,tg
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,ns,10):: eg

    !ctar       real tv(m,ns),ev(m,ns,10),rv(m,ns,10)

    REAL(KIND=r8), INTENT(IN), DIMENSION(m,ns):: tv
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,ns,10):: ev,rv

    !ctar       logical vege

    LOGICAL :: vege 

    !c---- output parameters -----

    !ctar       real bs(m),dbs(m),rflxs(m)

    REAL(KIND=r8), INTENT(OUT), DIMENSION(m):: bs,dbs,rflxs

    !c---- temporary arrays -----

    !ctar       integer i,j

    INTEGER :: i,j

    !ctar       real bg(m),dbg(m),bv(m),dbv(m),tx(m),ty(m),xx,yy,zz

    REAL(KIND=r8), DIMENSION(m):: bg,dbg,bv,dbv,tx,ty
    REAL(KIND=r8):: xx,yy,zz 

    !c*********************************************************

    IF (ns.EQ.1) THEN

       IF (.NOT.vege) THEN

          !c-----for homogeneous surface without vegetation
          !c     following Eqs. (9.4), (9.5), and (3.13)

          DO i=1,m
             tx(i)=tg(i,1)
          ENDDO

          CALL planck(ibn,m,tx,bg)
          CALL plancd(ibn,m,tx,dbg)

          DO i=1,m
             bs(i) =eg(i,1,ibn)*bg(i)
             dbs(i)=eg(i,1,ibn)*dbg(i)
             rflxs(i)=1.0_r8-eg(i,1,ibn)
          ENDDO

       ELSE

          !c-----With vegetation, following Eqs. (9.1), (9.3), and (9.13)

          DO i=1,m
             tx(i)=tg(i,1)
             ty(i)=tv(i,1)
          ENDDO

          CALL planck(ibn,m,tx,bg)
          CALL planck(ibn,m,ty,bv)
          CALL plancd(ibn,m,tx,dbg)
          CALL plancd(ibn,m,ty,dbv)

          DO i=1,m
             xx=ev(i,1,ibn)*bv(i)
             yy=1.0_r8-ev(i,1,ibn)-rv(i,1,ibn)
             zz=1.0_r8-eg(i,1,ibn)
             bs(i)=yy*(eg(i,1,ibn)*bg(i)+zz*xx)+xx

             xx=ev(i,1,ibn)*dbv(i)
             dbs(i)=yy*(eg(i,1,ibn)*dbg(i)+zz*xx)+xx

             rflxs(i)=rv(i,1,ibn)+zz*yy*yy/(1.0_r8-rv(i,1,ibn)*zz)
          ENDDO

       ENDIF

    ELSE

       !c-----for nonhomogeneous surface

       DO i=1,m
          bs(i)=0.0_r8
          dbs(i)=0.0_r8
          rflxs(i)=0.0_r8
       ENDDO

       IF(.NOT.vege) THEN

          !c-----No vegetation, following Eqs. (9.9), (9.10), and (9.13)

          DO j=1,ns

             DO i=1,m
                tx(i)=tg(i,j)
             ENDDO

             CALL planck(ibn,m,tx,bg)
             CALL plancd(ibn,m,tx,dbg)

             DO i=1,m
                bs(i)=bs(i)+fs(i,j)*eg(i,j,ibn)*bg(i)
                dbs(i)=dbs(i)+fs(i,j)*eg(i,j,ibn)*dbg(i)
                rflxs(i)=rflxs(i)+fs(i,j)*(1.0_r8-eg(i,j,ibn))
             ENDDO

          ENDDO

       ELSE

          !c-----With vegetation, following Eqs. (9.6), (9.7), and (9.13)

          DO j=1,ns
             DO i=1,m
                tx(i)=tg(i,j)
                ty(i)=tv(i,j)
             ENDDO

             CALL planck(ibn,m,tx,bg)
             CALL planck(ibn,m,ty,bv)
             CALL plancd(ibn,m,tx,dbg)
             CALL plancd(ibn,m,ty,dbv)

             DO i=1,m
                xx=ev(i,j,ibn)*bv(i)
                yy=1.0_r8-ev(i,j,ibn)-rv(i,j,ibn)
                zz=1.0_r8-eg(i,j,ibn)
                bs(i)=bs(i)+fs(i,j)*(yy*(eg(i,j,ibn)*bg(i)+zz*xx)+xx)

                xx=ev(i,j,ibn)*dbv(i)
                dbs(i)=dbs(i)+fs(i,j)*(yy*(eg(i,j,ibn)*dbg(i)+zz*xx)+xx)

                rflxs(i)=rflxs(i)+fs(i,j)*(rv(i,j,ibn)+zz*yy*yy  &
                     /(1.0_r8-rv(i,j,ibn)*zz))
             ENDDO
          ENDDO

       ENDIF

    ENDIF

  END SUBROUTINE sfcflux
  
  
  
  
     SUBROUTINE read_table()
      INTEGER  :: it ,k,j,i

      REAL(KIND=r8) :: data1(nx*3*nc)
      REAL(KIND=r8) :: data2(nx*3*no)
      REAL(KIND=r8) :: data3(nx*3*nh)
      REAL(KIND=r8) :: data4(nx*3*nh)
      REAL(KIND=r8) :: data5(nx*3*nh)

      data1(1:nx*3*nc)=RESHAPE(SOURCE=(/&
      ! data ((c1(ip,iw),iw=1,30), ip= 1, 1)/  &
        0.99985647_r8,  0.99976432_r8,  0.99963892_r8,  0.99948031_r8,  0.99927652_r8,  &
        0.99899602_r8,  0.99860001_r8,  0.99804801_r8,  0.99732202_r8,  0.99640399_r8,  &
        0.99526399_r8,  0.99384302_r8,  0.99204999_r8,  0.98979002_r8,  0.98694998_r8,  &
        0.98334998_r8,  0.97878999_r8,  0.97307003_r8,  0.96592999_r8,  0.95722002_r8,  &
        0.94660002_r8,  0.93366003_r8,  0.91777998_r8,  0.89819998_r8,  0.87419999_r8,  &
        0.84500003_r8,  0.81029999_r8,  0.76989996_r8,  0.72440004_r8,  0.67490000_r8,  &
      ! data ((c2(ip,iw),iw=1,30), ip= 1, 1)/  &
       -0.1841E-06_r8, -0.4666E-06_r8, -0.1050E-05_r8, -0.2069E-05_r8, -0.3601E-05_r8,  &
       -0.5805E-05_r8, -0.8863E-05_r8, -0.1291E-04_r8, -0.1806E-04_r8, -0.2460E-04_r8,  &
       -0.3317E-04_r8, -0.4452E-04_r8, -0.5944E-04_r8, -0.7884E-04_r8, -0.1036E-03_r8,  &
       -0.1346E-03_r8, -0.1727E-03_r8, -0.2186E-03_r8, -0.2728E-03_r8, -0.3364E-03_r8,  &
       -0.4102E-03_r8, -0.4948E-03_r8, -0.5890E-03_r8, -0.6900E-03_r8, -0.7930E-03_r8,  &
       -0.8921E-03_r8, -0.9823E-03_r8, -0.1063E-02_r8, -0.1138E-02_r8, -0.1214E-02_r8,  &
      ! data ((c3(ip,iw),iw=1,30), ip= 1, 1)/  &
        0.5821E-10_r8,  0.5821E-10_r8, -0.3201E-09_r8, -0.1804E-08_r8, -0.4336E-08_r8,  &
       -0.7829E-08_r8, -0.1278E-07_r8, -0.1847E-07_r8, -0.2827E-07_r8, -0.4495E-07_r8,  &
       -0.7126E-07_r8, -0.1071E-06_r8, -0.1524E-06_r8, -0.2160E-06_r8, -0.3014E-06_r8,  &
       -0.4097E-06_r8, -0.5349E-06_r8, -0.6718E-06_r8, -0.8125E-06_r8, -0.9755E-06_r8,  &
       -0.1157E-05_r8, -0.1339E-05_r8, -0.1492E-05_r8, -0.1563E-05_r8, -0.1485E-05_r8,  &
       -0.1210E-05_r8, -0.7280E-06_r8, -0.1107E-06_r8,  0.5369E-06_r8,  0.1154E-05_r8,  &
      ! data ((c1(ip,iw),iw=1,30), ip= 2, 2)/  &
        0.99985647_r8,  0.99976432_r8,  0.99963868_r8,  0.99947977_r8,  0.99927580_r8,  &
        0.99899501_r8,  0.99859601_r8,  0.99804401_r8,  0.99731201_r8,  0.99638498_r8,  &
        0.99523097_r8,  0.99378198_r8,  0.99194402_r8,  0.98961002_r8,  0.98664999_r8,  &
        0.98286998_r8,  0.97807002_r8,  0.97200000_r8,  0.96439999_r8,  0.95503998_r8,  &
        0.94352001_r8,  0.92931998_r8,  0.91175002_r8,  0.88989997_r8,  0.86300004_r8,  &
        0.83039999_r8,  0.79159999_r8,  0.74710000_r8,  0.69790000_r8,  0.64579999_r8,  &
      ! data ((c2(ip,iw),iw=1,30), ip= 2, 2)/  &
       -0.1831E-06_r8, -0.4642E-06_r8, -0.1048E-05_r8, -0.2067E-05_r8, -0.3596E-05_r8,  &
       -0.5797E-05_r8, -0.8851E-05_r8, -0.1289E-04_r8, -0.1802E-04_r8, -0.2454E-04_r8,  &
       -0.3307E-04_r8, -0.4435E-04_r8, -0.5916E-04_r8, -0.7842E-04_r8, -0.1031E-03_r8,  &
       -0.1342E-03_r8, -0.1725E-03_r8, -0.2189E-03_r8, -0.2739E-03_r8, -0.3386E-03_r8,  &
       -0.4138E-03_r8, -0.5003E-03_r8, -0.5968E-03_r8, -0.7007E-03_r8, -0.8076E-03_r8,  &
       -0.9113E-03_r8, -0.1007E-02_r8, -0.1096E-02_r8, -0.1181E-02_r8, -0.1271E-02_r8,  &
      ! data ((c3(ip,iw),iw=1,30)_r8, ip= 2_r8, 2)/  &
        0.5821E-10_r8,  0.5821E-10_r8, -0.3347E-09_r8, -0.1746E-08_r8, -0.4366E-08_r8,  &
       -0.7858E-08_r8, -0.1262E-07_r8, -0.1866E-07_r8, -0.2849E-07_r8, -0.4524E-07_r8,  &
       -0.7176E-07_r8, -0.1077E-06_r8, -0.1531E-06_r8, -0.2166E-06_r8, -0.3018E-06_r8,  &
       -0.4090E-06_r8, -0.5327E-06_r8, -0.6670E-06_r8, -0.8088E-06_r8, -0.9714E-06_r8,  &
       -0.1151E-05_r8, -0.1333E-05_r8, -0.1483E-05_r8, -0.1548E-05_r8, -0.1467E-05_r8,  &
       -0.1192E-05_r8, -0.7159E-06_r8, -0.1032E-06_r8,  0.5571E-06_r8,  0.1217E-05_r8,  &
      ! data ((c1(ip,iw),iw=1,30)_r8, ip= 3_r8, 3)/  &
        0.99985671_r8,  0.99976432_r8,  0.99963838_r8,  0.99947912_r8,  0.99927449_r8,  &
        0.99899203_r8,  0.99859202_r8,  0.99803501_r8,  0.99729699_r8,  0.99635702_r8,  &
        0.99518001_r8,  0.99369103_r8,  0.99178600_r8,  0.98935002_r8,  0.98623002_r8,  &
        0.98223001_r8,  0.97711003_r8,  0.97060001_r8,  0.96243000_r8,  0.95222998_r8,  &
        0.93957001_r8,  0.92379999_r8,  0.90411001_r8,  0.87959999_r8,  0.84930003_r8,  &
        0.81270003_r8,  0.76980001_r8,  0.72140002_r8,  0.66909999_r8,  0.61539996_r8,  &
      ! data ((c2(ip,iw),iw=1,30)_r8, ip= 3_r8, 3)/  &
       -0.1831E-06_r8, -0.4623E-06_r8, -0.1048E-05_r8, -0.2065E-05_r8, -0.3589E-05_r8,  &
       -0.5789E-05_r8, -0.8833E-05_r8, -0.1286E-04_r8, -0.1797E-04_r8, -0.2446E-04_r8,  &
       -0.3292E-04_r8, -0.4412E-04_r8, -0.5880E-04_r8, -0.7795E-04_r8, -0.1027E-03_r8,  &
       -0.1340E-03_r8, -0.1728E-03_r8, -0.2199E-03_r8, -0.2759E-03_r8, -0.3419E-03_r8,  &
       -0.4194E-03_r8, -0.5081E-03_r8, -0.6078E-03_r8, -0.7156E-03_r8, -0.8270E-03_r8,  &
       -0.9365E-03_r8, -0.1040E-02_r8, -0.1137E-02_r8, -0.1235E-02_r8, -0.1339E-02_r8,  &
      ! data ((c3(ip,iw),iw=1,30)_r8, ip= 3_r8, 3)/  &
        0.2910E-10_r8,  0.5821E-10_r8, -0.3201E-09_r8, -0.1732E-08_r8, -0.4307E-08_r8,  &
       -0.7843E-08_r8, -0.1270E-07_r8, -0.1882E-07_r8, -0.2862E-07_r8, -0.4571E-07_r8,  &
       -0.7225E-07_r8, -0.1082E-06_r8, -0.1535E-06_r8, -0.2171E-06_r8, -0.3021E-06_r8,  &
       -0.4084E-06_r8, -0.5302E-06_r8, -0.6615E-06_r8, -0.8059E-06_r8, -0.9668E-06_r8,  &
       -0.1146E-05_r8, -0.1325E-05_r8, -0.1468E-05_r8, -0.1530E-05_r8, -0.1448E-05_r8,  &
       -0.1168E-05_r8, -0.6907E-06_r8, -0.7148E-07_r8,  0.6242E-06_r8,  0.1357E-05_r8,  &
      ! data ((c1(ip,iw),iw=1,30)_r8, ip= 4_r8, 4)/  &
        0.99985629_r8,  0.99976349_r8,  0.99963838_r8,  0.99947798_r8,  0.99927282_r8,  &
        0.99898797_r8,  0.99858499_r8,  0.99802202_r8,  0.99727303_r8,  0.99631298_r8,  &
        0.99510002_r8,  0.99355298_r8,  0.99155599_r8,  0.98898000_r8,  0.98566002_r8,  &
        0.98136997_r8,  0.97584999_r8,  0.96880001_r8,  0.95986998_r8,  0.94862998_r8,  &
        0.93452001_r8,  0.91681999_r8,  0.89459997_r8,  0.86680001_r8,  0.83270001_r8,  &
        0.79189998_r8,  0.74479997_r8,  0.69290000_r8,  0.63839996_r8,  0.58410001_r8,  &
      ! data ((c2(ip,iw),iw=1,30)_r8, ip= 4_r8, 4)/  &
       -0.1808E-06_r8, -0.4642E-06_r8, -0.1045E-05_r8, -0.2058E-05_r8, -0.3581E-05_r8,  &
       -0.5776E-05_r8, -0.8801E-05_r8, -0.1281E-04_r8, -0.1789E-04_r8, -0.2433E-04_r8,  &
       -0.3273E-04_r8, -0.4382E-04_r8, -0.5840E-04_r8, -0.7755E-04_r8, -0.1024E-03_r8,  &
       -0.1342E-03_r8, -0.1737E-03_r8, -0.2217E-03_r8, -0.2791E-03_r8, -0.3473E-03_r8,  &
       -0.4272E-03_r8, -0.5191E-03_r8, -0.6227E-03_r8, -0.7354E-03_r8, -0.8526E-03_r8,  &
       -0.9688E-03_r8, -0.1081E-02_r8, -0.1189E-02_r8, -0.1300E-02_r8, -0.1417E-02_r8,  &
      ! data ((c3(ip,iw),iw=1,30)_r8, ip= 4_r8, 4)/  &
        0.1019E-09_r8,  0.1601E-09_r8, -0.4075E-09_r8, -0.1746E-08_r8, -0.4366E-08_r8,  &
       -0.7960E-08_r8, -0.1294E-07_r8, -0.1898E-07_r8, -0.2899E-07_r8, -0.4594E-07_r8,  &
       -0.7267E-07_r8, -0.1088E-06_r8, -0.1536E-06_r8, -0.2164E-06_r8, -0.3002E-06_r8,  &
       -0.4055E-06_r8, -0.5260E-06_r8, -0.6571E-06_r8, -0.8022E-06_r8, -0.9624E-06_r8,  &
       -0.1139E-05_r8, -0.1315E-05_r8, -0.1456E-05_r8, -0.1512E-05_r8, -0.1420E-05_r8,  &
       -0.1137E-05_r8, -0.6483E-06_r8,  0.6679E-08_r8,  0.7652E-06_r8,  0.1574E-05_r8,  &
      ! data ((c1(ip,iw),iw=1,30)_r8, ip= 5_r8, 5)/  &
        0.99985641_r8,  0.99976403_r8,  0.99963748_r8,  0.99947661_r8,  0.99926913_r8,  &
        0.99898303_r8,  0.99857402_r8,  0.99800003_r8,  0.99723399_r8,  0.99624503_r8,  &
        0.99498397_r8,  0.99335301_r8,  0.99123502_r8,  0.98847997_r8,  0.98488998_r8,  &
        0.98023999_r8,  0.97421998_r8,  0.96648002_r8,  0.95659000_r8,  0.94404000_r8,  &
        0.92815000_r8,  0.90802002_r8,  0.88270003_r8,  0.85119998_r8,  0.81290001_r8,  &
        0.76770002_r8,  0.71679997_r8,  0.66219997_r8,  0.60670000_r8,  0.55250001_r8,  &
      ! data ((c2(ip,iw),iw=1,30)_r8, ip= 5_r8, 5)/  &
       -0.1827E-06_r8, -0.4608E-06_r8, -0.1042E-05_r8, -0.2053E-05_r8, -0.3565E-05_r8,  &
       -0.5745E-05_r8, -0.8758E-05_r8, -0.1273E-04_r8, -0.1778E-04_r8, -0.2417E-04_r8,  &
       -0.3250E-04_r8, -0.4347E-04_r8, -0.5801E-04_r8, -0.7729E-04_r8, -0.1025E-03_r8,  &
       -0.1349E-03_r8, -0.1755E-03_r8, -0.2249E-03_r8, -0.2842E-03_r8, -0.3549E-03_r8,  &
       -0.4380E-03_r8, -0.5340E-03_r8, -0.6428E-03_r8, -0.7613E-03_r8, -0.8854E-03_r8,  &
       -0.1009E-02_r8, -0.1131E-02_r8, -0.1252E-02_r8, -0.1376E-02_r8, -0.1502E-02_r8,  &
      ! data ((c3(ip,iw),iw=1,30)_r8, ip= 5_r8, 5)/  &
        0.4366E-10_r8, -0.1455E-10_r8, -0.4075E-09_r8, -0.1804E-08_r8, -0.4293E-08_r8,  &
       -0.8178E-08_r8, -0.1301E-07_r8, -0.1915E-07_r8, -0.2938E-07_r8, -0.4664E-07_r8,  &
       -0.7365E-07_r8, -0.1090E-06_r8, -0.1539E-06_r8, -0.2158E-06_r8, -0.2992E-06_r8,  &
       -0.4033E-06_r8, -0.5230E-06_r8, -0.6537E-06_r8, -0.7976E-06_r8, -0.9601E-06_r8,  &
       -0.1135E-05_r8, -0.1305E-05_r8, -0.1440E-05_r8, -0.1490E-05_r8, -0.1389E-05_r8,  &
       -0.1087E-05_r8, -0.5646E-06_r8,  0.1475E-06_r8,  0.9852E-06_r8,  0.1853E-05_r8,  &
      ! data ((c1(ip,iw),iw=1,30)_r8, ip= 6_r8, 6)/  &
        0.99985617_r8,  0.99976331_r8,  0.99963629_r8,  0.99947429_r8,  0.99926388_r8,  &
        0.99897301_r8,  0.99855602_r8,  0.99796802_r8,  0.99717802_r8,  0.99614400_r8,  &
        0.99480897_r8,  0.99306899_r8,  0.99078500_r8,  0.98778999_r8,  0.98387998_r8,  &
        0.97876000_r8,  0.97211999_r8,  0.96350002_r8,  0.95240998_r8,  0.93821001_r8,  &
        0.92009002_r8,  0.89709997_r8,  0.86820000_r8,  0.83249998_r8,  0.78970003_r8,  &
        0.74039996_r8,  0.68630004_r8,  0.63010001_r8,  0.57459998_r8,  0.52069998_r8,  &
      ! data ((c2(ip,iw),iw=1,30)_r8, ip= 6_r8, 6)/  &
       -0.1798E-06_r8, -0.4580E-06_r8, -0.1033E-05_r8, -0.2039E-05_r8, -0.3544E-05_r8,  &
       -0.5709E-05_r8, -0.8696E-05_r8, -0.1264E-04_r8, -0.1763E-04_r8, -0.2395E-04_r8,  &
       -0.3220E-04_r8, -0.4311E-04_r8, -0.5777E-04_r8, -0.7732E-04_r8, -0.1032E-03_r8,  &
       -0.1365E-03_r8, -0.1784E-03_r8, -0.2295E-03_r8, -0.2914E-03_r8, -0.3653E-03_r8,  &
       -0.4527E-03_r8, -0.5541E-03_r8, -0.6689E-03_r8, -0.7947E-03_r8, -0.9265E-03_r8,  &
       -0.1060E-02_r8, -0.1192E-02_r8, -0.1326E-02_r8, -0.1460E-02_r8, -0.1586E-02_r8,  &
      ! data ((c3(ip,iw),iw=1,30)_r8, ip= 6_r8, 6)/  &
        0.8731E-10_r8,  0.0000E+00_r8, -0.3492E-09_r8, -0.1892E-08_r8, -0.4322E-08_r8,  &
       -0.8367E-08_r8, -0.1318E-07_r8, -0.1962E-07_r8, -0.3024E-07_r8, -0.4708E-07_r8,  &
       -0.7359E-07_r8, -0.1087E-06_r8, -0.1534E-06_r8, -0.2152E-06_r8, -0.2978E-06_r8,  &
       -0.4008E-06_r8, -0.5207E-06_r8, -0.6509E-06_r8, -0.7968E-06_r8, -0.9584E-06_r8,  &
       -0.1128E-05_r8, -0.1297E-05_r8, -0.1425E-05_r8, -0.1461E-05_r8, -0.1342E-05_r8,  &
       -0.1009E-05_r8, -0.4283E-06_r8,  0.3666E-06_r8,  0.1272E-05_r8,  0.2171E-05_r8,  &
      ! data ((c1(ip,iw),iw=1,30)_r8, ip= 7_r8, 7)/  &
        0.99985600_r8,  0.99976230_r8,  0.99963462_r8,  0.99947017_r8,  0.99925607_r8,  &
        0.99895698_r8,  0.99852800_r8,  0.99791902_r8,  0.99709100_r8,  0.99599499_r8,  &
        0.99456000_r8,  0.99267203_r8,  0.99017102_r8,  0.98688000_r8,  0.98255002_r8,  &
        0.97685999_r8,  0.96941000_r8,  0.95969999_r8,  0.94709998_r8,  0.93085998_r8,  &
        0.91001999_r8,  0.88360000_r8,  0.85060000_r8,  0.81040001_r8,  0.76319999_r8,  &
        0.71029997_r8,  0.65400004_r8,  0.59740001_r8,  0.54229999_r8,  0.48839998_r8,  &
      ! data ((c2(ip,iw),iw=1,30)_r8, ip= 7_r8, 7)/  &
       -0.1784E-06_r8, -0.4551E-06_r8, -0.1023E-05_r8, -0.2019E-05_r8, -0.3507E-05_r8,  &
       -0.5651E-05_r8, -0.8608E-05_r8, -0.1250E-04_r8, -0.1744E-04_r8, -0.2370E-04_r8,  &
       -0.3189E-04_r8, -0.4289E-04_r8, -0.5777E-04_r8, -0.7787E-04_r8, -0.1045E-03_r8,  &
       -0.1392E-03_r8, -0.1828E-03_r8, -0.2365E-03_r8, -0.3015E-03_r8, -0.3797E-03_r8,  &
       -0.4723E-03_r8, -0.5803E-03_r8, -0.7026E-03_r8, -0.8365E-03_r8, -0.9772E-03_r8,  &
       -0.1120E-02_r8, -0.1265E-02_r8, -0.1409E-02_r8, -0.1547E-02_r8, -0.1665E-02_r8,  &
      ! data ((c3(ip,iw),iw=1,30)_r8, ip= 7_r8, 7)/ &
        0.5821E-10_r8,  0.8731E-10_r8, -0.4366E-09_r8, -0.1935E-08_r8, -0.4555E-08_r8, &
       -0.8455E-08_r8, -0.1356E-07_r8, -0.2024E-07_r8, -0.3079E-07_r8, -0.4758E-07_r8, &
       -0.7352E-07_r8, -0.1078E-06_r8, -0.1520E-06_r8, -0.2139E-06_r8, -0.2964E-06_r8, &
       -0.3997E-06_r8, -0.5185E-06_r8, -0.6493E-06_r8, -0.7943E-06_r8, -0.9568E-06_r8, &
       -0.1127E-05_r8, -0.1288E-05_r8, -0.1405E-05_r8, -0.1425E-05_r8, -0.1275E-05_r8, &
       -0.8809E-06_r8, -0.2158E-06_r8,  0.6597E-06_r8,  0.1610E-05_r8,  0.2524E-05_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8, ip= 8_r8, 8)/ &
        0.99985582_r8,  0.99976122_r8,  0.99963123_r8,  0.99946368_r8,  0.99924308_r8, &
        0.99893397_r8,  0.99848598_r8,  0.99784499_r8,  0.99696398_r8,  0.99577999_r8, &
        0.99421299_r8,  0.99212801_r8,  0.98935997_r8,  0.98569000_r8,  0.98083001_r8, &
        0.97442001_r8,  0.96595001_r8,  0.95486999_r8,  0.94040000_r8,  0.92163002_r8, &
        0.89760000_r8,  0.86720002_r8,  0.82969999_r8,  0.78499997_r8,  0.73370004_r8, &
        0.67799997_r8,  0.62070000_r8,  0.56439996_r8,  0.50960004_r8,  0.45539999_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8, ip= 8_r8, 8)/ &
       -0.1760E-06_r8, -0.4451E-06_r8, -0.1004E-05_r8, -0.1989E-05_r8, -0.3457E-05_r8, &
       -0.5574E-05_r8, -0.8470E-05_r8, -0.1230E-04_r8, -0.1721E-04_r8, -0.2344E-04_r8, &
       -0.3168E-04_r8, -0.4286E-04_r8, -0.5815E-04_r8, -0.7898E-04_r8, -0.1070E-03_r8, &
       -0.1434E-03_r8, -0.1892E-03_r8, -0.2460E-03_r8, -0.3152E-03_r8, -0.3985E-03_r8, &
       -0.4981E-03_r8, -0.6139E-03_r8, -0.7448E-03_r8, -0.8878E-03_r8, -0.1038E-02_r8, &
       -0.1193E-02_r8, -0.1348E-02_r8, -0.1499E-02_r8, -0.1631E-02_r8, -0.1735E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30), ip= 8, 8)/ &
       -0.1455E-10_r8, 0.4366E-10_r8,-0.3929E-09_r8,-0.2081E-08_r8,-0.4700E-08_r8, &
       -0.8804E-08_r8,-0.1417E-07_r8,-0.2068E-07_r8,-0.3143E-07_r8,-0.4777E-07_r8, &
       -0.7336E-07_r8,-0.1070E-06_r8,-0.1517E-06_r8,-0.2134E-06_r8,-0.2967E-06_r8, &
       -0.3991E-06_r8,-0.5164E-06_r8,-0.6510E-06_r8,-0.7979E-06_r8,-0.9575E-06_r8, &
       -0.1123E-05_r8,-0.1279E-05_r8,-0.1382E-05_r8,-0.1374E-05_r8,-0.1166E-05_r8, &
       -0.6893E-06_r8, 0.7339E-07_r8, 0.1013E-05_r8, 0.1982E-05_r8, 0.2896E-05_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8,ip= 9_r8,9)/ &
        0.99985498_r8, 0.99975908_r8, 0.99962622_r8, 0.99945402_r8, 0.99922228_r8, &
        0.99889803_r8, 0.99842203_r8, 0.99773699_r8, 0.99677801_r8, 0.99547797_r8, &
        0.99373603_r8, 0.99140298_r8, 0.98829001_r8, 0.98413998_r8, 0.97863001_r8, &
        0.97127002_r8, 0.96156001_r8, 0.94875997_r8, 0.93197000_r8, 0.91017997_r8, &
        0.88230002_r8, 0.84749997_r8, 0.80540001_r8, 0.75620002_r8, 0.70159996_r8, &
        0.64429998_r8, 0.58710003_r8, 0.53130001_r8, 0.47640002_r8, 0.42189997_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8,ip= 9_r8,9)/ &
       -0.1717E-06_r8,-0.4327E-06_r8,-0.9759E-06_r8,-0.1943E-05_r8,-0.3391E-05_r8, &
       -0.5454E-05_r8,-0.8297E-05_r8,-0.1209E-04_r8,-0.1697E-04_r8,-0.2322E-04_r8, &
       -0.3163E-04_r8,-0.4318E-04_r8,-0.5910E-04_r8,-0.8111E-04_r8,-0.1108E-03_r8, &
       -0.1493E-03_r8,-0.1982E-03_r8,-0.2588E-03_r8,-0.3333E-03_r8,-0.4237E-03_r8, &
       -0.5312E-03_r8,-0.6562E-03_r8,-0.7968E-03_r8,-0.9496E-03_r8,-0.1110E-02_r8, &
       -0.1276E-02_r8,-0.1439E-02_r8,-0.1588E-02_r8,-0.1708E-02_r8,-0.1796E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30)_r8,ip= 9_r8,9)/ &
        0.0000E+00_r8, 0.1455E-10_r8,-0.3638E-09_r8,-0.2299E-08_r8,-0.4744E-08_r8, &
       -0.9284E-08_r8,-0.1445E-07_r8,-0.2141E-07_r8,-0.3162E-07_r8,-0.4761E-07_r8, &
       -0.7248E-07_r8,-0.1065E-06_r8,-0.1501E-06_r8,-0.2140E-06_r8,-0.2981E-06_r8, &
       -0.3994E-06_r8,-0.5201E-06_r8,-0.6549E-06_r8,-0.8009E-06_r8,-0.9627E-06_r8, &
       -0.1125E-05_r8,-0.1266E-05_r8,-0.1348E-05_r8,-0.1292E-05_r8,-0.1005E-05_r8, &
       -0.4166E-06_r8, 0.4279E-06_r8, 0.1401E-05_r8, 0.2379E-05_r8, 0.3278E-05_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8,ip=10,10)/ &
        0.99985462_r8, 0.99975640_r8, 0.99961889_r8, 0.99943668_r8, 0.99919188_r8, &
        0.99884301_r8, 0.99832898_r8, 0.99757999_r8, 0.99651998_r8, 0.99506402_r8, &
        0.99309200_r8, 0.99044400_r8, 0.98689002_r8, 0.98215997_r8, 0.97579002_r8, &
        0.96730000_r8, 0.95603001_r8, 0.94110000_r8, 0.92149001_r8, 0.89609998_r8, &
        0.86399996_r8, 0.82449996_r8, 0.77759999_r8, 0.72459996_r8, 0.66769999_r8, &
        0.61000001_r8, 0.55340004_r8, 0.49769998_r8, 0.44250000_r8, 0.38810003_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8,ip=10,10)/ &
       -0.1607E-06_r8,-0.4160E-06_r8,-0.9320E-06_r8,-0.1872E-05_r8,-0.3281E-05_r8, &
       -0.5286E-05_r8,-0.8097E-05_r8,-0.1187E-04_r8,-0.1677E-04_r8,-0.2320E-04_r8, &
       -0.3190E-04_r8,-0.4402E-04_r8,-0.6081E-04_r8,-0.8441E-04_r8,-0.1162E-03_r8, &
       -0.1576E-03_r8,-0.2102E-03_r8,-0.2760E-03_r8,-0.3571E-03_r8,-0.4558E-03_r8, &
       -0.5730E-03_r8,-0.7082E-03_r8,-0.8591E-03_r8,-0.1022E-02_r8,-0.1194E-02_r8, &
       -0.1368E-02_r8,-0.1533E-02_r8,-0.1671E-02_r8,-0.1775E-02_r8,-0.1843E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30)_r8,ip=10,10)/ &
       -0.1164E-09_r8,-0.7276E-10_r8,-0.5530E-09_r8,-0.2270E-08_r8,-0.5093E-08_r8, &
       -0.9517E-08_r8,-0.1502E-07_r8,-0.2219E-07_r8,-0.3171E-07_r8,-0.4712E-07_r8, &
       -0.7123E-07_r8,-0.1042E-06_r8,-0.1493E-06_r8,-0.2156E-06_r8,-0.2999E-06_r8, &
       -0.4027E-06_r8,-0.5243E-06_r8,-0.6616E-06_r8,-0.8125E-06_r8,-0.9691E-06_r8, &
       -0.1126E-05_r8,-0.1251E-05_r8,-0.1294E-05_r8,-0.1163E-05_r8,-0.7639E-06_r8, &
       -0.7395E-07_r8, 0.8279E-06_r8, 0.1819E-05_r8, 0.2795E-05_r8, 0.3647E-05_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8,ip=11,11)/ &
        0.99985212_r8, 0.99975210_r8, 0.99960798_r8, 0.99941242_r8, 0.99914628_r8, &
        0.99876302_r8, 0.99819702_r8, 0.99736100_r8, 0.99616700_r8, 0.99450397_r8, &
        0.99225003_r8, 0.98920000_r8, 0.98510998_r8, 0.97961998_r8, 0.97220999_r8, &
        0.96231002_r8, 0.94909000_r8, 0.93155003_r8, 0.90856999_r8, 0.87910002_r8, &
        0.84219998_r8, 0.79790002_r8, 0.74669999_r8, 0.69080001_r8, 0.63300002_r8, &
        0.57570004_r8, 0.51950002_r8, 0.46359998_r8, 0.40829998_r8, 0.35450000_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8,ip=11,11)/ &
       -0.1531E-06_r8,-0.3864E-06_r8,-0.8804E-06_r8,-0.1776E-05_r8,-0.3131E-05_r8, &
       -0.5082E-05_r8,-0.7849E-05_r8,-0.1164E-04_r8,-0.1669E-04_r8,-0.2340E-04_r8, &
       -0.3261E-04_r8,-0.4546E-04_r8,-0.6380E-04_r8,-0.8932E-04_r8,-0.1237E-03_r8, &
       -0.1687E-03_r8,-0.2262E-03_r8,-0.2984E-03_r8,-0.3880E-03_r8,-0.4964E-03_r8, &
       -0.6244E-03_r8,-0.7705E-03_r8,-0.9325E-03_r8,-0.1107E-02_r8,-0.1288E-02_r8, &
       -0.1466E-02_r8,-0.1623E-02_r8,-0.1746E-02_r8,-0.1831E-02_r8,-0.1875E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30)_r8,ip=11,11)/ &
        0.1019E-09_r8,-0.2037E-09_r8,-0.8004E-09_r8,-0.2387E-08_r8,-0.5326E-08_r8, &
       -0.9764E-08_r8,-0.1576E-07_r8,-0.2256E-07_r8,-0.3180E-07_r8,-0.4616E-07_r8, &
       -0.7026E-07_r8,-0.1031E-06_r8,-0.1520E-06_r8,-0.2181E-06_r8,-0.3037E-06_r8, &
       -0.4109E-06_r8,-0.5354E-06_r8,-0.6740E-06_r8,-0.8241E-06_r8,-0.9810E-06_r8, &
       -0.1126E-05_r8,-0.1221E-05_r8,-0.1200E-05_r8,-0.9678E-06_r8,-0.4500E-06_r8, &
        0.3236E-06_r8, 0.1256E-05_r8, 0.2259E-05_r8, 0.3206E-05_r8, 0.3978E-05_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8,ip=12,12)/ &
        0.99985027_r8, 0.99974507_r8, 0.99959022_r8, 0.99937689_r8, 0.99907988_r8, &
        0.99865198_r8, 0.99801201_r8, 0.99706602_r8, 0.99569201_r8, 0.99377203_r8, &
        0.99115402_r8, 0.98762000_r8, 0.98286003_r8, 0.97640002_r8, 0.96771997_r8, &
        0.95604998_r8, 0.94045001_r8, 0.91979003_r8, 0.89289999_r8, 0.85879999_r8, &
        0.81700003_r8, 0.76800001_r8, 0.71340001_r8, 0.65579998_r8, 0.59810001_r8, &
        0.54139996_r8, 0.48519999_r8, 0.42909998_r8, 0.37410003_r8, 0.32190001_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8,ip=12,12)/ &
       -0.1340E-06_r8,-0.3478E-06_r8,-0.8189E-06_r8,-0.1653E-05_r8,-0.2944E-05_r8, &
       -0.4852E-05_r8,-0.7603E-05_r8,-0.1150E-04_r8,-0.1682E-04_r8,-0.2400E-04_r8, &
       -0.3390E-04_r8,-0.4799E-04_r8,-0.6807E-04_r8,-0.9596E-04_r8,-0.1338E-03_r8, &
       -0.1833E-03_r8,-0.2471E-03_r8,-0.3275E-03_r8,-0.4268E-03_r8,-0.5466E-03_r8, &
       -0.6862E-03_r8,-0.8439E-03_r8,-0.1017E-02_r8,-0.1201E-02_r8,-0.1389E-02_r8, &
       -0.1563E-02_r8,-0.1706E-02_r8,-0.1809E-02_r8,-0.1872E-02_r8,-0.1890E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30)_r8,ip=12,12)/ &
       -0.1455E-10_r8,-0.1892E-09_r8,-0.8295E-09_r8,-0.2547E-08_r8,-0.5544E-08_r8, &
       -0.1014E-07_r8,-0.1605E-07_r8,-0.2341E-07_r8,-0.3156E-07_r8,-0.4547E-07_r8, &
       -0.6749E-07_r8,-0.1034E-06_r8,-0.1550E-06_r8,-0.2230E-06_r8,-0.3130E-06_r8, &
       -0.4219E-06_r8,-0.5469E-06_r8,-0.6922E-06_r8,-0.8448E-06_r8,-0.9937E-06_r8, &
       -0.1118E-05_r8,-0.1166E-05_r8,-0.1054E-05_r8,-0.6926E-06_r8,-0.7180E-07_r8, &
        0.7515E-06_r8, 0.1709E-05_r8, 0.2703E-05_r8, 0.3593E-05_r8, 0.4232E-05_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8,ip=13,13)/ &
        0.99984729_r8, 0.99973530_r8, 0.99956691_r8, 0.99932659_r8, 0.99898797_r8, &
        0.99849701_r8, 0.99776399_r8, 0.99667102_r8, 0.99507397_r8, 0.99283201_r8, &
        0.98977000_r8, 0.98563999_r8, 0.98001999_r8, 0.97241002_r8, 0.96213001_r8, &
        0.94830000_r8, 0.92980999_r8, 0.90546000_r8, 0.87409997_r8, 0.83510000_r8, &
        0.78850001_r8, 0.73549998_r8, 0.67850000_r8, 0.62049997_r8, 0.56340003_r8, &
        0.50699997_r8, 0.45050001_r8, 0.39450002_r8, 0.34060001_r8, 0.29079998_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8,ip=13,13)/ &
       -0.1163E-06_r8,-0.3048E-06_r8,-0.7186E-06_r8,-0.1495E-05_r8,-0.2726E-05_r8, &
       -0.4588E-05_r8,-0.7396E-05_r8,-0.1152E-04_r8,-0.1725E-04_r8,-0.2514E-04_r8, &
       -0.3599E-04_r8,-0.5172E-04_r8,-0.7403E-04_r8,-0.1051E-03_r8,-0.1469E-03_r8, &
       -0.2023E-03_r8,-0.2735E-03_r8,-0.3637E-03_r8,-0.4746E-03_r8,-0.6067E-03_r8, &
       -0.7586E-03_r8,-0.9281E-03_r8,-0.1112E-02_r8,-0.1304E-02_r8,-0.1491E-02_r8, &
       -0.1653E-02_r8,-0.1777E-02_r8,-0.1860E-02_r8,-0.1896E-02_r8,-0.1891E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30)_r8,ip=13,13)/  &
       -0.1455E-09_r8,-0.2765E-09_r8,-0.9750E-09_r8,-0.2794E-08_r8,-0.5413E-08_r8, &
       -0.1048E-07_r8,-0.1625E-07_r8,-0.2344E-07_r8,-0.3105E-07_r8,-0.4304E-07_r8, &
       -0.6608E-07_r8,-0.1057E-06_r8,-0.1587E-06_r8,-0.2308E-06_r8,-0.3235E-06_r8, &
       -0.4373E-06_r8,-0.5687E-06_r8,-0.7156E-06_r8,-0.8684E-06_r8,-0.1007E-05_r8, &
       -0.1094E-05_r8,-0.1062E-05_r8,-0.8273E-06_r8,-0.3485E-06_r8, 0.3463E-06_r8, &
        0.1206E-05_r8, 0.2173E-05_r8, 0.3132E-05_r8, 0.3919E-05_r8, 0.4370E-05_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8,ip=14,14)/  &
        0.99984348_r8, 0.99972272_r8, 0.99953479_r8, 0.99926043_r8, 0.99886698_r8, &
        0.99829400_r8, 0.99744201_r8, 0.99615997_r8, 0.99429500_r8, 0.99166000_r8, &
        0.98806000_r8, 0.98316997_r8, 0.97649997_r8, 0.96748000_r8, 0.95525998_r8, &
        0.93878001_r8, 0.91687000_r8, 0.88830000_r8, 0.85220003_r8, 0.80820000_r8, &
        0.75699997_r8, 0.70099998_r8, 0.64300001_r8, 0.58550000_r8, 0.52890003_r8, &
        0.47219998_r8, 0.41560000_r8, 0.36040002_r8, 0.30849999_r8, 0.26169997_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8,ip=14,14)/  &
       -0.8581E-07_r8,-0.2557E-06_r8,-0.6103E-06_r8,-0.1305E-05_r8,-0.2472E-05_r8, &
       -0.4334E-05_r8,-0.7233E-05_r8,-0.1167E-04_r8,-0.1806E-04_r8,-0.2679E-04_r8, &
       -0.3933E-04_r8,-0.5705E-04_r8,-0.8194E-04_r8,-0.1165E-03_r8,-0.1637E-03_r8, &
       -0.2259E-03_r8,-0.3068E-03_r8,-0.4082E-03_r8,-0.5318E-03_r8,-0.6769E-03_r8, &  
       -0.8415E-03_r8,-0.1023E-02_r8,-0.1216E-02_r8,-0.1410E-02_r8,-0.1588E-02_r8, &
       -0.1733E-02_r8,-0.1837E-02_r8,-0.1894E-02_r8,-0.1904E-02_r8,-0.1881E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30)_r8,ip=14,14)/  &
       -0.2037E-09_r8,-0.4220E-09_r8,-0.1091E-08_r8,-0.2896E-08_r8,-0.5821E-08_r8, &
       -0.1052E-07_r8,-0.1687E-07_r8,-0.2353E-07_r8,-0.3193E-07_r8,-0.4254E-07_r8, &
       -0.6685E-07_r8,-0.1072E-06_r8,-0.1638E-06_r8,-0.2427E-06_r8,-0.3421E-06_r8, &
       -0.4600E-06_r8,-0.5946E-06_r8,-0.7472E-06_r8,-0.8958E-06_r8,-0.1009E-05_r8, &
       -0.1032E-05_r8,-0.8919E-06_r8,-0.5224E-06_r8, 0.5218E-07_r8, 0.7886E-06_r8, &
        0.1672E-05_r8, 0.2626E-05_r8, 0.3513E-05_r8, 0.4138E-05_r8, 0.4379E-05_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8,ip=15,15)/  &
        0.99983788_r8, 0.99970680_r8, 0.99949580_r8, 0.99917668_r8, 0.99871200_r8, &
        0.99803603_r8, 0.99703097_r8, 0.99552703_r8, 0.99333203_r8, 0.99023402_r8, &
        0.98597997_r8, 0.98013997_r8, 0.97223002_r8, 0.96145999_r8, 0.94686002_r8, &
        0.92727000_r8, 0.90142000_r8, 0.86820000_r8, 0.82700002_r8, 0.77820003_r8, &  
        0.72350001_r8, 0.66569996_r8, 0.60769999_r8, 0.55089998_r8, 0.49430001_r8, &
        0.43739998_r8, 0.38110000_r8, 0.32749999_r8, 0.27840000_r8, 0.23479998_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8,ip=15,15)/  &
       -0.8246E-07_r8,-0.2070E-06_r8,-0.4895E-06_r8,-0.1106E-05_r8,-0.2216E-05_r8, &
       -0.4077E-05_r8,-0.7150E-05_r8,-0.1202E-04_r8,-0.1920E-04_r8,-0.2938E-04_r8, &
       -0.4380E-04_r8,-0.6390E-04_r8,-0.9209E-04_r8,-0.1310E-03_r8,-0.1843E-03_r8, &
       -0.2554E-03_r8,-0.3468E-03_r8,-0.4611E-03_r8,-0.5982E-03_r8,-0.7568E-03_r8, &
       -0.9340E-03_r8,-0.1126E-02_r8,-0.1324E-02_r8,-0.1514E-02_r8,-0.1676E-02_r8, &
       -0.1801E-02_r8,-0.1881E-02_r8,-0.1911E-02_r8,-0.1900E-02_r8,-0.1867E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30)_r8,ip=15,15)/  &
       -0.1601E-09_r8,-0.3492E-09_r8,-0.1019E-08_r8,-0.2634E-08_r8,-0.5632E-08_r8, &
       -0.1065E-07_r8,-0.1746E-07_r8,-0.2542E-07_r8,-0.3206E-07_r8,-0.4390E-07_r8, &
       -0.6956E-07_r8,-0.1093E-06_r8,-0.1729E-06_r8,-0.2573E-06_r8,-0.3612E-06_r8, &
       -0.4904E-06_r8,-0.6342E-06_r8,-0.7834E-06_r8,-0.9175E-06_r8,-0.9869E-06_r8, &
       -0.9164E-06_r8,-0.6386E-06_r8,-0.1544E-06_r8, 0.4798E-06_r8, 0.1252E-05_r8, &
        0.2137E-05_r8, 0.3043E-05_r8, 0.3796E-05_r8, 0.4211E-05_r8, 0.4332E-05_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8,ip=16,16)/  &
        0.99983227_r8, 0.99968958_r8, 0.99945217_r8, 0.99907941_r8, 0.99852598_r8, &
        0.99772000_r8, 0.99652398_r8, 0.99475902_r8, 0.99218899_r8, 0.98856002_r8, &
        0.98348999_r8, 0.97653997_r8, 0.96708000_r8, 0.95420998_r8, 0.93677002_r8, &
        0.91352999_r8, 0.88330001_r8, 0.84509999_r8, 0.79900002_r8, 0.74599999_r8, &
        0.68879998_r8, 0.63049996_r8, 0.57319999_r8, 0.51660001_r8, 0.45969999_r8, &
        0.40289998_r8, 0.34780002_r8, 0.29650003_r8, 0.25070000_r8, 0.20959997_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8,ip=16,16)/  &
       -0.7004E-07_r8,-0.1592E-06_r8,-0.3936E-06_r8,-0.9145E-06_r8,-0.1958E-05_r8, &
       -0.3850E-05_r8,-0.7093E-05_r8,-0.1252E-04_r8,-0.2066E-04_r8,-0.3271E-04_r8, &
       -0.4951E-04_r8,-0.7268E-04_r8,-0.1045E-03_r8,-0.1487E-03_r8,-0.2092E-03_r8, &
       -0.2899E-03_r8,-0.3936E-03_r8,-0.5215E-03_r8,-0.6729E-03_r8,-0.8454E-03_r8, &
       -0.1035E-02_r8,-0.1235E-02_r8,-0.1432E-02_r8,-0.1608E-02_r8,-0.1751E-02_r8, &
       -0.1854E-02_r8,-0.1907E-02_r8,-0.1913E-02_r8,-0.1888E-02_r8,-0.1857E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30)_r8,ip=16,16)/  &
       -0.2328E-09_r8,-0.3347E-09_r8,-0.9750E-09_r8,-0.2314E-08_r8,-0.5166E-08_r8, &
       -0.1052E-07_r8,-0.1726E-07_r8,-0.2605E-07_r8,-0.3532E-07_r8,-0.4949E-07_r8, &
       -0.7229E-07_r8,-0.1133E-06_r8,-0.1799E-06_r8,-0.2725E-06_r8,-0.3881E-06_r8, &
       -0.5249E-06_r8,-0.6763E-06_r8,-0.8227E-06_r8,-0.9279E-06_r8,-0.9205E-06_r8, &
       -0.7228E-06_r8,-0.3109E-06_r8, 0.2583E-06_r8, 0.9390E-06_r8, 0.1726E-05_r8, &
        0.2579E-05_r8, 0.3376E-05_r8, 0.3931E-05_r8, 0.4161E-05_r8, 0.4369E-05_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8,ip=17,17)/  &
        0.99982637_r8, 0.99967217_r8, 0.99940813_r8, 0.99897701_r8, 0.99831802_r8, &
        0.99734300_r8, 0.99592501_r8, 0.99385202_r8, 0.99086499_r8, 0.98659998_r8, &
        0.98057997_r8, 0.97229999_r8, 0.96098000_r8, 0.94555002_r8, 0.92479002_r8, &
        0.89740002_r8, 0.86240000_r8, 0.81919998_r8, 0.76859999_r8, 0.71249998_r8, &
        0.65419996_r8, 0.59630001_r8, 0.53950000_r8, 0.48259997_r8, 0.42549998_r8, &
        0.36940002_r8, 0.31629997_r8, 0.26810002_r8, 0.22520000_r8, 0.18580002_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8,ip=17,17)/  &
       -0.6526E-07_r8,-0.1282E-06_r8,-0.3076E-06_r8,-0.7454E-06_r8,-0.1685E-05_r8, &
       -0.3600E-05_r8,-0.7071E-05_r8,-0.1292E-04_r8,-0.2250E-04_r8,-0.3665E-04_r8, &
       -0.5623E-04_r8,-0.8295E-04_r8,-0.1195E-03_r8,-0.1696E-03_r8,-0.2385E-03_r8, &
       -0.3298E-03_r8,-0.4465E-03_r8,-0.5887E-03_r8,-0.7546E-03_r8,-0.9408E-03_r8, &
       -0.1141E-02_r8,-0.1345E-02_r8,-0.1533E-02_r8,-0.1691E-02_r8,-0.1813E-02_r8, &
       -0.1889E-02_r8,-0.1916E-02_r8,-0.1904E-02_r8,-0.1877E-02_r8,-0.1850E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30)_r8,ip=17,17)/  &
       -0.1746E-09_r8,-0.2037E-09_r8,-0.8149E-09_r8,-0.2095E-08_r8,-0.4889E-08_r8, &
       -0.9517E-08_r8,-0.1759E-07_r8,-0.2740E-07_r8,-0.4147E-07_r8,-0.5774E-07_r8, &
       -0.7909E-07_r8,-0.1199E-06_r8,-0.1877E-06_r8,-0.2859E-06_r8,-0.4137E-06_r8, &
       -0.5649E-06_r8,-0.7218E-06_r8,-0.8516E-06_r8,-0.9022E-06_r8,-0.7905E-06_r8, &
       -0.4531E-06_r8, 0.6917E-07_r8, 0.7009E-06_r8, 0.1416E-05_r8, 0.2194E-05_r8, &
        0.2963E-05_r8, 0.3578E-05_r8, 0.3900E-05_r8, 0.4094E-05_r8, 0.4642E-05_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8,ip=18,18)/  &
        0.99982101_r8, 0.99965781_r8, 0.99936712_r8, 0.99887502_r8, 0.99809802_r8, &
        0.99692702_r8, 0.99523401_r8, 0.99281400_r8, 0.98935997_r8, 0.98435003_r8, &
        0.97728002_r8, 0.96740997_r8, 0.95381999_r8, 0.93539000_r8, 0.91082001_r8, &
        0.87889999_r8, 0.83889997_r8, 0.79100001_r8, 0.73660004_r8, 0.67879999_r8, &
        0.62049997_r8, 0.56330001_r8, 0.50629997_r8, 0.44900000_r8, 0.39209998_r8, &
        0.33749998_r8, 0.28729999_r8, 0.24229997_r8, 0.20150000_r8, 0.16280001_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8,ip=18,18)/  &
       -0.6477E-07_r8,-0.1243E-06_r8,-0.2536E-06_r8,-0.6173E-06_r8,-0.1495E-05_r8, &
       -0.3353E-05_r8,-0.6919E-05_r8,-0.1337E-04_r8,-0.2418E-04_r8,-0.4049E-04_r8, &
       -0.6354E-04_r8,-0.9455E-04_r8,-0.1367E-03_r8,-0.1942E-03_r8,-0.2717E-03_r8, &
       -0.3744E-03_r8,-0.5042E-03_r8,-0.6609E-03_r8,-0.8416E-03_r8,-0.1041E-02_r8, &
       -0.1249E-02_r8,-0.1448E-02_r8,-0.1622E-02_r8,-0.1760E-02_r8,-0.1857E-02_r8, &
       -0.1906E-02_r8,-0.1911E-02_r8,-0.1892E-02_r8,-0.1870E-02_r8,-0.1844E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30)_r8,ip=18,18)/  &
       -0.5821E-10_r8,-0.2328E-09_r8,-0.6985E-09_r8,-0.1368E-08_r8,-0.4351E-08_r8, &
       -0.8993E-08_r8,-0.1579E-07_r8,-0.2916E-07_r8,-0.4904E-07_r8,-0.7010E-07_r8, &
       -0.9623E-07_r8,-0.1332E-06_r8,-0.1928E-06_r8,-0.2977E-06_r8,-0.4371E-06_r8, &
       -0.5992E-06_r8,-0.7586E-06_r8,-0.8580E-06_r8,-0.8238E-06_r8,-0.5811E-06_r8, &
       -0.1298E-06_r8, 0.4702E-06_r8, 0.1162E-05_r8, 0.1905E-05_r8, 0.2632E-05_r8, &
        0.3247E-05_r8, 0.3609E-05_r8, 0.3772E-05_r8, 0.4166E-05_r8, 0.5232E-05_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8,ip=19,19)/  &
        0.99981648_r8, 0.99964571_r8, 0.99933147_r8, 0.99878597_r8, 0.99787998_r8, &
        0.99649400_r8, 0.99448699_r8, 0.99166602_r8, 0.98762000_r8, 0.98181999_r8, &
        0.97352999_r8, 0.96183002_r8, 0.94558001_r8, 0.92363000_r8, 0.89480001_r8, &
        0.85799998_r8, 0.81309998_r8, 0.76100004_r8, 0.70420003_r8, 0.64590001_r8, &
        0.58840001_r8, 0.53139997_r8, 0.47380000_r8, 0.41619998_r8, 0.36030000_r8, &
        0.30809999_r8, 0.26109999_r8, 0.21880001_r8, 0.17909998_r8, 0.14080000_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8,ip=19,19)/  &
       -0.7906E-07_r8,-0.1291E-06_r8,-0.2430E-06_r8,-0.5145E-06_r8,-0.1327E-05_r8, &
       -0.3103E-05_r8,-0.6710E-05_r8,-0.1371E-04_r8,-0.2561E-04_r8,-0.4405E-04_r8, &
       -0.7051E-04_r8,-0.1070E-03_r8,-0.1560E-03_r8,-0.2217E-03_r8,-0.3090E-03_r8, &
       -0.4228E-03_r8,-0.5657E-03_r8,-0.7371E-03_r8,-0.9322E-03_r8,-0.1142E-02_r8, &
       -0.1352E-02_r8,-0.1541E-02_r8,-0.1697E-02_r8,-0.1813E-02_r8,-0.1883E-02_r8, &
       -0.1906E-02_r8,-0.1898E-02_r8,-0.1882E-02_r8,-0.1866E-02_r8,-0.1832E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30)_r8,ip=19,19)/  &
        0.2910E-10_r8, 0.1455E-10_r8,-0.2765E-09_r8,-0.1426E-08_r8,-0.2576E-08_r8, &
       -0.5923E-08_r8,-0.1429E-07_r8,-0.3159E-07_r8,-0.5441E-07_r8,-0.8367E-07_r8, &
       -0.1161E-06_r8,-0.1526E-06_r8,-0.2060E-06_r8,-0.3007E-06_r8,-0.4450E-06_r8, &
       -0.6182E-06_r8,-0.7683E-06_r8,-0.8170E-06_r8,-0.6754E-06_r8,-0.3122E-06_r8, &
        0.2234E-06_r8, 0.8828E-06_r8, 0.1632E-05_r8, 0.2373E-05_r8, 0.3002E-05_r8, &
        0.3384E-05_r8, 0.3499E-05_r8, 0.3697E-05_r8, 0.4517E-05_r8, 0.6117E-05_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8,ip=20,20)/  &
        0.99981302_r8, 0.99963689_r8, 0.99930489_r8, 0.99870700_r8, 0.99768901_r8, &
        0.99608499_r8, 0.99373102_r8, 0.99039900_r8, 0.98566997_r8, 0.97895002_r8, &
        0.96930999_r8, 0.95548999_r8, 0.93621999_r8, 0.91029000_r8, 0.87669998_r8, &
        0.83490002_r8, 0.78549999_r8, 0.73019999_r8, 0.67240000_r8, 0.61469996_r8, &
        0.55779999_r8, 0.50029999_r8, 0.44220001_r8, 0.38489997_r8, 0.33069998_r8, &
        0.28149998_r8, 0.23760003_r8, 0.19690001_r8, 0.15759999_r8, 0.11989999_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8,ip=20,20)/  &
       -0.7762E-07_r8,-0.1319E-06_r8,-0.2315E-06_r8,-0.4780E-06_r8,-0.1187E-05_r8, &
       -0.2750E-05_r8,-0.6545E-05_r8,-0.1393E-04_r8,-0.2645E-04_r8,-0.4652E-04_r8, &
       -0.7657E-04_r8,-0.1190E-03_r8,-0.1766E-03_r8,-0.2520E-03_r8,-0.3499E-03_r8, &
       -0.4751E-03_r8,-0.6307E-03_r8,-0.8160E-03_r8,-0.1024E-02_r8,-0.1240E-02_r8, &
       -0.1443E-02_r8,-0.1619E-02_r8,-0.1757E-02_r8,-0.1849E-02_r8,-0.1892E-02_r8, &
       -0.1896E-02_r8,-0.1886E-02_r8,-0.1878E-02_r8,-0.1861E-02_r8,-0.1807E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30)_r8,ip=20,20)/ &
        0.8731E-10_r8,-0.7276E-10_r8,-0.2328E-09_r8,-0.6403E-09_r8,-0.1455E-08_r8, &
       -0.3827E-08_r8,-0.1270E-07_r8,-0.3014E-07_r8,-0.5594E-07_r8,-0.9677E-07_r8, &
       -0.1422E-06_r8,-0.1823E-06_r8,-0.2296E-06_r8,-0.3094E-06_r8,-0.4399E-06_r8, &
       -0.6008E-06_r8,-0.7239E-06_r8,-0.7014E-06_r8,-0.4562E-06_r8,-0.7778E-08_r8, &
        0.5785E-06_r8, 0.1291E-05_r8, 0.2072E-05_r8, 0.2783E-05_r8, 0.3247E-05_r8, &
        0.3358E-05_r8, 0.3364E-05_r8, 0.3847E-05_r8, 0.5194E-05_r8, 0.7206E-05_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8,ip=21,21)/ &
        0.99981070_r8, 0.99962878_r8, 0.99928439_r8, 0.99864298_r8, 0.99752903_r8, &
        0.99573100_r8, 0.99301797_r8, 0.98905998_r8, 0.98354000_r8, 0.97570997_r8, &
        0.96449000_r8, 0.94837999_r8, 0.92576003_r8, 0.89539999_r8, 0.85680002_r8, &
        0.81000000_r8, 0.75660002_r8, 0.69949996_r8, 0.64199996_r8, 0.58529997_r8, &
        0.52829999_r8, 0.47020000_r8, 0.41200000_r8, 0.35570002_r8, 0.30400002_r8, &
        0.25800002_r8, 0.21609998_r8, 0.17610002_r8, 0.13709998_r8, 0.10020000_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8,ip=21,21)/ &
       -0.1010E-06_r8,-0.1533E-06_r8,-0.2347E-06_r8,-0.4535E-06_r8,-0.1029E-05_r8, &
       -0.2530E-05_r8,-0.6335E-05_r8,-0.1381E-04_r8,-0.2681E-04_r8,-0.4777E-04_r8, &
       -0.8083E-04_r8,-0.1296E-03_r8,-0.1966E-03_r8,-0.2836E-03_r8,-0.3937E-03_r8, &
       -0.5313E-03_r8,-0.6995E-03_r8,-0.8972E-03_r8,-0.1113E-02_r8,-0.1327E-02_r8, &
       -0.1520E-02_r8,-0.1681E-02_r8,-0.1800E-02_r8,-0.1867E-02_r8,-0.1887E-02_r8, &
       -0.1884E-02_r8,-0.1881E-02_r8,-0.1879E-02_r8,-0.1849E-02_r8,-0.1764E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30)_r8,ip=21,21)/ &
        0.8731E-10_r8, 0.1310E-09_r8,-0.2474E-09_r8,-0.2619E-09_r8, 0.8295E-09_r8, &
       -0.1979E-08_r8,-0.1141E-07_r8,-0.2621E-07_r8,-0.5799E-07_r8,-0.1060E-06_r8, &
       -0.1621E-06_r8,-0.2281E-06_r8,-0.2793E-06_r8,-0.3335E-06_r8,-0.4277E-06_r8, &
       -0.5429E-06_r8,-0.5970E-06_r8,-0.4872E-06_r8,-0.1775E-06_r8, 0.3028E-06_r8, &
        0.9323E-06_r8, 0.1680E-05_r8, 0.2452E-05_r8, 0.3063E-05_r8, 0.3299E-05_r8, &
        0.3219E-05_r8, 0.3369E-05_r8, 0.4332E-05_r8, 0.6152E-05_r8, 0.8413E-05_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8,ip=22,22)/ &
        0.99980962_r8, 0.99962330_r8, 0.99926400_r8, 0.99858999_r8, 0.99741602_r8, &
        0.99547201_r8, 0.99236798_r8, 0.98776001_r8, 0.98124999_r8, 0.97210997_r8, &
        0.95902997_r8, 0.94033003_r8, 0.91415000_r8, 0.87919998_r8, 0.83529997_r8, &
        0.78380001_r8, 0.72749996_r8, 0.66990000_r8, 0.61339998_r8, 0.55720001_r8, &
        0.49980003_r8, 0.44129997_r8, 0.38360000_r8, 0.32929999_r8, 0.28070003_r8, &
        0.23710001_r8, 0.19620001_r8, 0.15619999_r8, 0.11769998_r8, 0.08200002_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8,ip=22,22)/ &
       -0.1258E-06_r8,-0.1605E-06_r8,-0.2581E-06_r8,-0.4286E-06_r8,-0.8321E-06_r8, &
       -0.2392E-05_r8,-0.6163E-05_r8,-0.1358E-04_r8,-0.2646E-04_r8,-0.4792E-04_r8, &
       -0.8284E-04_r8,-0.1369E-03_r8,-0.2138E-03_r8,-0.3141E-03_r8,-0.4393E-03_r8, &
       -0.5917E-03_r8,-0.7731E-03_r8,-0.9796E-03_r8,-0.1195E-02_r8,-0.1399E-02_r8, &
       -0.1579E-02_r8,-0.1725E-02_r8,-0.1822E-02_r8,-0.1867E-02_r8,-0.1877E-02_r8, &
       -0.1879E-02_r8,-0.1886E-02_r8,-0.1879E-02_r8,-0.1825E-02_r8,-0.1706E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30)_r8,ip=22,22)/ &
       -0.8731E-10_r8, 0.2910E-10_r8, 0.7276E-10_r8, 0.1281E-08_r8, 0.1222E-08_r8, &
       -0.1935E-08_r8,-0.8004E-08_r8,-0.2258E-07_r8,-0.5428E-07_r8,-0.1085E-06_r8, &
       -0.1835E-06_r8,-0.2716E-06_r8,-0.3446E-06_r8,-0.3889E-06_r8,-0.4203E-06_r8, &
       -0.4394E-06_r8,-0.3716E-06_r8,-0.1677E-06_r8, 0.1622E-06_r8, 0.6327E-06_r8, &
        0.1275E-05_r8, 0.2018E-05_r8, 0.2716E-05_r8, 0.3137E-05_r8, 0.3136E-05_r8, &
        0.3078E-05_r8, 0.3649E-05_r8, 0.5152E-05_r8, 0.7315E-05_r8, 0.9675E-05_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8,ip=23,23)/ &
        0.99980921_r8, 0.99961692_r8, 0.99924570_r8, 0.99854898_r8, 0.99734801_r8, &
        0.99527103_r8, 0.99182302_r8, 0.98655999_r8, 0.97895002_r8, 0.96814001_r8, &
        0.95284998_r8, 0.93124998_r8, 0.90130001_r8, 0.86170000_r8, 0.81290001_r8, &
        0.75740004_r8, 0.69920003_r8, 0.64199996_r8, 0.58640003_r8, 0.53020000_r8, &
        0.47240001_r8, 0.41399997_r8, 0.35780001_r8, 0.30650002_r8, 0.26069999_r8, &
        0.21850002_r8, 0.17750001_r8, 0.13739997_r8, 0.09950000_r8, 0.06540000_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8,ip=23,23)/ &
       -0.1434E-06_r8,-0.1676E-06_r8,-0.2699E-06_r8,-0.2859E-06_r8,-0.7542E-06_r8, &
       -0.2273E-05_r8,-0.5898E-05_r8,-0.1292E-04_r8,-0.2538E-04_r8,-0.4649E-04_r8, &
       -0.8261E-04_r8,-0.1405E-03_r8,-0.2259E-03_r8,-0.3407E-03_r8,-0.4845E-03_r8, &
       -0.6561E-03_r8,-0.8524E-03_r8,-0.1062E-02_r8,-0.1266E-02_r8,-0.1456E-02_r8, &
       -0.1621E-02_r8,-0.1748E-02_r8,-0.1823E-02_r8,-0.1854E-02_r8,-0.1868E-02_r8, &
       -0.1886E-02_r8,-0.1899E-02_r8,-0.1876E-02_r8,-0.1790E-02_r8,-0.1636E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30)_r8,ip=23,23)/ &
       -0.1892E-09_r8,-0.2474E-09_r8, 0.1892E-09_r8, 0.2561E-08_r8, 0.4366E-09_r8, &
       -0.1499E-08_r8,-0.4336E-08_r8,-0.1740E-07_r8,-0.5233E-07_r8,-0.1055E-06_r8, &
       -0.1940E-06_r8,-0.3113E-06_r8,-0.4161E-06_r8,-0.4620E-06_r8,-0.4316E-06_r8, &
       -0.3031E-06_r8,-0.5438E-07_r8, 0.2572E-06_r8, 0.5773E-06_r8, 0.1008E-05_r8, &
        0.1609E-05_r8, 0.2290E-05_r8, 0.2817E-05_r8, 0.2940E-05_r8, 0.2803E-05_r8, &
        0.3061E-05_r8, 0.4235E-05_r8, 0.6225E-05_r8, 0.8615E-05_r8, 0.1095E-04_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8,ip=24,24)/ &
        0.99980992_r8, 0.99961102_r8, 0.99922198_r8, 0.99852699_r8, 0.99732202_r8, &
        0.99510902_r8, 0.99140203_r8, 0.98550999_r8, 0.97672999_r8, 0.96399999_r8, &
        0.94602001_r8, 0.92101002_r8, 0.88709998_r8, 0.84310001_r8, 0.79020000_r8, &
        0.73189998_r8, 0.67299998_r8, 0.61619997_r8, 0.56060004_r8, 0.50400001_r8, &
        0.44610000_r8, 0.38880002_r8, 0.33530003_r8, 0.28740001_r8, 0.24390000_r8, &
        0.20179999_r8, 0.16009998_r8, 0.11979997_r8, 0.08260000_r8, 0.05049998_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8,ip=24,24)/ &
       -0.1529E-06_r8,-0.2005E-06_r8,-0.2861E-06_r8,-0.1652E-06_r8,-0.6334E-06_r8, &
       -0.1965E-05_r8,-0.5437E-05_r8,-0.1182E-04_r8,-0.2344E-04_r8,-0.4384E-04_r8, &
       -0.7982E-04_r8,-0.1398E-03_r8,-0.2321E-03_r8,-0.3616E-03_r8,-0.5274E-03_r8, &
       -0.7239E-03_r8,-0.9363E-03_r8,-0.1142E-02_r8,-0.1328E-02_r8,-0.1499E-02_r8, &
       -0.1645E-02_r8,-0.1748E-02_r8,-0.1804E-02_r8,-0.1834E-02_r8,-0.1867E-02_r8, &
       -0.1903E-02_r8,-0.1914E-02_r8,-0.1866E-02_r8,-0.1746E-02_r8,-0.1558E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30)_r8,ip=24,24)/ &
       -0.3638E-09_r8,-0.9313E-09_r8, 0.1703E-08_r8, 0.2081E-08_r8,-0.1251E-08_r8, &
       -0.1208E-08_r8,-0.6883E-08_r8,-0.1608E-07_r8,-0.4559E-07_r8,-0.1047E-06_r8, &
       -0.2040E-06_r8,-0.3312E-06_r8,-0.4624E-06_r8,-0.5198E-06_r8,-0.4326E-06_r8, &
       -0.1452E-06_r8, 0.3003E-06_r8, 0.7455E-06_r8, 0.1102E-05_r8, 0.1470E-05_r8, &
        0.1957E-05_r8, 0.2474E-05_r8, 0.2691E-05_r8, 0.2484E-05_r8, 0.2414E-05_r8, &
        0.3232E-05_r8, 0.5050E-05_r8, 0.7455E-05_r8, 0.9997E-05_r8, 0.1217E-04_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8,ip=25,25)/ &
        0.99980998_r8, 0.99960178_r8, 0.99920201_r8, 0.99852800_r8, 0.99729002_r8, &
        0.99498200_r8, 0.99102801_r8, 0.98461998_r8, 0.97465998_r8, 0.95982999_r8, &
        0.93866003_r8, 0.90968001_r8, 0.87140000_r8, 0.82340002_r8, 0.76770002_r8, &
        0.70860004_r8, 0.64999998_r8, 0.59290004_r8, 0.53610003_r8, 0.47860003_r8, &
        0.42110002_r8, 0.36610001_r8, 0.31639999_r8, 0.27200001_r8, 0.22960001_r8, &
        0.18690002_r8, 0.14429998_r8, 0.10380000_r8, 0.06739998_r8, 0.03740001_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8,ip=25,25)/ &
       -0.1453E-06_r8,-0.2529E-06_r8,-0.1807E-06_r8,-0.1109E-06_r8,-0.4469E-06_r8, &
       -0.1885E-05_r8,-0.4590E-05_r8,-0.1043E-04_r8,-0.2057E-04_r8,-0.3951E-04_r8, &
       -0.7466E-04_r8,-0.1356E-03_r8,-0.2341E-03_r8,-0.3783E-03_r8,-0.5688E-03_r8, &
       -0.7935E-03_r8,-0.1021E-02_r8,-0.1219E-02_r8,-0.1388E-02_r8,-0.1535E-02_r8, &
       -0.1653E-02_r8,-0.1726E-02_r8,-0.1768E-02_r8,-0.1813E-02_r8,-0.1874E-02_r8, &
       -0.1925E-02_r8,-0.1927E-02_r8,-0.1851E-02_r8,-0.1697E-02_r8,-0.1478E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30)_r8,ip=25,25)/ &
       -0.6257E-09_r8,-0.1382E-08_r8, 0.2095E-08_r8, 0.1863E-08_r8,-0.1834E-08_r8, &
       -0.2125E-08_r8,-0.6985E-08_r8,-0.1634E-07_r8,-0.4128E-07_r8,-0.9924E-07_r8, &
       -0.1938E-06_r8,-0.3275E-06_r8,-0.4556E-06_r8,-0.5046E-06_r8,-0.3633E-06_r8, &
        0.2484E-07_r8, 0.6195E-06_r8, 0.1249E-05_r8, 0.1731E-05_r8, 0.2053E-05_r8, &
        0.2358E-05_r8, 0.2569E-05_r8, 0.2342E-05_r8, 0.1883E-05_r8, 0.2103E-05_r8, &
        0.3570E-05_r8, 0.5973E-05_r8, 0.8752E-05_r8, 0.1140E-04_r8, 0.1328E-04_r8, &
      ! data ((c1(ip,iw),iw=1,30)_r8,ip=26,26)/ &
        0.99980712_r8, 0.99958581_r8, 0.99919039_r8, 0.99854302_r8, 0.99724799_r8, &
        0.99486500_r8, 0.99071401_r8, 0.98379999_r8, 0.97279000_r8, 0.95585001_r8, &
        0.93112999_r8, 0.89749998_r8, 0.85460001_r8, 0.80320001_r8, 0.74660003_r8, &
        0.68869996_r8, 0.63100004_r8, 0.57249999_r8, 0.51320004_r8, 0.45450002_r8, &
        0.39810002_r8, 0.34649998_r8, 0.30119997_r8, 0.25950003_r8, 0.21740001_r8, &
        0.17379999_r8, 0.13029999_r8, 0.08950001_r8, 0.05400002_r8, 0.02640003_r8, &
      ! data ((c2(ip,iw),iw=1,30)_r8,ip=26,26)/ &
       -0.1257E-06_r8,-0.2495E-06_r8,-0.1334E-06_r8,-0.8414E-07_r8,-0.1698E-06_r8, &
       -0.1346E-05_r8,-0.3692E-05_r8,-0.8625E-05_r8,-0.1750E-04_r8,-0.3483E-04_r8, &
       -0.6843E-04_r8,-0.1305E-03_r8,-0.2362E-03_r8,-0.3971E-03_r8,-0.6127E-03_r8, &
       -0.8621E-03_r8,-0.1101E-02_r8,-0.1297E-02_r8,-0.1452E-02_r8,-0.1570E-02_r8, &
       -0.1647E-02_r8,-0.1688E-02_r8,-0.1727E-02_r8,-0.1797E-02_r8,-0.1887E-02_r8, &
       -0.1947E-02_r8,-0.1935E-02_r8,-0.1833E-02_r8,-0.1647E-02_r8,-0.1401E-02_r8, &
      ! data ((c3(ip,iw),iw=1,30)_r8,ip=26,26)/ &
       -0.1222E-08_r8,-0.1164E-09_r8, 0.2285E-08_r8, 0.2037E-09_r8, 0.5675E-09_r8, &
       -0.5239E-08_r8,-0.9211E-08_r8,-0.1483E-07_r8,-0.3981E-07_r8,-0.9641E-07_r8, &
       -0.1717E-06_r8,-0.2796E-06_r8,-0.3800E-06_r8,-0.3762E-06_r8,-0.1936E-06_r8, &
        0.1920E-06_r8, 0.8335E-06_r8, 0.1691E-05_r8, 0.2415E-05_r8, 0.2767E-05_r8, &
        0.2823E-05_r8, 0.2551E-05_r8, 0.1839E-05_r8, 0.1314E-05_r8, 0.1960E-05_r8, &
        0.4003E-05_r8, 0.6909E-05_r8, 0.1004E-04_r8, 0.1273E-04_r8, 0.1423E-04_r8/),SHAPE=(/nx*3*nc/))


    it=0
    DO k=1,nx
       DO i=1,3
          DO j=1,nc
             it=it+1              
              IF(i==1)THEN
                ! WRITE(*,'(a5,2e16.9) ' )'c1',data1(it),c1(k,j)
                 c1(k,j)=data1(it)
              END IF 
              IF(i==2) THEN
                 !WRITE(*,'(a5,2e16.9) ' )'c2',data1(it),c2(k,j)
                 c2(k,j)=data1(it)
              END IF  
              IF(i==3) THEN
                 !WRITE(*,'(a5,2e16.9) ' )'c3',data1(it),c3(k,j)
                 c3(k,j) = data1(it)
              END IF
          END DO
       END DO
    END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      data2(1:nx*3*no)=RESHAPE(SOURCE=(/&
       ! data ((o1(ip,iw),iw=1,21)_r8, ip= 1_r8, 1)/ & 
        0.99999344_r8,  0.99998689_r8,  0.99997336_r8,  0.99994606_r8,  0.99989170_r8, & 
        0.99978632_r8,  0.99957907_r8,  0.99918377_r8,  0.99844402_r8,  0.99712098_r8, & 
        0.99489498_r8,  0.99144602_r8,  0.98655999_r8,  0.98008001_r8,  0.97165000_r8, & 
        0.96043998_r8,  0.94527000_r8,  0.92462999_r8,  0.89709997_r8,  0.86180001_r8, & 
        0.81800002_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip= 1_r8, 1)/  &
        0.6531E-10_r8,  0.5926E-10_r8, -0.1646E-09_r8, -0.1454E-08_r8, -0.7376E-08_r8, &
       -0.2968E-07_r8, -0.1071E-06_r8, -0.3584E-06_r8, -0.1125E-05_r8, -0.3289E-05_r8, &
       -0.8760E-05_r8, -0.2070E-04_r8, -0.4259E-04_r8, -0.7691E-04_r8, -0.1264E-03_r8, &
       -0.1957E-03_r8, -0.2895E-03_r8, -0.4107E-03_r8, -0.5588E-03_r8, -0.7300E-03_r8, &
       -0.9199E-03_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip= 1_r8, 1)/ &
       -0.2438E-10_r8, -0.4826E-10_r8, -0.9474E-10_r8, -0.1828E-09_r8, -0.3406E-09_r8, &
       -0.6223E-09_r8, -0.1008E-08_r8, -0.1412E-08_r8, -0.1244E-08_r8,  0.8485E-09_r8, &
        0.6343E-08_r8,  0.1201E-07_r8,  0.2838E-08_r8, -0.4024E-07_r8, -0.1257E-06_r8, &
       -0.2566E-06_r8, -0.4298E-06_r8, -0.6184E-06_r8, -0.7657E-06_r8, -0.8153E-06_r8, &
       -0.7552E-06_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip= 2_r8, 2)/ &
        0.99999344_r8,  0.99998689_r8,  0.99997348_r8,  0.99994606_r8,  0.99989170_r8, &
        0.99978632_r8,  0.99957907_r8,  0.99918377_r8,  0.99844402_r8,  0.99712098_r8, &
        0.99489498_r8,  0.99144298_r8,  0.98654997_r8,  0.98006999_r8,  0.97162998_r8, &
        0.96042001_r8,  0.94520003_r8,  0.92449999_r8,  0.89690000_r8,  0.86140001_r8, &
        0.81739998_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip= 2_r8, 2)/ &
        0.6193E-10_r8,  0.5262E-10_r8, -0.1774E-09_r8, -0.1478E-08_r8, -0.7416E-08_r8, &
       -0.2985E-07_r8, -0.1071E-06_r8, -0.3584E-06_r8, -0.1124E-05_r8, -0.3287E-05_r8, &
       -0.8753E-05_r8, -0.2069E-04_r8, -0.4256E-04_r8, -0.7686E-04_r8, -0.1264E-03_r8, &
       -0.1956E-03_r8, -0.2893E-03_r8, -0.4103E-03_r8, -0.5580E-03_r8, -0.7285E-03_r8, &
       -0.9171E-03_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip= 2_r8, 2)/ &
       -0.2436E-10_r8, -0.4822E-10_r8, -0.9466E-10_r8, -0.1827E-09_r8, -0.3404E-09_r8, &
       -0.6220E-09_r8, -0.1008E-08_r8, -0.1414E-08_r8, -0.1247E-08_r8,  0.8360E-09_r8, &
        0.6312E-08_r8,  0.1194E-07_r8,  0.2753E-08_r8, -0.4040E-07_r8, -0.1260E-06_r8, &
       -0.2571E-06_r8, -0.4307E-06_r8, -0.6202E-06_r8, -0.7687E-06_r8, -0.8204E-06_r8, &
       -0.7636E-06_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip= 3_r8, 3)/ &
        0.99999344_r8,  0.99998689_r8,  0.99997348_r8,  0.99994606_r8,  0.99989170_r8, &
        0.99978632_r8,  0.99957907_r8,  0.99918377_r8,  0.99844402_r8,  0.99712098_r8, &
        0.99489301_r8,  0.99143898_r8,  0.98654997_r8,  0.98005998_r8,  0.97158998_r8, &
        0.96035999_r8,  0.94509000_r8,  0.92431998_r8,  0.89660001_r8,  0.86080003_r8, &
        0.81639999_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip= 3_r8, 3)/ &
        0.5658E-10_r8,  0.4212E-10_r8, -0.1977E-09_r8, -0.1516E-08_r8, -0.7481E-08_r8, &
       -0.2995E-07_r8, -0.1072E-06_r8, -0.3583E-06_r8, -0.1123E-05_r8, -0.3283E-05_r8, &
       -0.8744E-05_r8, -0.2067E-04_r8, -0.4252E-04_r8, -0.7679E-04_r8, -0.1262E-03_r8, &
       -0.1953E-03_r8, -0.2889E-03_r8, -0.4096E-03_r8, -0.5567E-03_r8, -0.7263E-03_r8, &
       -0.9130E-03_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip= 3_r8, 3)/ &
       -0.2433E-10_r8, -0.4815E-10_r8, -0.9453E-10_r8, -0.1825E-09_r8, -0.3400E-09_r8, &
       -0.6215E-09_r8, -0.1007E-08_r8, -0.1415E-08_r8, -0.1253E-08_r8,  0.8143E-09_r8, &
        0.6269E-08_r8,  0.1186E-07_r8,  0.2604E-08_r8, -0.4067E-07_r8, -0.1264E-06_r8, &
       -0.2579E-06_r8, -0.4321E-06_r8, -0.6229E-06_r8, -0.7732E-06_r8, -0.8277E-06_r8, &
       -0.7752E-06_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip= 4_r8, 4)/ &
        0.99999344_r8,  0.99998689_r8,  0.99997348_r8,  0.99994606_r8,  0.99989200_r8, &
        0.99978632_r8,  0.99957907_r8,  0.99918377_r8,  0.99844402_r8,  0.99711901_r8, &
        0.99489301_r8,  0.99143499_r8,  0.98653001_r8,  0.98003000_r8,  0.97153997_r8, &
        0.96026999_r8,  0.94493997_r8,  0.92404002_r8,  0.89609998_r8,  0.85990000_r8, &
        0.81480002_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip= 4_r8, 4)/ &
        0.4814E-10_r8,  0.2552E-10_r8, -0.2298E-09_r8, -0.1576E-08_r8, -0.7579E-08_r8, &
       -0.3009E-07_r8, -0.1074E-06_r8, -0.3581E-06_r8, -0.1122E-05_r8, -0.3278E-05_r8, &
       -0.8729E-05_r8, -0.2063E-04_r8, -0.4245E-04_r8, -0.7667E-04_r8, -0.1260E-03_r8, &
       -0.1950E-03_r8, -0.2883E-03_r8, -0.4086E-03_r8, -0.5549E-03_r8, -0.7229E-03_r8, &
       -0.9071E-03_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip= 4_r8, 4)/ &
       -0.2428E-10_r8, -0.4805E-10_r8, -0.9433E-10_r8, -0.1821E-09_r8, -0.3394E-09_r8, &
       -0.6206E-09_r8, -0.1008E-08_r8, -0.1416E-08_r8, -0.1261E-08_r8,  0.7860E-09_r8, &
        0.6188E-08_r8,  0.1171E-07_r8,  0.2389E-08_r8, -0.4109E-07_r8, -0.1271E-06_r8, &
       -0.2591E-06_r8, -0.4344E-06_r8, -0.6267E-06_r8, -0.7797E-06_r8, -0.8378E-06_r8, &
       -0.7901E-06_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip= 5_r8, 5)/ &
        0.99999344_r8,  0.99998689_r8,  0.99997348_r8,  0.99994606_r8,  0.99989200_r8, &
        0.99978638_r8,  0.99957907_r8,  0.99918377_r8,  0.99844402_r8,  0.99711901_r8, &
        0.99488801_r8,  0.99142599_r8,  0.98650998_r8,  0.97999001_r8,  0.97148001_r8, &
        0.96011001_r8,  0.94467002_r8,  0.92356998_r8,  0.89530003_r8,  0.85860002_r8, &
        0.81250000_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip= 5_r8, 5)/ &
        0.3482E-10_r8, -0.6492E-12_r8, -0.2805E-09_r8, -0.1671E-08_r8, -0.7740E-08_r8, &
       -0.3032E-07_r8, -0.1076E-06_r8, -0.3582E-06_r8, -0.1120E-05_r8, -0.3270E-05_r8, &
       -0.8704E-05_r8, -0.2058E-04_r8, -0.4235E-04_r8, -0.7649E-04_r8, -0.1257E-03_r8, &
       -0.1945E-03_r8, -0.2874E-03_r8, -0.4070E-03_r8, -0.5521E-03_r8, -0.7181E-03_r8, &
       -0.8990E-03_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip= 5_r8, 5)/ &
       -0.2419E-10_r8, -0.4788E-10_r8, -0.9401E-10_r8, -0.1815E-09_r8, -0.3385E-09_r8, &
       -0.6192E-09_r8, -0.1006E-08_r8, -0.1417E-08_r8, -0.1273E-08_r8,  0.7404E-09_r8, &
        0.6068E-08_r8,  0.1148E-07_r8,  0.2021E-08_r8, -0.4165E-07_r8, -0.1281E-06_r8, &
       -0.2609E-06_r8, -0.4375E-06_r8, -0.6323E-06_r8, -0.7887E-06_r8, -0.8508E-06_r8, &
       -0.8067E-06_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip= 6_r8, 6)/ &
        0.99999344_r8,  0.99998689_r8,  0.99997348_r8,  0.99994606_r8,  0.99989200_r8, &
        0.99978638_r8,  0.99957931_r8,  0.99918377_r8,  0.99844301_r8,  0.99711698_r8, &
        0.99488401_r8,  0.99141300_r8,  0.98648000_r8,  0.97992003_r8,  0.97135001_r8, &
        0.95989001_r8,  0.94428003_r8,  0.92286998_r8,  0.89410001_r8,  0.85640001_r8, &
        0.80890000_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip= 6_r8, 6)/ &
        0.1388E-10_r8, -0.4180E-10_r8, -0.3601E-09_r8, -0.1820E-08_r8, -0.7993E-08_r8, &
       -0.3068E-07_r8, -0.1081E-06_r8, -0.3580E-06_r8, -0.1117E-05_r8, -0.3257E-05_r8, &
       -0.8667E-05_r8, -0.2049E-04_r8, -0.4218E-04_r8, -0.7620E-04_r8, -0.1253E-03_r8, &
       -0.1937E-03_r8, -0.2860E-03_r8, -0.4047E-03_r8, -0.5481E-03_r8, -0.7115E-03_r8, &
       -0.8885E-03_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip= 6_r8, 6)/ &
       -0.2406E-10_r8, -0.4762E-10_r8, -0.9351E-10_r8, -0.1806E-09_r8, -0.3370E-09_r8, &
       -0.6170E-09_r8, -0.1004E-08_r8, -0.1417E-08_r8, -0.1297E-08_r8,  0.6738E-09_r8, &
        0.5895E-08_r8,  0.1113E-07_r8,  0.1466E-08_r8, -0.4265E-07_r8, -0.1298E-06_r8, &
       -0.2636E-06_r8, -0.4423E-06_r8, -0.6402E-06_r8, -0.8005E-06_r8, -0.8658E-06_r8, &
       -0.8222E-06_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip= 7_r8, 7)/ &
        0.99999344_r8,  0.99998689_r8,  0.99997348_r8,  0.99994630_r8,  0.99989200_r8, &
        0.99978638_r8,  0.99957931_r8,  0.99918360_r8,  0.99844301_r8,  0.99711502_r8, &
        0.99487501_r8,  0.99138802_r8,  0.98642999_r8,  0.97982001_r8,  0.97114998_r8, &
        0.95954001_r8,  0.94363999_r8,  0.92176998_r8,  0.89219999_r8,  0.85329998_r8, &
        0.80379999_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip= 7_r8, 7)/ &
       -0.1889E-10_r8, -0.1062E-09_r8, -0.4847E-09_r8, -0.2053E-08_r8, -0.8389E-08_r8, &
       -0.3140E-07_r8, -0.1089E-06_r8, -0.3577E-06_r8, -0.1112E-05_r8, -0.3236E-05_r8, &
       -0.8607E-05_r8, -0.2035E-04_r8, -0.4192E-04_r8, -0.7576E-04_r8, -0.1245E-03_r8, &
       -0.1925E-03_r8, -0.2840E-03_r8, -0.4013E-03_r8, -0.5427E-03_r8, -0.7029E-03_r8, &
       -0.8756E-03_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip= 7_r8, 7)/ &
       -0.2385E-10_r8, -0.4722E-10_r8, -0.9273E-10_r8, -0.1791E-09_r8, -0.3348E-09_r8, &
       -0.6121E-09_r8, -0.9974E-09_r8, -0.1422E-08_r8, -0.1326E-08_r8,  0.5603E-09_r8, &
        0.5604E-08_r8,  0.1061E-07_r8,  0.6106E-09_r8, -0.4398E-07_r8, -0.1321E-06_r8, &
       -0.2676E-06_r8, -0.4490E-06_r8, -0.6507E-06_r8, -0.8145E-06_r8, -0.8801E-06_r8, &
       -0.8311E-06_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip= 8_r8, 8)/ &
        0.99999344_r8,  0.99998689_r8,  0.99997348_r8,  0.99994630_r8,  0.99989229_r8, &
        0.99978650_r8,  0.99957931_r8,  0.99918288_r8,  0.99844098_r8,  0.99711001_r8, &
        0.99486202_r8,  0.99135500_r8,  0.98635000_r8,  0.97965997_r8,  0.97083998_r8, &
        0.95898998_r8,  0.94266999_r8,  0.92009997_r8,  0.88929999_r8,  0.84860003_r8, &
        0.79640001_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip= 8_r8, 8)/ &
       -0.6983E-10_r8, -0.2063E-09_r8, -0.6785E-09_r8, -0.2416E-08_r8, -0.9000E-08_r8, &
       -0.3243E-07_r8, -0.1100E-06_r8, -0.3574E-06_r8, -0.1104E-05_r8, -0.3205E-05_r8, &
       -0.8516E-05_r8, -0.2014E-04_r8, -0.4151E-04_r8, -0.7508E-04_r8, -0.1234E-03_r8, &
       -0.1907E-03_r8, -0.2811E-03_r8, -0.3966E-03_r8, -0.5355E-03_r8, -0.6924E-03_r8, &
       -0.8613E-03_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip= 8_r8, 8)/ &
       -0.2353E-10_r8, -0.4659E-10_r8, -0.9153E-10_r8, -0.1769E-09_r8, -0.3313E-09_r8, &
       -0.6054E-09_r8, -0.9899E-09_r8, -0.1430E-08_r8, -0.1375E-08_r8,  0.3874E-09_r8, &
        0.5171E-08_r8,  0.9807E-08_r8, -0.7345E-09_r8, -0.4604E-07_r8, -0.1356E-06_r8, &
       -0.2731E-06_r8, -0.4577E-06_r8, -0.6632E-06_r8, -0.8284E-06_r8, -0.8894E-06_r8, &
       -0.8267E-06_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip= 9_r8, 9)/ &
        0.99999344_r8,  0.99998689_r8,  0.99997360_r8,  0.99994630_r8,  0.99989229_r8, &
        0.99978650_r8,  0.99957961_r8,  0.99918252_r8,  0.99843901_r8,  0.99710202_r8, &
        0.99484003_r8,  0.99130303_r8,  0.98623002_r8,  0.97940999_r8,  0.97038001_r8, &
        0.95815003_r8,  0.94119000_r8,  0.91755998_r8,  0.88510001_r8,  0.84189999_r8, &
        0.78610003_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip= 9_r8, 9)/ &
       -0.1481E-09_r8, -0.3601E-09_r8, -0.9762E-09_r8, -0.2973E-08_r8, -0.1014E-07_r8, &
       -0.3421E-07_r8, -0.1121E-06_r8, -0.3569E-06_r8, -0.1092E-05_r8, -0.3156E-05_r8, &
       -0.8375E-05_r8, -0.1981E-04_r8, -0.4090E-04_r8, -0.7405E-04_r8, -0.1218E-03_r8, &
       -0.1881E-03_r8, -0.2770E-03_r8, -0.3906E-03_r8, -0.5269E-03_r8, -0.6810E-03_r8, &
       -0.8471E-03_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip= 9_r8, 9)/ &
       -0.2304E-10_r8, -0.4564E-10_r8, -0.8969E-10_r8, -0.1735E-09_r8, -0.3224E-09_r8, &
       -0.5933E-09_r8, -0.9756E-09_r8, -0.1428E-08_r8, -0.1446E-08_r8,  0.1156E-09_r8, &
        0.4499E-08_r8,  0.8469E-08_r8, -0.2720E-08_r8, -0.4904E-07_r8, -0.1401E-06_r8, &
       -0.2801E-06_r8, -0.4681E-06_r8, -0.6761E-06_r8, -0.8387E-06_r8, -0.8879E-06_r8, &
       -0.8040E-06_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip=10,10)/ &
        0.99999344_r8,  0.99998689_r8,  0.99997360_r8,  0.99994630_r8,  0.99989259_r8, &
        0.99978650_r8,  0.99957931_r8,  0.99918163_r8,  0.99843597_r8,  0.99709100_r8, &
        0.99480897_r8,  0.99122101_r8,  0.98604000_r8,  0.97902000_r8,  0.96965003_r8, &
        0.95684999_r8,  0.93896997_r8,  0.91386002_r8,  0.87910002_r8,  0.83249998_r8, &
        0.77200001_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip=10,10)/ &
       -0.2661E-09_r8, -0.5923E-09_r8, -0.1426E-08_r8, -0.3816E-08_r8, -0.1159E-07_r8, &
       -0.3654E-07_r8, -0.1143E-06_r8, -0.3559E-06_r8, -0.1074E-05_r8, -0.3083E-05_r8, &
       -0.8159E-05_r8, -0.1932E-04_r8, -0.3998E-04_r8, -0.7253E-04_r8, -0.1194E-03_r8, &
       -0.1845E-03_r8, -0.2718E-03_r8, -0.3833E-03_r8, -0.5176E-03_r8, -0.6701E-03_r8, &
       -0.8354E-03_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip=10,10)/ &
       -0.2232E-10_r8, -0.4421E-10_r8, -0.8695E-10_r8, -0.1684E-09_r8, -0.3141E-09_r8, &
       -0.5765E-09_r8, -0.9606E-09_r8, -0.1434E-08_r8, -0.1551E-08_r8, -0.2663E-09_r8, &
        0.3515E-08_r8,  0.6549E-08_r8, -0.5479E-08_r8, -0.5312E-07_r8, -0.1460E-06_r8, &
       -0.2883E-06_r8, -0.4787E-06_r8, -0.6863E-06_r8, -0.8399E-06_r8, -0.8703E-06_r8, &
       -0.7602E-06_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip=11,11)/ &
        0.99999356_r8,  0.99998701_r8,  0.99997360_r8,  0.99994630_r8,  0.99989289_r8, &
        0.99978679_r8,  0.99957907_r8,  0.99917960_r8,  0.99843001_r8,  0.99707502_r8, &
        0.99475998_r8,  0.99109501_r8,  0.98575002_r8,  0.97843999_r8,  0.96855003_r8, &
        0.95494002_r8,  0.93572998_r8,  0.90853000_r8,  0.87070000_r8,  0.81970000_r8, &
        0.75380003_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip=11,11)/ &
       -0.4394E-09_r8, -0.9330E-09_r8, -0.2086E-08_r8, -0.5054E-08_r8, -0.1373E-07_r8, &
       -0.3971E-07_r8, -0.1178E-06_r8, -0.3546E-06_r8, -0.1049E-05_r8, -0.2976E-05_r8, &
       -0.7847E-05_r8, -0.1860E-04_r8, -0.3864E-04_r8, -0.7038E-04_r8, -0.1162E-03_r8, &
       -0.1798E-03_r8, -0.2654E-03_r8, -0.3754E-03_r8, -0.5091E-03_r8, -0.6621E-03_r8, &
       -0.8286E-03_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip=11,11)/ &
       -0.2127E-10_r8, -0.4216E-10_r8, -0.8300E-10_r8, -0.1611E-09_r8, -0.3019E-09_r8, &
       -0.5597E-09_r8, -0.9431E-09_r8, -0.1450E-08_r8, -0.1694E-08_r8, -0.7913E-09_r8, &
        0.2144E-08_r8,  0.3990E-08_r8, -0.9282E-08_r8, -0.5810E-07_r8, -0.1525E-06_r8, &
       -0.2965E-06_r8, -0.4869E-06_r8, -0.6894E-06_r8, -0.8281E-06_r8, -0.8350E-06_r8, &
       -0.6956E-06_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip=12,12)/  &
        0.99999368_r8,  0.99998701_r8,  0.99997377_r8,  0.99994630_r8,  0.99989259_r8,  &
        0.99978709_r8,  0.99957848_r8,  0.99917740_r8,  0.99842203_r8,  0.99704897_r8,  &
        0.99468797_r8,  0.99090999_r8,  0.98532999_r8,  0.97758001_r8,  0.96693999_r8,  &
        0.95213997_r8,  0.93109000_r8,  0.90110999_r8,  0.85930002_r8,  0.80290002_r8,  &
        0.73019999_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip=12,12)/  &
       -0.6829E-09_r8, -0.1412E-08_r8, -0.3014E-08_r8, -0.6799E-08_r8, -0.1675E-07_r8,  &
       -0.4450E-07_r8, -0.1235E-06_r8, -0.3538E-06_r8, -0.1014E-05_r8, -0.2827E-05_r8,  &
       -0.7407E-05_r8, -0.1759E-04_r8, -0.3676E-04_r8, -0.6744E-04_r8, -0.1120E-03_r8,  &
       -0.1742E-03_r8, -0.2585E-03_r8, -0.3683E-03_r8, -0.5034E-03_r8, -0.6594E-03_r8,  &
       -0.8290E-03_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip=12,12)/  &
       -0.1985E-10_r8, -0.3937E-10_r8, -0.7761E-10_r8, -0.1511E-09_r8, -0.2855E-09_r8,  &
       -0.5313E-09_r8, -0.9251E-09_r8, -0.1470E-08_r8, -0.1898E-08_r8, -0.1519E-08_r8,  &
        0.2914E-09_r8,  0.5675E-09_r8, -0.1405E-07_r8, -0.6359E-07_r8, -0.1584E-06_r8,  &
       -0.3020E-06_r8, -0.4893E-06_r8, -0.6821E-06_r8, -0.8021E-06_r8, -0.7834E-06_r8,  &
       -0.6105E-06_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip=13,13)/  &
        0.99999368_r8,  0.99998701_r8,  0.99997389_r8,  0.99994695_r8,  0.99989289_r8,  &
        0.99978721_r8,  0.99957782_r8,  0.99917412_r8,  0.99840999_r8,  0.99701297_r8,  &
        0.99458599_r8,  0.99064600_r8,  0.98471999_r8,  0.97632003_r8,  0.96464998_r8,  &
        0.94819999_r8,  0.92467999_r8,  0.89109999_r8,  0.84430003_r8,  0.78139997_r8,  &
        0.70070004_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip=13,13)/  &
       -0.1004E-08_r8, -0.2043E-08_r8, -0.4239E-08_r8, -0.9104E-08_r8, -0.2075E-07_r8,  &
       -0.5096E-07_r8, -0.1307E-06_r8, -0.3520E-06_r8, -0.9671E-06_r8, -0.2630E-05_r8,  &
       -0.6825E-05_r8, -0.1624E-04_r8, -0.3429E-04_r8, -0.6369E-04_r8, -0.1069E-03_r8,  &
       -0.1680E-03_r8, -0.2520E-03_r8, -0.3635E-03_r8, -0.5029E-03_r8, -0.6647E-03_r8,  &
       -0.8390E-03_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip=13,13)/  &
       -0.1807E-10_r8, -0.3587E-10_r8, -0.7085E-10_r8, -0.1385E-09_r8, -0.2648E-09_r8,  &
       -0.4958E-09_r8, -0.8900E-09_r8, -0.1473E-08_r8, -0.2112E-08_r8, -0.2399E-08_r8,  &
       -0.2002E-08_r8, -0.3646E-08_r8, -0.1931E-07_r8, -0.6852E-07_r8, -0.1618E-06_r8,  &
       -0.3021E-06_r8, -0.4828E-06_r8, -0.6634E-06_r8, -0.7643E-06_r8, -0.7177E-06_r8,  &
       -0.5054E-06_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip=14,14)/ &
        0.99999368_r8,  0.99998713_r8,  0.99997389_r8,  0.99994725_r8,  0.99989289_r8, &
        0.99978679_r8,  0.99957597_r8,  0.99916971_r8,  0.99839503_r8,  0.99696702_r8, &
        0.99444997_r8,  0.99028301_r8,  0.98387003_r8,  0.97457999_r8,  0.96148002_r8, &
        0.94284999_r8,  0.91613001_r8,  0.87809998_r8,  0.82520002_r8,  0.75489998_r8, &
        0.66520000_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip=14,14)/ &
       -0.1387E-08_r8, -0.2798E-08_r8, -0.5706E-08_r8, -0.1187E-07_r8, -0.2564E-07_r8, &
       -0.5866E-07_r8, -0.1398E-06_r8, -0.3516E-06_r8, -0.9148E-06_r8, -0.2398E-05_r8, &
       -0.6122E-05_r8, -0.1459E-04_r8, -0.3125E-04_r8, -0.5923E-04_r8, -0.1013E-03_r8, &
       -0.1620E-03_r8, -0.2473E-03_r8, -0.3631E-03_r8, -0.5098E-03_r8, -0.6800E-03_r8, &
       -0.8603E-03_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip=14,14)/ &
       -0.1610E-10_r8, -0.3200E-10_r8, -0.6337E-10_r8, -0.1245E-09_r8, -0.2408E-09_r8, &
       -0.4533E-09_r8, -0.8405E-09_r8, -0.1464E-08_r8, -0.2337E-08_r8, -0.3341E-08_r8, &
       -0.4467E-08_r8, -0.8154E-08_r8, -0.2436E-07_r8, -0.7128E-07_r8, -0.1604E-06_r8, &
       -0.2945E-06_r8, -0.4666E-06_r8, -0.6357E-06_r8, -0.7187E-06_r8, -0.6419E-06_r8, &
       -0.3795E-06_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip=15,15)/ &
        0.99999410_r8,  0.99998724_r8,  0.99997455_r8,  0.99994725_r8,  0.99989331_r8, &
        0.99978632_r8,  0.99957472_r8,  0.99916393_r8,  0.99837703_r8,  0.99690801_r8, &
        0.99427801_r8,  0.98982000_r8,  0.98277998_r8,  0.97232002_r8,  0.95731997_r8, &
        0.93585998_r8,  0.90521002_r8,  0.86180001_r8,  0.80190003_r8,  0.72290003_r8, &
        0.62380004_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip=15,15)/ &
       -0.1788E-08_r8, -0.3588E-08_r8, -0.7244E-08_r8, -0.1479E-07_r8, -0.3083E-07_r8, &
       -0.6671E-07_r8, -0.1497E-06_r8, -0.3519E-06_r8, -0.8607E-06_r8, -0.2154E-05_r8, &
       -0.5364E-05_r8, -0.1276E-04_r8, -0.2785E-04_r8, -0.5435E-04_r8, -0.9573E-04_r8, &
       -0.1570E-03_r8, -0.2455E-03_r8, -0.3682E-03_r8, -0.5253E-03_r8, -0.7065E-03_r8, &
       -0.8938E-03_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip=15,15)/ &
       -0.1429E-10_r8, -0.2843E-10_r8, -0.5645E-10_r8, -0.1115E-09_r8, -0.2181E-09_r8, &
       -0.4200E-09_r8, -0.7916E-09_r8, -0.1460E-08_r8, -0.2542E-08_r8, -0.4168E-08_r8, &
       -0.6703E-08_r8, -0.1215E-07_r8, -0.2821E-07_r8, -0.7073E-07_r8, -0.1530E-06_r8, &
       -0.2791E-06_r8, -0.4426E-06_r8, -0.6027E-06_r8, -0.6707E-06_r8, -0.5591E-06_r8, &
       -0.2328E-06_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip=16,16)/ &
        0.99999434_r8,  0.99998778_r8,  0.99997467_r8,  0.99994761_r8,  0.99989331_r8, &
        0.99978602_r8,  0.99957269_r8,  0.99915779_r8,  0.99835497_r8,  0.99684399_r8, &
        0.99408400_r8,  0.98929000_r8,  0.98148000_r8,  0.96954000_r8,  0.95212001_r8, &
        0.92719001_r8,  0.89170003_r8,  0.84200001_r8,  0.77420002_r8,  0.68620002_r8, &
        0.57780004_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip=16,16)/ &
       -0.2141E-08_r8, -0.4286E-08_r8, -0.8603E-08_r8, -0.1737E-07_r8, -0.3548E-07_r8, &
       -0.7410E-07_r8, -0.1590E-06_r8, -0.3537E-06_r8, -0.8142E-06_r8, -0.1935E-05_r8, &
       -0.4658E-05_r8, -0.1099E-04_r8, -0.2444E-04_r8, -0.4948E-04_r8, -0.9067E-04_r8, &
       -0.1538E-03_r8, -0.2474E-03_r8, -0.3793E-03_r8, -0.5495E-03_r8, -0.7439E-03_r8, &
       -0.9383E-03_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip=16,16)/ &
       -0.1295E-10_r8, -0.2581E-10_r8, -0.5136E-10_r8, -0.1019E-09_r8, -0.2011E-09_r8, &
       -0.3916E-09_r8, -0.7585E-09_r8, -0.1439E-08_r8, -0.2648E-08_r8, -0.4747E-08_r8, &
       -0.8301E-08_r8, -0.1499E-07_r8, -0.3024E-07_r8, -0.6702E-07_r8, -0.1399E-06_r8, &
       -0.2564E-06_r8, -0.4117E-06_r8, -0.5669E-06_r8, -0.6239E-06_r8, -0.4748E-06_r8, &
       -0.7013E-07_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip=17,17)/ &
        0.99999434_r8,  0.99998778_r8,  0.99997479_r8,  0.99994791_r8,  0.99989331_r8, &
        0.99978608_r8,  0.99957120_r8,  0.99915212_r8,  0.99833500_r8,  0.99677801_r8, &
        0.99388403_r8,  0.98873001_r8,  0.98005998_r8,  0.96639001_r8,  0.94606000_r8, &
        0.91689998_r8,  0.87580001_r8,  0.81889999_r8,  0.74280000_r8,  0.64559996_r8, &
        0.52869999_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip=17,17)/ &
       -0.2400E-08_r8, -0.4796E-08_r8, -0.9599E-08_r8, -0.1927E-07_r8, -0.3892E-07_r8, &
       -0.7954E-07_r8, -0.1661E-06_r8, -0.3540E-06_r8, -0.7780E-06_r8, -0.1763E-05_r8, &
       -0.4092E-05_r8, -0.9512E-05_r8, -0.2142E-04_r8, -0.4502E-04_r8, -0.8640E-04_r8, &
       -0.1525E-03_r8, -0.2526E-03_r8, -0.3955E-03_r8, -0.5805E-03_r8, -0.7897E-03_r8, &
       -0.9899E-03_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip=17,17)/ &
       -0.1220E-10_r8, -0.2432E-10_r8, -0.4845E-10_r8, -0.9640E-10_r8, -0.1912E-09_r8, &
       -0.3771E-09_r8, -0.7392E-09_r8, -0.1420E-08_r8, -0.2702E-08_r8, -0.5049E-08_r8, &
       -0.9214E-08_r8, -0.1659E-07_r8, -0.3101E-07_r8, -0.6162E-07_r8, -0.1235E-06_r8, &
       -0.2287E-06_r8, -0.3755E-06_r8, -0.5274E-06_r8, -0.5790E-06_r8, -0.3947E-06_r8, &
        0.1003E-06_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip=18,18)/  &
        0.99999464_r8,  0.99998808_r8,  0.99997497_r8,  0.99994791_r8,  0.99989331_r8,  &
        0.99978518_r8,  0.99957031_r8,  0.99914658_r8,  0.99831802_r8,  0.99671799_r8,  &
        0.99370098_r8,  0.98821002_r8,  0.97867000_r8,  0.96313000_r8,  0.93948001_r8,  &
        0.90534002_r8,  0.85769999_r8,  0.79310000_r8,  0.70840001_r8,  0.60290003_r8,  &
        0.47930002_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip=18,18)/  &
       -0.2557E-08_r8, -0.5106E-08_r8, -0.1020E-07_r8, -0.2043E-07_r8, -0.4103E-07_r8,  &
       -0.8293E-07_r8, -0.1697E-06_r8, -0.3531E-06_r8, -0.7531E-06_r8, -0.1645E-05_r8,  &
       -0.3690E-05_r8, -0.8411E-05_r8, -0.1902E-04_r8, -0.4118E-04_r8, -0.8276E-04_r8,  &
       -0.1525E-03_r8, -0.2601E-03_r8, -0.4147E-03_r8, -0.6149E-03_r8, -0.8384E-03_r8,  &
       -0.1042E-02_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip=18,18)/  &
       -0.1189E-10_r8, -0.2372E-10_r8, -0.4729E-10_r8, -0.9421E-10_r8, -0.1873E-09_r8,  &
       -0.3713E-09_r8, -0.7317E-09_r8, -0.1437E-08_r8, -0.2764E-08_r8, -0.5243E-08_r8,  &
       -0.9691E-08_r8, -0.1751E-07_r8, -0.3122E-07_r8, -0.5693E-07_r8, -0.1076E-06_r8,  &
       -0.1981E-06_r8, -0.3324E-06_r8, -0.4785E-06_r8, -0.5280E-06_r8, -0.3174E-06_r8,  &
        0.2672E-06_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip=19,19)/  &
        0.99999464_r8,  0.99998820_r8,  0.99997509_r8,  0.99994779_r8,  0.99989331_r8,  &
        0.99978518_r8,  0.99956989_r8,  0.99914283_r8,  0.99830401_r8,  0.99667197_r8,  &
        0.99355298_r8,  0.98776001_r8,  0.97741997_r8,  0.96007001_r8,  0.93285000_r8,  &
        0.89310002_r8,  0.83819997_r8,  0.76520002_r8,  0.67250001_r8,  0.56000000_r8,  &
        0.43199998_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip=19,19)/  &
       -0.2630E-08_r8, -0.5249E-08_r8, -0.1048E-07_r8, -0.2096E-07_r8, -0.4198E-07_r8,  &
       -0.8440E-07_r8, -0.1710E-06_r8, -0.3513E-06_r8, -0.7326E-06_r8, -0.1562E-05_r8,  &
       -0.3416E-05_r8, -0.7637E-05_r8, -0.1719E-04_r8, -0.3795E-04_r8, -0.7926E-04_r8,  &
       -0.1524E-03_r8, -0.2680E-03_r8, -0.4344E-03_r8, -0.6486E-03_r8, -0.8838E-03_r8,  &
       -0.1089E-02_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip=19,19)/  &
       -0.1188E-10_r8, -0.2369E-10_r8, -0.4725E-10_r8, -0.9417E-10_r8, -0.1875E-09_r8,  &
       -0.3725E-09_r8, -0.7365E-09_r8, -0.1445E-08_r8, -0.2814E-08_r8, -0.5384E-08_r8,  &
       -0.1008E-07_r8, -0.1816E-07_r8, -0.3179E-07_r8, -0.5453E-07_r8, -0.9500E-07_r8,  &
       -0.1679E-06_r8, -0.2819E-06_r8, -0.4109E-06_r8, -0.4555E-06_r8, -0.2283E-06_r8,  &
        0.4283E-06_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip=20,20)/ &
        0.99999487_r8,  0.99998832_r8,  0.99997520_r8,  0.99994791_r8,  0.99989331_r8, &
        0.99978459_r8,  0.99956900_r8,  0.99913990_r8,  0.99829400_r8,  0.99663699_r8, &
        0.99344099_r8,  0.98741001_r8,  0.97643000_r8,  0.95743001_r8,  0.92672002_r8, &
        0.88099998_r8,  0.81809998_r8,  0.73660004_r8,  0.63620001_r8,  0.51880002_r8, &
        0.38880002_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip=20,20)/ &
       -0.2651E-08_r8, -0.5291E-08_r8, -0.1056E-07_r8, -0.2110E-07_r8, -0.4221E-07_r8, &
       -0.8462E-07_r8, -0.1705E-06_r8, -0.3466E-06_r8, -0.7155E-06_r8, -0.1501E-05_r8, &
       -0.3223E-05_r8, -0.7079E-05_r8, -0.1581E-04_r8, -0.3517E-04_r8, -0.7553E-04_r8, &
       -0.1510E-03_r8, -0.2746E-03_r8, -0.4528E-03_r8, -0.6789E-03_r8, -0.9214E-03_r8, &
       -0.1124E-02_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip=20,20)/ &
       -0.1193E-10_r8, -0.2380E-10_r8, -0.4748E-10_r8, -0.9465E-10_r8, -0.1886E-09_r8, &
       -0.3751E-09_r8, -0.7436E-09_r8, -0.1466E-08_r8, -0.2872E-08_r8, -0.5508E-08_r8, &
       -0.1038E-07_r8, -0.1891E-07_r8, -0.3279E-07_r8, -0.5420E-07_r8, -0.8711E-07_r8, &
       -0.1403E-06_r8, -0.2248E-06_r8, -0.3221E-06_r8, -0.3459E-06_r8, -0.1066E-06_r8, &
        0.5938E-06_r8, & 
       ! data ((o1(ip,iw),iw=1,21)_r8, ip=21,21)/ &
        0.99999487_r8,  0.99998873_r8,  0.99997509_r8,  0.99994779_r8,  0.99989349_r8, &
        0.99978501_r8,  0.99956918_r8,  0.99913877_r8,  0.99828798_r8,  0.99661303_r8, &
        0.99335998_r8,  0.98715001_r8,  0.97566003_r8,  0.95530999_r8,  0.92153001_r8, &
        0.87000000_r8,  0.79869998_r8,  0.70819998_r8,  0.60109997_r8,  0.48110002_r8, &
        0.35140002_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip=21,21)/ &
       -0.2654E-08_r8, -0.5296E-08_r8, -0.1057E-07_r8, -0.2111E-07_r8, -0.4219E-07_r8, &
       -0.8445E-07_r8, -0.1696E-06_r8, -0.3428E-06_r8, -0.7013E-06_r8, -0.1458E-05_r8, &
       -0.3084E-05_r8, -0.6678E-05_r8, -0.1476E-04_r8, -0.3284E-04_r8, -0.7173E-04_r8, &
       -0.1481E-03_r8, -0.2786E-03_r8, -0.4688E-03_r8, -0.7052E-03_r8, -0.9506E-03_r8, &
       -0.1148E-02_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip=21,21)/ &
       -0.1195E-10_r8, -0.2384E-10_r8, -0.4755E-10_r8, -0.9482E-10_r8, -0.1890E-09_r8, &
       -0.3761E-09_r8, -0.7469E-09_r8, -0.1476E-08_r8, -0.2892E-08_r8, -0.5603E-08_r8, &
       -0.1060E-07_r8, -0.1942E-07_r8, -0.3393E-07_r8, -0.5508E-07_r8, -0.8290E-07_r8, &
       -0.1182E-06_r8, -0.1657E-06_r8, -0.2170E-06_r8, -0.1997E-06_r8,  0.6227E-07_r8, &
        0.7847E-06_r8, & 
       ! data ((o1(ip,iw),iw=1,21)_r8, ip=22,22)/ &
        0.99999541_r8,  0.99998873_r8,  0.99997497_r8,  0.99994737_r8,  0.99989349_r8, &
        0.99978501_r8,  0.99956882_r8,  0.99913770_r8,  0.99828303_r8,  0.99659699_r8, &
        0.99330199_r8,  0.98697001_r8,  0.97510999_r8,  0.95372999_r8,  0.91742998_r8, &
        0.86080003_r8,  0.78139997_r8,  0.68220001_r8,  0.56920004_r8,  0.44809997_r8, &
        0.32080001_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip=22,22)/ &
       -0.2653E-08_r8, -0.5295E-08_r8, -0.1057E-07_r8, -0.2110E-07_r8, -0.4215E-07_r8, &
       -0.8430E-07_r8, -0.1690E-06_r8, -0.3403E-06_r8, -0.6919E-06_r8, -0.1427E-05_r8, &
       -0.2991E-05_r8, -0.6399E-05_r8, -0.1398E-04_r8, -0.3099E-04_r8, -0.6824E-04_r8, &
       -0.1441E-03_r8, -0.2795E-03_r8, -0.4814E-03_r8, -0.7282E-03_r8, -0.9739E-03_r8, &
       -0.1163E-02_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip=22,22)/ &
       -0.1195E-10_r8, -0.2384E-10_r8, -0.4756E-10_r8, -0.9485E-10_r8, -0.1891E-09_r8, &
       -0.3765E-09_r8, -0.7483E-09_r8, -0.1481E-08_r8, -0.2908E-08_r8, -0.5660E-08_r8, &
       -0.1075E-07_r8, -0.1980E-07_r8, -0.3472E-07_r8, -0.5626E-07_r8, -0.8149E-07_r8, &
       -0.1027E-06_r8, -0.1136E-06_r8, -0.1071E-06_r8, -0.2991E-07_r8,  0.2743E-06_r8, &
        0.1017E-05_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip=23,23)/ &
        0.99999595_r8,  0.99998885_r8,  0.99997479_r8,  0.99994725_r8,  0.99989331_r8, &
        0.99978518_r8,  0.99956882_r8,  0.99913692_r8,  0.99827999_r8,  0.99658602_r8, &
        0.99326497_r8,  0.98685002_r8,  0.97474003_r8,  0.95260000_r8,  0.91441000_r8, &
        0.85360003_r8,  0.76719999_r8,  0.65990001_r8,  0.54190004_r8,  0.42119998_r8, &
        0.29699999_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip=23,23)/ &
       -0.2653E-08_r8, -0.5294E-08_r8, -0.1057E-07_r8, -0.2109E-07_r8, -0.4212E-07_r8, &
       -0.8420E-07_r8, -0.1686E-06_r8, -0.3388E-06_r8, -0.6858E-06_r8, -0.1406E-05_r8, &
       -0.2928E-05_r8, -0.6206E-05_r8, -0.1344E-04_r8, -0.2961E-04_r8, -0.6533E-04_r8, &
       -0.1399E-03_r8, -0.2780E-03_r8, -0.4904E-03_r8, -0.7488E-03_r8, -0.9953E-03_r8, &
       -0.1175E-02_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip=23,23)/ &
       -0.1195E-10_r8, -0.2384E-10_r8, -0.4756E-10_r8, -0.9485E-10_r8, -0.1891E-09_r8, &
       -0.3767E-09_r8, -0.7492E-09_r8, -0.1485E-08_r8, -0.2924E-08_r8, -0.5671E-08_r8, &
       -0.1084E-07_r8, -0.2009E-07_r8, -0.3549E-07_r8, -0.5773E-07_r8, -0.8208E-07_r8, &
       -0.9394E-07_r8, -0.7270E-07_r8, -0.3947E-08_r8,  0.1456E-06_r8,  0.5083E-06_r8, &
        0.1270E-05_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip=24,24)/ &
        0.99999630_r8,  0.99998873_r8,  0.99997401_r8,  0.99994725_r8,  0.99989349_r8, &
        0.99978501_r8,  0.99956959_r8,  0.99913663_r8,  0.99827701_r8,  0.99658000_r8, &
        0.99324101_r8,  0.98676997_r8,  0.97447002_r8,  0.95185000_r8,  0.91232002_r8, &
        0.84850001_r8,  0.75660002_r8,  0.64230001_r8,  0.52030003_r8,  0.40090001_r8, &
        0.27980000_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip=24,24)/ &
       -0.2653E-08_r8, -0.5294E-08_r8, -0.1056E-07_r8, -0.2109E-07_r8, -0.4210E-07_r8, &
       -0.8413E-07_r8, -0.1684E-06_r8, -0.3379E-06_r8, -0.6820E-06_r8, -0.1393E-05_r8, &
       -0.2889E-05_r8, -0.6080E-05_r8, -0.1307E-04_r8, -0.2861E-04_r8, -0.6310E-04_r8, &
       -0.1363E-03_r8, -0.2758E-03_r8, -0.4969E-03_r8, -0.7681E-03_r8, -0.1017E-02_r8, &
       -0.1186E-02_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip=24,24)/ &
       -0.1195E-10_r8, -0.2384E-10_r8, -0.4756E-10_r8, -0.9485E-10_r8, -0.1891E-09_r8, &
       -0.3768E-09_r8, -0.7497E-09_r8, -0.1487E-08_r8, -0.2933E-08_r8, -0.5710E-08_r8, &
       -0.1089E-07_r8, -0.2037E-07_r8, -0.3616E-07_r8, -0.5907E-07_r8, -0.8351E-07_r8, &
       -0.8925E-07_r8, -0.4122E-07_r8,  0.8779E-07_r8,  0.3143E-06_r8,  0.7281E-06_r8, &
        0.1500E-05_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip=25,25)/ &
        0.99999648_r8,  0.99998897_r8,  0.99997377_r8,  0.99994749_r8,  0.99989331_r8, &
        0.99978501_r8,  0.99956989_r8,  0.99913692_r8,  0.99827600_r8,  0.99657297_r8, &
        0.99322498_r8,  0.98672003_r8,  0.97431999_r8,  0.95137000_r8,  0.91095001_r8, &
        0.84500003_r8,  0.74909997_r8,  0.62979996_r8,  0.50510001_r8,  0.38679999_r8, &
        0.26789999_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip=25,25)/ &
       -0.2653E-08_r8, -0.5293E-08_r8, -0.1056E-07_r8, -0.2108E-07_r8, -0.4209E-07_r8, &
       -0.8409E-07_r8, -0.1682E-06_r8, -0.3373E-06_r8, -0.6797E-06_r8, -0.1383E-05_r8, &
       -0.2862E-05_r8, -0.5993E-05_r8, -0.1283E-04_r8, -0.2795E-04_r8, -0.6158E-04_r8, &
       -0.1338E-03_r8, -0.2743E-03_r8, -0.5030E-03_r8, -0.7863E-03_r8, -0.1038E-02_r8, &
       -0.1196E-02_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip=25,25)/ &
       -0.1195E-10_r8, -0.2383E-10_r8, -0.4755E-10_r8, -0.9484E-10_r8, -0.1891E-09_r8, &
       -0.3768E-09_r8, -0.7499E-09_r8, -0.1489E-08_r8, -0.2939E-08_r8, -0.5741E-08_r8, &
       -0.1100E-07_r8, -0.2066E-07_r8, -0.3660E-07_r8, -0.6002E-07_r8, -0.8431E-07_r8, &
       -0.8556E-07_r8, -0.1674E-07_r8,  0.1638E-06_r8,  0.4525E-06_r8,  0.8949E-06_r8, &
        0.1669E-05_r8, &
       ! data ((o1(ip,iw),iw=1,21)_r8, ip=26,26)/ &
        0.99999672_r8,  0.99998909_r8,  0.99997377_r8,  0.99994695_r8,  0.99989349_r8, &
        0.99978518_r8,  0.99956989_r8,  0.99913692_r8,  0.99827498_r8,  0.99657100_r8, &
        0.99321902_r8,  0.98668998_r8,  0.97421002_r8,  0.95106000_r8,  0.91009998_r8, &
        0.84280002_r8,  0.74430001_r8,  0.62180001_r8,  0.49519998_r8,  0.37800002_r8, &
        0.25999999_r8, &
       ! data ((o2(ip,iw),iw=1,21)_r8, ip=26,26)/ &
       -0.2652E-08_r8, -0.5292E-08_r8, -0.1056E-07_r8, -0.2108E-07_r8, -0.4208E-07_r8, &
       -0.8406E-07_r8, -0.1681E-06_r8, -0.3369E-06_r8, -0.6784E-06_r8, -0.1378E-05_r8, &
       -0.2843E-05_r8, -0.5944E-05_r8, -0.1269E-04_r8, -0.2759E-04_r8, -0.6078E-04_r8, &
       -0.1326E-03_r8, -0.2742E-03_r8, -0.5088E-03_r8, -0.8013E-03_r8, -0.1054E-02_r8, &
       -0.1202E-02_r8, &
       ! data ((o3(ip,iw),iw=1,21)_r8, ip=26,26)/ &
       -0.1194E-10_r8, -0.2383E-10_r8, -0.4754E-10_r8, -0.9482E-10_r8, -0.1891E-09_r8, &
       -0.3768E-09_r8, -0.7499E-09_r8, -0.1489E-08_r8, -0.2941E-08_r8, -0.5752E-08_r8, &
       -0.1104E-07_r8, -0.2069E-07_r8, -0.3661E-07_r8, -0.6012E-07_r8, -0.8399E-07_r8, &
       -0.8183E-07_r8,  0.1930E-08_r8,  0.2167E-06_r8,  0.5434E-06_r8,  0.9990E-06_r8, &
        0.1787E-05_r8/),SHAPE=(/nx*3*no/))

   it=0
   DO k=1,nx
      DO i=1,3
         DO j=1,no
            it=it+1
             IF(i==1) THEN
                !WRITE(*,'(a5,2e17.9) ' )'o1   ',data2(it),o1(k,j)
                o1(k,j)=data2(it)
             END IF
             IF(i==2)THEN
                !WRITE(*,'(a5,2e17.9) ' )'o2   ',data2(it),o2(k,j)
                o2(k,j)=data2(it)
             END IF 
             IF(i==3) THEN
                !WRITE(*,'(a5,2e17.9) ' )'o3   ',data2(it),o3(k,j)
                o3(k,j) =data2(it)
             END IF
         END DO
      END DO
   END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      data3(1:nx*3*nh)=RESHAPE(SOURCE=(/&
       ! data ((h11(ip,iw),iw=1,31)_r8, ip= 1_r8, 1)/  &
        0.99993843_r8,  0.99990183_r8,  0.99985260_r8,  0.99979079_r8,  0.99971771_r8,  &
        0.99963379_r8,  0.99953848_r8,  0.99942899_r8,  0.99930018_r8,  0.99914461_r8,  &
        0.99895102_r8,  0.99870503_r8,  0.99838799_r8,  0.99797899_r8,  0.99745202_r8,  &
        0.99677002_r8,  0.99587703_r8,  0.99469399_r8,  0.99311298_r8,  0.99097902_r8,  &
        0.98807001_r8,  0.98409998_r8,  0.97864997_r8,  0.97114998_r8,  0.96086001_r8,  &
        0.94682997_r8,  0.92777002_r8,  0.90200001_r8,  0.86739999_r8,  0.82169998_r8,  &
        0.76270002_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip= 1_r8, 1)/  &
       -0.2021E-06_r8, -0.3628E-06_r8, -0.5891E-06_r8, -0.8735E-06_r8, -0.1204E-05_r8,  &
       -0.1579E-05_r8, -0.2002E-05_r8, -0.2494E-05_r8, -0.3093E-05_r8, -0.3852E-05_r8,  &
       -0.4835E-05_r8, -0.6082E-05_r8, -0.7591E-05_r8, -0.9332E-05_r8, -0.1128E-04_r8,  &
       -0.1347E-04_r8, -0.1596E-04_r8, -0.1890E-04_r8, -0.2241E-04_r8, -0.2672E-04_r8,  &
       -0.3208E-04_r8, -0.3884E-04_r8, -0.4747E-04_r8, -0.5854E-04_r8, -0.7272E-04_r8,  &
       -0.9092E-04_r8, -0.1146E-03_r8, -0.1458E-03_r8, -0.1877E-03_r8, -0.2435E-03_r8,  &
       -0.3159E-03_r8,  &  
       ! data ((h13(ip,iw),iw=1,31)_r8, ip= 1_r8, 1)/  &
        0.5907E-09_r8,  0.8541E-09_r8,  0.1095E-08_r8,  0.1272E-08_r8,  0.1297E-08_r8,  &
        0.1105E-08_r8,  0.6788E-09_r8, -0.5585E-10_r8, -0.1147E-08_r8, -0.2746E-08_r8,  &
       -0.5001E-08_r8, -0.7715E-08_r8, -0.1037E-07_r8, -0.1227E-07_r8, -0.1287E-07_r8,  &
       -0.1175E-07_r8, -0.8517E-08_r8, -0.2920E-08_r8,  0.4786E-08_r8,  0.1407E-07_r8,  &
        0.2476E-07_r8,  0.3781E-07_r8,  0.5633E-07_r8,  0.8578E-07_r8,  0.1322E-06_r8,  &
        0.2013E-06_r8,  0.3006E-06_r8,  0.4409E-06_r8,  0.6343E-06_r8,  0.8896E-06_r8,  &
        0.1216E-05_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip= 2_r8, 2)/  &
        0.99993837_r8,  0.99990171_r8,  0.99985230_r8,  0.99979031_r8,  0.99971670_r8,  &
        0.99963200_r8,  0.99953520_r8,  0.99942321_r8,  0.99928987_r8,  0.99912637_r8,  &
        0.99892002_r8,  0.99865198_r8,  0.99830002_r8,  0.99783802_r8,  0.99723297_r8,  &
        0.99643701_r8,  0.99537897_r8,  0.99396098_r8,  0.99204701_r8,  0.98944002_r8,  &
        0.98588002_r8,  0.98098999_r8,  0.97425997_r8,  0.96502000_r8,  0.95236999_r8,  &
        0.93515998_r8,  0.91184998_r8,  0.88040000_r8,  0.83859998_r8,  0.78429997_r8,  &
        0.71560001_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip= 2_r8, 2)/  &
       -0.2017E-06_r8, -0.3620E-06_r8, -0.5878E-06_r8, -0.8713E-06_r8, -0.1201E-05_r8,  &
       -0.1572E-05_r8, -0.1991E-05_r8, -0.2476E-05_r8, -0.3063E-05_r8, -0.3808E-05_r8,  &
       -0.4776E-05_r8, -0.6011E-05_r8, -0.7516E-05_r8, -0.9272E-05_r8, -0.1127E-04_r8,  &
       -0.1355E-04_r8, -0.1620E-04_r8, -0.1936E-04_r8, -0.2321E-04_r8, -0.2797E-04_r8,  &
       -0.3399E-04_r8, -0.4171E-04_r8, -0.5172E-04_r8, -0.6471E-04_r8, -0.8150E-04_r8,  &
       -0.1034E-03_r8, -0.1321E-03_r8, -0.1705E-03_r8, -0.2217E-03_r8, -0.2889E-03_r8,  &
       -0.3726E-03_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip= 2_r8, 2)/  &
        0.5894E-09_r8,  0.8519E-09_r8,  0.1092E-08_r8,  0.1267E-08_r8,  0.1289E-08_r8,  &
        0.1093E-08_r8,  0.6601E-09_r8, -0.7831E-10_r8, -0.1167E-08_r8, -0.2732E-08_r8,  &
       -0.4864E-08_r8, -0.7334E-08_r8, -0.9581E-08_r8, -0.1097E-07_r8, -0.1094E-07_r8,  &
       -0.8999E-08_r8, -0.4669E-08_r8,  0.2391E-08_r8,  0.1215E-07_r8,  0.2424E-07_r8,  &
        0.3877E-07_r8,  0.5711E-07_r8,  0.8295E-07_r8,  0.1218E-06_r8,  0.1793E-06_r8,  &
        0.2621E-06_r8,  0.3812E-06_r8,  0.5508E-06_r8,  0.7824E-06_r8,  0.1085E-05_r8,  &
        0.1462E-05_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip= 3_r8, 3)/ &
        0.99993825_r8,  0.99990153_r8,  0.99985188_r8,  0.99978942_r8,  0.99971509_r8, &
        0.99962920_r8,  0.99953020_r8,  0.99941432_r8,  0.99927431_r8,  0.99909937_r8, &
        0.99887401_r8,  0.99857497_r8,  0.99817699_r8,  0.99764699_r8,  0.99694097_r8, &
        0.99599802_r8,  0.99473000_r8,  0.99301600_r8,  0.99068397_r8,  0.98749000_r8, &
        0.98311001_r8,  0.97707999_r8,  0.96877003_r8,  0.95738000_r8,  0.94186002_r8, &
        0.92079002_r8,  0.89230001_r8,  0.85420001_r8,  0.80430001_r8,  0.74049997_r8, &
        0.66200000_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip= 3_r8, 3)/ &
       -0.2011E-06_r8, -0.3609E-06_r8, -0.5859E-06_r8, -0.8680E-06_r8, -0.1195E-05_r8, &
       -0.1563E-05_r8, -0.1975E-05_r8, -0.2450E-05_r8, -0.3024E-05_r8, -0.3755E-05_r8, &
       -0.4711E-05_r8, -0.5941E-05_r8, -0.7455E-05_r8, -0.9248E-05_r8, -0.1132E-04_r8, &
       -0.1373E-04_r8, -0.1659E-04_r8, -0.2004E-04_r8, -0.2431E-04_r8, -0.2966E-04_r8, &
       -0.3653E-04_r8, -0.4549E-04_r8, -0.5724E-04_r8, -0.7259E-04_r8, -0.9265E-04_r8, &
       -0.1191E-03_r8, -0.1543E-03_r8, -0.2013E-03_r8, -0.2633E-03_r8, -0.3421E-03_r8, &
       -0.4350E-03_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip= 3_r8, 3)/ &
        0.5872E-09_r8,  0.8484E-09_r8,  0.1087E-08_r8,  0.1259E-08_r8,  0.1279E-08_r8, &
        0.1077E-08_r8,  0.6413E-09_r8, -0.9334E-10_r8, -0.1161E-08_r8, -0.2644E-08_r8, &
       -0.4588E-08_r8, -0.6709E-08_r8, -0.8474E-08_r8, -0.9263E-08_r8, -0.8489E-08_r8, &
       -0.5553E-08_r8,  0.1203E-09_r8,  0.9035E-08_r8,  0.2135E-07_r8,  0.3689E-07_r8, &
        0.5610E-07_r8,  0.8097E-07_r8,  0.1155E-06_r8,  0.1649E-06_r8,  0.2350E-06_r8, &
        0.3353E-06_r8,  0.4806E-06_r8,  0.6858E-06_r8,  0.9617E-06_r8,  0.1315E-05_r8, &
        0.1741E-05_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip= 4_r8, 4)/ &
        0.99993813_r8,  0.99990118_r8,  0.99985123_r8,  0.99978811_r8,  0.99971271_r8, &
        0.99962479_r8,  0.99952239_r8,  0.99940068_r8,  0.99925101_r8,  0.99905968_r8, &
        0.99880803_r8,  0.99846900_r8,  0.99800998_r8,  0.99738997_r8,  0.99655402_r8, &
        0.99542397_r8,  0.99389100_r8,  0.99180400_r8,  0.98895001_r8,  0.98501998_r8, &
        0.97961003_r8,  0.97215003_r8,  0.96191001_r8,  0.94791001_r8,  0.92887998_r8, &
        0.90311998_r8,  0.86849999_r8,  0.82270002_r8,  0.76370001_r8,  0.69000000_r8, &
        0.60240000_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip= 4_r8, 4)/ &
       -0.2001E-06_r8, -0.3592E-06_r8, -0.5829E-06_r8, -0.8631E-06_r8, -0.1187E-05_r8, &
       -0.1549E-05_r8, -0.1953E-05_r8, -0.2415E-05_r8, -0.2975E-05_r8, -0.3694E-05_r8, &
       -0.4645E-05_r8, -0.5882E-05_r8, -0.7425E-05_r8, -0.9279E-05_r8, -0.1147E-04_r8, &
       -0.1406E-04_r8, -0.1717E-04_r8, -0.2100E-04_r8, -0.2580E-04_r8, -0.3191E-04_r8, &
       -0.3989E-04_r8, -0.5042E-04_r8, -0.6432E-04_r8, -0.8261E-04_r8, -0.1068E-03_r8, &
       -0.1389E-03_r8, -0.1820E-03_r8, -0.2391E-03_r8, -0.3127E-03_r8, -0.4021E-03_r8, &
       -0.5002E-03_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip= 4_r8, 4)/ &
        0.5838E-09_r8,  0.8426E-09_r8,  0.1081E-08_r8,  0.1249E-08_r8,  0.1267E-08_r8, &
        0.1062E-08_r8,  0.6313E-09_r8, -0.8241E-10_r8, -0.1094E-08_r8, -0.2436E-08_r8, &
       -0.4100E-08_r8, -0.5786E-08_r8, -0.6992E-08_r8, -0.7083E-08_r8, -0.5405E-08_r8, &
       -0.1259E-08_r8,  0.6099E-08_r8,  0.1732E-07_r8,  0.3276E-07_r8,  0.5256E-07_r8, &
        0.7756E-07_r8,  0.1103E-06_r8,  0.1547E-06_r8,  0.2159E-06_r8,  0.3016E-06_r8, &
        0.4251E-06_r8,  0.6033E-06_r8,  0.8499E-06_r8,  0.1175E-05_r8,  0.1579E-05_r8, &
        0.2044E-05_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip= 5_r8, 5)/ &
        0.99993789_r8,  0.99990070_r8,  0.99985009_r8,  0.99978602_r8,  0.99970889_r8, &
        0.99961799_r8,  0.99951053_r8,  0.99938041_r8,  0.99921662_r8,  0.99900270_r8, &
        0.99871498_r8,  0.99832201_r8,  0.99778402_r8,  0.99704897_r8,  0.99604702_r8, &
        0.99468100_r8,  0.99281400_r8,  0.99025702_r8,  0.98673999_r8,  0.98189002_r8, &
        0.97521001_r8,  0.96600002_r8,  0.95337999_r8,  0.93620998_r8,  0.91292000_r8, &
        0.88150001_r8,  0.83969998_r8,  0.78530002_r8,  0.71650004_r8,  0.63330001_r8, &
        0.53799999_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip= 5_r8, 5)/ &
       -0.1987E-06_r8, -0.3565E-06_r8, -0.5784E-06_r8, -0.8557E-06_r8, -0.1175E-05_r8, &
       -0.1530E-05_r8, -0.1923E-05_r8, -0.2372E-05_r8, -0.2919E-05_r8, -0.3631E-05_r8, &
       -0.4587E-05_r8, -0.5848E-05_r8, -0.7442E-05_r8, -0.9391E-05_r8, -0.1173E-04_r8, &
       -0.1455E-04_r8, -0.1801E-04_r8, -0.2232E-04_r8, -0.2779E-04_r8, -0.3489E-04_r8, &
       -0.4428E-04_r8, -0.5678E-04_r8, -0.7333E-04_r8, -0.9530E-04_r8, -0.1246E-03_r8, &
       -0.1639E-03_r8, -0.2164E-03_r8, -0.2848E-03_r8, -0.3697E-03_r8, -0.4665E-03_r8, &
       -0.5646E-03_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip= 5_r8, 5)/ &
        0.5785E-09_r8,  0.8338E-09_r8,  0.1071E-08_r8,  0.1239E-08_r8,  0.1256E-08_r8, &
        0.1057E-08_r8,  0.6480E-09_r8, -0.1793E-10_r8, -0.9278E-09_r8, -0.2051E-08_r8, &
       -0.3337E-08_r8, -0.4514E-08_r8, -0.5067E-08_r8, -0.4328E-08_r8, -0.1545E-08_r8, &
        0.4100E-08_r8,  0.1354E-07_r8,  0.2762E-07_r8,  0.4690E-07_r8,  0.7190E-07_r8, &
        0.1040E-06_r8,  0.1459E-06_r8,  0.2014E-06_r8,  0.2764E-06_r8,  0.3824E-06_r8, &
        0.5359E-06_r8,  0.7532E-06_r8,  0.1047E-05_r8,  0.1424E-05_r8,  0.1873E-05_r8, &
        0.2356E-05_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip= 6_r8, 6)/ &
        0.99993753_r8,  0.99989992_r8,  0.99984848_r8,  0.99978292_r8,  0.99970299_r8, &
        0.99960762_r8,  0.99949282_r8,  0.99935049_r8,  0.99916708_r8,  0.99892199_r8, &
        0.99858701_r8,  0.99812400_r8,  0.99748403_r8,  0.99660099_r8,  0.99538797_r8, &
        0.99372399_r8,  0.99143797_r8,  0.98829001_r8,  0.98395002_r8,  0.97794998_r8, &
        0.96968001_r8,  0.95832998_r8,  0.94283003_r8,  0.92179000_r8,  0.89330000_r8, &
        0.85530001_r8,  0.80519998_r8,  0.74140000_r8,  0.66280001_r8,  0.57099998_r8, &
        0.47049999_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip= 6_r8, 6)/ &
       -0.1964E-06_r8, -0.3526E-06_r8, -0.5717E-06_r8, -0.8451E-06_r8, -0.1158E-05_r8, &
       -0.1504E-05_r8, -0.1886E-05_r8, -0.2322E-05_r8, -0.2861E-05_r8, -0.3576E-05_r8, &
       -0.4552E-05_r8, -0.5856E-05_r8, -0.7529E-05_r8, -0.9609E-05_r8, -0.1216E-04_r8, &
       -0.1528E-04_r8, -0.1916E-04_r8, -0.2408E-04_r8, -0.3043E-04_r8, -0.3880E-04_r8, &
       -0.4997E-04_r8, -0.6488E-04_r8, -0.8474E-04_r8, -0.1113E-03_r8, -0.1471E-03_r8, &
       -0.1950E-03_r8, -0.2583E-03_r8, -0.3384E-03_r8, -0.4326E-03_r8, -0.5319E-03_r8, &
       -0.6244E-03_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip= 6_r8, 6)/ &
        0.5713E-09_r8,  0.8263E-09_r8,  0.1060E-08_r8,  0.1226E-08_r8,  0.1252E-08_r8, &
        0.1076E-08_r8,  0.7149E-09_r8,  0.1379E-09_r8, -0.6043E-09_r8, -0.1417E-08_r8, &
       -0.2241E-08_r8, -0.2830E-08_r8, -0.2627E-08_r8, -0.8950E-09_r8,  0.3231E-08_r8, &
        0.1075E-07_r8,  0.2278E-07_r8,  0.4037E-07_r8,  0.6439E-07_r8,  0.9576E-07_r8, &
        0.1363E-06_r8,  0.1886E-06_r8,  0.2567E-06_r8,  0.3494E-06_r8,  0.4821E-06_r8, &
        0.6719E-06_r8,  0.9343E-06_r8,  0.1280E-05_r8,  0.1705E-05_r8,  0.2184E-05_r8, &
        0.2651E-05_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip= 7_r8, 7)/ &
        0.99993700_r8,  0.99989867_r8,  0.99984592_r8,  0.99977797_r8,  0.99969423_r8, &
        0.99959219_r8,  0.99946660_r8,  0.99930722_r8,  0.99909681_r8,  0.99880999_r8, &
        0.99841303_r8,  0.99786001_r8,  0.99708802_r8,  0.99601799_r8,  0.99453998_r8, &
        0.99250001_r8,  0.98969001_r8,  0.98580003_r8,  0.98041999_r8,  0.97299999_r8, &
        0.96279001_r8,  0.94881999_r8,  0.92980999_r8,  0.90407002_r8,  0.86949998_r8, &
        0.82370001_r8,  0.76459998_r8,  0.69089997_r8,  0.60310000_r8,  0.50479996_r8, &
        0.40219998_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip= 7_r8, 7)/ &
       -0.1932E-06_r8, -0.3467E-06_r8, -0.5623E-06_r8, -0.8306E-06_r8, -0.1136E-05_r8, &
       -0.1472E-05_r8, -0.1842E-05_r8, -0.2269E-05_r8, -0.2807E-05_r8, -0.3539E-05_r8, &
       -0.4553E-05_r8, -0.5925E-05_r8, -0.7710E-05_r8, -0.9968E-05_r8, -0.1278E-04_r8, &
       -0.1629E-04_r8, -0.2073E-04_r8, -0.2644E-04_r8, -0.3392E-04_r8, -0.4390E-04_r8, &
       -0.5727E-04_r8, -0.7516E-04_r8, -0.9916E-04_r8, -0.1315E-03_r8, -0.1752E-03_r8, &
       -0.2333E-03_r8, -0.3082E-03_r8, -0.3988E-03_r8, -0.4982E-03_r8, -0.5947E-03_r8, &
       -0.6764E-03_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip= 7_r8, 7)/ &
        0.5612E-09_r8,  0.8116E-09_r8,  0.1048E-08_r8,  0.1222E-08_r8,  0.1270E-08_r8, &
        0.1141E-08_r8,  0.8732E-09_r8,  0.4336E-09_r8, -0.6548E-10_r8, -0.4774E-09_r8, &
       -0.7556E-09_r8, -0.6577E-09_r8,  0.4377E-09_r8,  0.3359E-08_r8,  0.9159E-08_r8, &
        0.1901E-07_r8,  0.3422E-07_r8,  0.5616E-07_r8,  0.8598E-07_r8,  0.1251E-06_r8, &
        0.1752E-06_r8,  0.2392E-06_r8,  0.3228E-06_r8,  0.4389E-06_r8,  0.6049E-06_r8, &
        0.8370E-06_r8,  0.1150E-05_r8,  0.1547E-05_r8,  0.2012E-05_r8,  0.2493E-05_r8, &
        0.2913E-05_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip= 8_r8, 8)/  &
        0.99993622_r8,  0.99989682_r8,  0.99984211_r8,  0.99977070_r8,  0.99968100_r8,  &
        0.99956948_r8,  0.99942881_r8,  0.99924588_r8,  0.99899900_r8,  0.99865800_r8,  &
        0.99818099_r8,  0.99751103_r8,  0.99657297_r8,  0.99526602_r8,  0.99345201_r8,  &
        0.99094099_r8,  0.98746002_r8,  0.98264998_r8,  0.97599000_r8,  0.96682000_r8,  &
        0.95423001_r8,  0.93708003_r8,  0.91380000_r8,  0.88239998_r8,  0.84060001_r8,  &
        0.78610003_r8,  0.71730000_r8,  0.63400000_r8,  0.53859997_r8,  0.43660003_r8,  &
        0.33510000_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip= 8_r8, 8)/  &
       -0.1885E-06_r8, -0.3385E-06_r8, -0.5493E-06_r8, -0.8114E-06_r8, -0.1109E-05_r8,  &
       -0.1436E-05_r8, -0.1796E-05_r8, -0.2219E-05_r8, -0.2770E-05_r8, -0.3535E-05_r8,  &
       -0.4609E-05_r8, -0.6077E-05_r8, -0.8016E-05_r8, -0.1051E-04_r8, -0.1367E-04_r8,  &
       -0.1768E-04_r8, -0.2283E-04_r8, -0.2955E-04_r8, -0.3849E-04_r8, -0.5046E-04_r8,  &
       -0.6653E-04_r8, -0.8813E-04_r8, -0.1173E-03_r8, -0.1569E-03_r8, -0.2100E-03_r8,  &
       -0.2794E-03_r8, -0.3656E-03_r8, -0.4637E-03_r8, -0.5629E-03_r8, -0.6512E-03_r8,  &
       -0.7167E-03_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip= 8_r8, 8)/  &
        0.5477E-09_r8,  0.8000E-09_r8,  0.1039E-08_r8,  0.1234E-08_r8,  0.1331E-08_r8,  &
        0.1295E-08_r8,  0.1160E-08_r8,  0.9178E-09_r8,  0.7535E-09_r8,  0.8301E-09_r8,  &
        0.1184E-08_r8,  0.2082E-08_r8,  0.4253E-08_r8,  0.8646E-08_r8,  0.1650E-07_r8,  &
        0.2920E-07_r8,  0.4834E-07_r8,  0.7564E-07_r8,  0.1125E-06_r8,  0.1606E-06_r8,  &
        0.2216E-06_r8,  0.2992E-06_r8,  0.4031E-06_r8,  0.5493E-06_r8,  0.7549E-06_r8,  &
        0.1035E-05_r8,  0.1400E-05_r8,  0.1843E-05_r8,  0.2327E-05_r8,  0.2774E-05_r8,  &
        0.3143E-05_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip= 9_r8, 9)/  &
        0.99993503_r8,  0.99989408_r8,  0.99983650_r8,  0.99975997_r8,  0.99966192_r8,  &
        0.99953687_r8,  0.99937540_r8,  0.99916059_r8,  0.99886602_r8,  0.99845397_r8,  &
        0.99787402_r8,  0.99705601_r8,  0.99590701_r8,  0.99430102_r8,  0.99206603_r8,  &
        0.98896003_r8,  0.98465002_r8,  0.97869003_r8,  0.97044003_r8,  0.95911002_r8,  &
        0.94363999_r8,  0.92260998_r8,  0.89419997_r8,  0.85609996_r8,  0.80610001_r8,  &
        0.74220002_r8,  0.66359997_r8,  0.57169998_r8,  0.47100002_r8,  0.36860001_r8,  &
        0.27079999_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip= 9_r8, 9)/  &
       -0.1822E-06_r8, -0.3274E-06_r8, -0.5325E-06_r8, -0.7881E-06_r8, -0.1079E-05_r8,  &
       -0.1398E-05_r8, -0.1754E-05_r8, -0.2184E-05_r8, -0.2763E-05_r8, -0.3581E-05_r8,  &
       -0.4739E-05_r8, -0.6341E-05_r8, -0.8484E-05_r8, -0.1128E-04_r8, -0.1490E-04_r8,  &
       -0.1955E-04_r8, -0.2561E-04_r8, -0.3364E-04_r8, -0.4438E-04_r8, -0.5881E-04_r8,  &
       -0.7822E-04_r8, -0.1045E-03_r8, -0.1401E-03_r8, -0.1884E-03_r8, -0.2523E-03_r8,  &
       -0.3335E-03_r8, -0.4289E-03_r8, -0.5296E-03_r8, -0.6231E-03_r8, -0.6980E-03_r8,  &
       -0.7406E-03_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip= 9_r8, 9)/ &
        0.5334E-09_r8,  0.7859E-09_r8,  0.1043E-08_r8,  0.1279E-08_r8,  0.1460E-08_r8, &
        0.1560E-08_r8,  0.1618E-08_r8,  0.1657E-08_r8,  0.1912E-08_r8,  0.2569E-08_r8, &
        0.3654E-08_r8,  0.5509E-08_r8,  0.8964E-08_r8,  0.1518E-07_r8,  0.2560E-07_r8, &
        0.4178E-07_r8,  0.6574E-07_r8,  0.9958E-07_r8,  0.1449E-06_r8,  0.2031E-06_r8, &
        0.2766E-06_r8,  0.3718E-06_r8,  0.5022E-06_r8,  0.6849E-06_r8,  0.9360E-06_r8, &
        0.1268E-05_r8,  0.1683E-05_r8,  0.2157E-05_r8,  0.2625E-05_r8,  0.3020E-05_r8, &
        0.3364E-05_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip=10,10)/ &
        0.99993336_r8,  0.99989021_r8,  0.99982840_r8,  0.99974459_r8,  0.99963468_r8, &
        0.99949121_r8,  0.99930137_r8,  0.99904430_r8,  0.99868703_r8,  0.99818403_r8, &
        0.99747300_r8,  0.99646801_r8,  0.99505299_r8,  0.99307102_r8,  0.99030602_r8, &
        0.98645997_r8,  0.98111999_r8,  0.97372001_r8,  0.96353000_r8,  0.94957000_r8, &
        0.93058997_r8,  0.90486002_r8,  0.87029999_r8,  0.82449996_r8,  0.76530004_r8, &
        0.69159997_r8,  0.60380000_r8,  0.50529999_r8,  0.40259999_r8,  0.30269998_r8, &
        0.21020001_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip=10,10)/ &
       -0.1742E-06_r8, -0.3134E-06_r8, -0.5121E-06_r8, -0.7619E-06_r8, -0.1048E-05_r8, &
       -0.1364E-05_r8, -0.1725E-05_r8, -0.2177E-05_r8, -0.2801E-05_r8, -0.3694E-05_r8, &
       -0.4969E-05_r8, -0.6748E-05_r8, -0.9161E-05_r8, -0.1236E-04_r8, -0.1655E-04_r8, &
       -0.2203E-04_r8, -0.2927E-04_r8, -0.3894E-04_r8, -0.5192E-04_r8, -0.6936E-04_r8, &
       -0.9294E-04_r8, -0.1250E-03_r8, -0.1686E-03_r8, -0.2271E-03_r8, -0.3027E-03_r8, &
       -0.3944E-03_r8, -0.4951E-03_r8, -0.5928E-03_r8, -0.6755E-03_r8, -0.7309E-03_r8, &
       -0.7417E-03_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip=10,10)/ &
        0.5179E-09_r8,  0.7789E-09_r8,  0.1071E-08_r8,  0.1382E-08_r8,  0.1690E-08_r8, &
        0.1979E-08_r8,  0.2297E-08_r8,  0.2704E-08_r8,  0.3466E-08_r8,  0.4794E-08_r8, &
        0.6746E-08_r8,  0.9739E-08_r8,  0.1481E-07_r8,  0.2331E-07_r8,  0.3679E-07_r8, &
        0.5726E-07_r8,  0.8716E-07_r8,  0.1289E-06_r8,  0.1837E-06_r8,  0.2534E-06_r8, &
        0.3424E-06_r8,  0.4609E-06_r8,  0.6245E-06_r8,  0.8495E-06_r8,  0.1151E-05_r8, &
        0.1536E-05_r8,  0.1991E-05_r8,  0.2468E-05_r8,  0.2891E-05_r8,  0.3245E-05_r8, &
        0.3580E-05_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip=11,11)/ &
        0.99993110_r8,  0.99988490_r8,  0.99981719_r8,  0.99972337_r8,  0.99959719_r8, &
        0.99942869_r8,  0.99920130_r8,  0.99888903_r8,  0.99845201_r8,  0.99783301_r8, &
        0.99695599_r8,  0.99571502_r8,  0.99396503_r8,  0.99150997_r8,  0.98808002_r8, &
        0.98329997_r8,  0.97667003_r8,  0.96750998_r8,  0.95494002_r8,  0.93779999_r8, &
        0.91453999_r8,  0.88319999_r8,  0.84130001_r8,  0.78689998_r8,  0.71799999_r8, &
        0.63470000_r8,  0.53909999_r8,  0.43699998_r8,  0.33550000_r8,  0.24010003_r8, &
        0.15420002_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip=11,11)/ &
       -0.1647E-06_r8, -0.2974E-06_r8, -0.4900E-06_r8, -0.7358E-06_r8, -0.1022E-05_r8, &
       -0.1344E-05_r8, -0.1721E-05_r8, -0.2212E-05_r8, -0.2901E-05_r8, -0.3896E-05_r8, &
       -0.5327E-05_r8, -0.7342E-05_r8, -0.1011E-04_r8, -0.1382E-04_r8, -0.1875E-04_r8, &
       -0.2530E-04_r8, -0.3403E-04_r8, -0.4573E-04_r8, -0.6145E-04_r8, -0.8264E-04_r8, &
       -0.1114E-03_r8, -0.1507E-03_r8, -0.2039E-03_r8, -0.2737E-03_r8, -0.3607E-03_r8, &
       -0.4599E-03_r8, -0.5604E-03_r8, -0.6497E-03_r8, -0.7161E-03_r8, -0.7443E-03_r8, &
       -0.7133E-03_r8,  & 
       ! data ((h13(ip,iw),iw=1,31)_r8, ip=11,11)/ &
        0.5073E-09_r8,  0.7906E-09_r8,  0.1134E-08_r8,  0.1560E-08_r8,  0.2046E-08_r8, &
        0.2589E-08_r8,  0.3254E-08_r8,  0.4107E-08_r8,  0.5481E-08_r8,  0.7602E-08_r8, &
        0.1059E-07_r8,  0.1501E-07_r8,  0.2210E-07_r8,  0.3334E-07_r8,  0.5055E-07_r8, &
        0.7629E-07_r8,  0.1134E-06_r8,  0.1642E-06_r8,  0.2298E-06_r8,  0.3133E-06_r8, &
        0.4225E-06_r8,  0.5709E-06_r8,  0.7739E-06_r8,  0.1047E-05_r8,  0.1401E-05_r8, &
        0.1833E-05_r8,  0.2308E-05_r8,  0.2753E-05_r8,  0.3125E-05_r8,  0.3467E-05_r8, &
        0.3748E-05_r8,  & 
       ! data ((h11(ip,iw),iw=1,31)_r8, ip=12,12)/ &
        0.99992824_r8,  0.99987793_r8,  0.99980247_r8,  0.99969512_r8,  0.99954712_r8, &
        0.99934530_r8,  0.99906880_r8,  0.99868500_r8,  0.99814498_r8,  0.99738002_r8, &
        0.99629498_r8,  0.99475700_r8,  0.99258602_r8,  0.98953998_r8,  0.98527998_r8, &
        0.97934997_r8,  0.97112000_r8,  0.95981002_r8,  0.94433999_r8,  0.92332000_r8, &
        0.89490002_r8,  0.85680002_r8,  0.80680001_r8,  0.74290001_r8,  0.66420001_r8, &
        0.57220000_r8,  0.47149998_r8,  0.36900002_r8,  0.27109998_r8,  0.18159997_r8, &
        0.10460001_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip=12,12)/ &
       -0.1548E-06_r8, -0.2808E-06_r8, -0.4683E-06_r8, -0.7142E-06_r8, -0.1008E-05_r8, &
       -0.1347E-05_r8, -0.1758E-05_r8, -0.2306E-05_r8, -0.3083E-05_r8, -0.4214E-05_r8, &
       -0.5851E-05_r8, -0.8175E-05_r8, -0.1140E-04_r8, -0.1577E-04_r8, -0.2166E-04_r8, &
       -0.2955E-04_r8, -0.4014E-04_r8, -0.5434E-04_r8, -0.7343E-04_r8, -0.9931E-04_r8, &
       -0.1346E-03_r8, -0.1826E-03_r8, -0.2467E-03_r8, -0.3283E-03_r8, -0.4246E-03_r8, &
       -0.5264E-03_r8, -0.6211E-03_r8, -0.6970E-03_r8, -0.7402E-03_r8, -0.7316E-03_r8, &
       -0.6486E-03_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip=12,12)/ &
        0.5078E-09_r8,  0.8244E-09_r8,  0.1255E-08_r8,  0.1826E-08_r8,  0.2550E-08_r8, &
        0.3438E-08_r8,  0.4532E-08_r8,  0.5949E-08_r8,  0.8041E-08_r8,  0.1110E-07_r8, &
        0.1534E-07_r8,  0.2157E-07_r8,  0.3116E-07_r8,  0.4570E-07_r8,  0.6747E-07_r8, &
        0.9961E-07_r8,  0.1451E-06_r8,  0.2061E-06_r8,  0.2843E-06_r8,  0.3855E-06_r8, &
        0.5213E-06_r8,  0.7060E-06_r8,  0.9544E-06_r8,  0.1280E-05_r8,  0.1684E-05_r8, &
        0.2148E-05_r8,  0.2609E-05_r8,  0.3002E-05_r8,  0.3349E-05_r8,  0.3670E-05_r8, &
        0.3780E-05_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip=13,13)/ &
        0.99992472_r8,  0.99986941_r8,  0.99978399_r8,  0.99965900_r8,  0.99948251_r8, &
        0.99923742_r8,  0.99889702_r8,  0.99842298_r8,  0.99775398_r8,  0.99680400_r8, &
        0.99545598_r8,  0.99354500_r8,  0.99084800_r8,  0.98706001_r8,  0.98176998_r8, &
        0.97439998_r8,  0.96423000_r8,  0.95029002_r8,  0.93129998_r8,  0.90557003_r8, &
        0.87099999_r8,  0.82520002_r8,  0.76600003_r8,  0.69220001_r8,  0.60440004_r8, &
        0.50580001_r8,  0.40310001_r8,  0.30299997_r8,  0.21039999_r8,  0.12860000_r8, &
        0.06360000_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip=13,13)/ &
       -0.1461E-06_r8, -0.2663E-06_r8, -0.4512E-06_r8, -0.7027E-06_r8, -0.1014E-05_r8, &
       -0.1387E-05_r8, -0.1851E-05_r8, -0.2478E-05_r8, -0.3373E-05_r8, -0.4682E-05_r8, &
       -0.6588E-05_r8, -0.9311E-05_r8, -0.1311E-04_r8, -0.1834E-04_r8, -0.2544E-04_r8, &
       -0.3502E-04_r8, -0.4789E-04_r8, -0.6515E-04_r8, -0.8846E-04_r8, -0.1202E-03_r8, &
       -0.1635E-03_r8, -0.2217E-03_r8, -0.2975E-03_r8, -0.3897E-03_r8, -0.4913E-03_r8, &
       -0.5902E-03_r8, -0.6740E-03_r8, -0.7302E-03_r8, -0.7415E-03_r8, -0.6858E-03_r8, &
       -0.5447E-03_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip=13,13)/ &
        0.5236E-09_r8,  0.8873E-09_r8,  0.1426E-08_r8,  0.2193E-08_r8,  0.3230E-08_r8, &
        0.4555E-08_r8,  0.6200E-08_r8,  0.8298E-08_r8,  0.1126E-07_r8,  0.1544E-07_r8, &
        0.2130E-07_r8,  0.2978E-07_r8,  0.4239E-07_r8,  0.6096E-07_r8,  0.8829E-07_r8, &
        0.1280E-06_r8,  0.1830E-06_r8,  0.2555E-06_r8,  0.3493E-06_r8,  0.4740E-06_r8, &
        0.6431E-06_r8,  0.8701E-06_r8,  0.1169E-05_r8,  0.1547E-05_r8,  0.1992E-05_r8, &
        0.2460E-05_r8,  0.2877E-05_r8,  0.3230E-05_r8,  0.3569E-05_r8,  0.3782E-05_r8, &
        0.3591E-05_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip=14,14)/ &
        0.99992090_r8,  0.99985969_r8,  0.99976218_r8,  0.99961531_r8,  0.99940270_r8, &
        0.99910218_r8,  0.99868101_r8,  0.99809098_r8,  0.99725902_r8,  0.99607700_r8, &
        0.99440002_r8,  0.99202299_r8,  0.98866999_r8,  0.98395997_r8,  0.97737998_r8, &
        0.96825999_r8,  0.95570999_r8,  0.93857002_r8,  0.91531003_r8,  0.88389999_r8, &
        0.84210002_r8,  0.78759998_r8,  0.71869999_r8,  0.63530004_r8,  0.53970003_r8, &
        0.43750000_r8,  0.33590001_r8,  0.24040002_r8,  0.15439999_r8,  0.08300000_r8, &
        0.03299999_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip=14,14)/ &
       -0.1402E-06_r8, -0.2569E-06_r8, -0.4428E-06_r8, -0.7076E-06_r8, -0.1051E-05_r8, &
       -0.1478E-05_r8, -0.2019E-05_r8, -0.2752E-05_r8, -0.3802E-05_r8, -0.5343E-05_r8, &
       -0.7594E-05_r8, -0.1082E-04_r8, -0.1536E-04_r8, -0.2166E-04_r8, -0.3028E-04_r8, &
       -0.4195E-04_r8, -0.5761E-04_r8, -0.7867E-04_r8, -0.1072E-03_r8, -0.1462E-03_r8, &
       -0.1990E-03_r8, -0.2687E-03_r8, -0.3559E-03_r8, -0.4558E-03_r8, -0.5572E-03_r8, &
       -0.6476E-03_r8, -0.7150E-03_r8, -0.7439E-03_r8, -0.7133E-03_r8, -0.6015E-03_r8, &
       -0.4089E-03_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip=14,14)/ &
        0.5531E-09_r8,  0.9757E-09_r8,  0.1644E-08_r8,  0.2650E-08_r8,  0.4074E-08_r8, &
        0.5957E-08_r8,  0.8314E-08_r8,  0.1128E-07_r8,  0.1528E-07_r8,  0.2087E-07_r8, &
        0.2874E-07_r8,  0.4002E-07_r8,  0.5631E-07_r8,  0.7981E-07_r8,  0.1139E-06_r8, &
        0.1621E-06_r8,  0.2275E-06_r8,  0.3136E-06_r8,  0.4280E-06_r8,  0.5829E-06_r8, &
        0.7917E-06_r8,  0.1067E-05_r8,  0.1419E-05_r8,  0.1844E-05_r8,  0.2310E-05_r8, &
        0.2747E-05_r8,  0.3113E-05_r8,  0.3455E-05_r8,  0.3739E-05_r8,  0.3715E-05_r8, &
        0.3125E-05_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip=15,15)/ &
        0.99991709_r8,  0.99984968_r8,  0.99973857_r8,  0.99956548_r8,  0.99930853_r8, &
        0.99893898_r8,  0.99841601_r8,  0.99768001_r8,  0.99664098_r8,  0.99516898_r8, &
        0.99308002_r8,  0.99012297_r8,  0.98594999_r8,  0.98009998_r8,  0.97194999_r8, &
        0.96066999_r8,  0.94523001_r8,  0.92421001_r8,  0.89579999_r8,  0.85769999_r8, &
        0.80760002_r8,  0.74360001_r8,  0.66490000_r8,  0.57290000_r8,  0.47200000_r8, &
        0.36940002_r8,  0.27139997_r8,  0.18180001_r8,  0.10479999_r8,  0.04699999_r8, &
        0.01359999_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip=15,15)/ &
       -0.1378E-06_r8, -0.2542E-06_r8, -0.4461E-06_r8, -0.7333E-06_r8, -0.1125E-05_r8, &
       -0.1630E-05_r8, -0.2281E-05_r8, -0.3159E-05_r8, -0.4410E-05_r8, -0.6246E-05_r8, &
       -0.8933E-05_r8, -0.1280E-04_r8, -0.1826E-04_r8, -0.2589E-04_r8, -0.3639E-04_r8, &
       -0.5059E-04_r8, -0.6970E-04_r8, -0.9552E-04_r8, -0.1307E-03_r8, -0.1784E-03_r8, &
       -0.2422E-03_r8, -0.3237E-03_r8, -0.4203E-03_r8, -0.5227E-03_r8, -0.6184E-03_r8, &
       -0.6953E-03_r8, -0.7395E-03_r8, -0.7315E-03_r8, -0.6487E-03_r8, -0.4799E-03_r8, &
       -0.2625E-03_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip=15,15)/ &
        0.5891E-09_r8,  0.1074E-08_r8,  0.1885E-08_r8,  0.3167E-08_r8,  0.5051E-08_r8, &
        0.7631E-08_r8,  0.1092E-07_r8,  0.1500E-07_r8,  0.2032E-07_r8,  0.2769E-07_r8, &
        0.3810E-07_r8,  0.5279E-07_r8,  0.7361E-07_r8,  0.1032E-06_r8,  0.1450E-06_r8, &
        0.2026E-06_r8,  0.2798E-06_r8,  0.3832E-06_r8,  0.5242E-06_r8,  0.7159E-06_r8, &
        0.9706E-06_r8,  0.1299E-05_r8,  0.1701E-05_r8,  0.2159E-05_r8,  0.2612E-05_r8, &
        0.2998E-05_r8,  0.3341E-05_r8,  0.3661E-05_r8,  0.3775E-05_r8,  0.3393E-05_r8, &
        0.2384E-05_r8,  & 
       ! data ((h11(ip,iw),iw=1,31)_r8, ip=16,16)/ &
        0.99991363_r8,  0.99984020_r8,  0.99971467_r8,  0.99951237_r8,  0.99920303_r8, &
        0.99874902_r8,  0.99809903_r8,  0.99717999_r8,  0.99588197_r8,  0.99404502_r8, &
        0.99144298_r8,  0.98776001_r8,  0.98258001_r8,  0.97533000_r8,  0.96524000_r8, &
        0.95135999_r8,  0.93241000_r8,  0.90667999_r8,  0.87199998_r8,  0.82620001_r8, &
        0.76700002_r8,  0.69309998_r8,  0.60510004_r8,  0.50650001_r8,  0.40359998_r8, &
        0.30350000_r8,  0.21069998_r8,  0.12870002_r8,  0.06370002_r8,  0.02200001_r8, &
        0.00389999_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip=16,16)/ &
       -0.1383E-06_r8, -0.2577E-06_r8, -0.4608E-06_r8, -0.7793E-06_r8, -0.1237E-05_r8, &
       -0.1850E-05_r8, -0.2652E-05_r8, -0.3728E-05_r8, -0.5244E-05_r8, -0.7451E-05_r8, &
       -0.1067E-04_r8, -0.1532E-04_r8, -0.2193E-04_r8, -0.3119E-04_r8, -0.4395E-04_r8, &
       -0.6126E-04_r8, -0.8466E-04_r8, -0.1164E-03_r8, -0.1596E-03_r8, -0.2177E-03_r8, &
       -0.2933E-03_r8, -0.3855E-03_r8, -0.4874E-03_r8, -0.5870E-03_r8, -0.6718E-03_r8, &
       -0.7290E-03_r8, -0.7411E-03_r8, -0.6859E-03_r8, -0.5450E-03_r8, -0.3353E-03_r8, &
       -0.1363E-03_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip=16,16)/ &
        0.6217E-09_r8,  0.1165E-08_r8,  0.2116E-08_r8,  0.3685E-08_r8,  0.6101E-08_r8, &
        0.9523E-08_r8,  0.1400E-07_r8,  0.1959E-07_r8,  0.2668E-07_r8,  0.3629E-07_r8, &
        0.4982E-07_r8,  0.6876E-07_r8,  0.9523E-07_r8,  0.1321E-06_r8,  0.1825E-06_r8, &
        0.2505E-06_r8,  0.3420E-06_r8,  0.4677E-06_r8,  0.6416E-06_r8,  0.8760E-06_r8, &
        0.1183E-05_r8,  0.1565E-05_r8,  0.2010E-05_r8,  0.2472E-05_r8,  0.2882E-05_r8, &
        0.3229E-05_r8,  0.3564E-05_r8,  0.3777E-05_r8,  0.3589E-05_r8,  0.2786E-05_r8, &
        0.1487E-05_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip=17,17)/ &
        0.99991077_r8,  0.99983180_r8,  0.99969262_r8,  0.99945968_r8,  0.99909151_r8, &
        0.99853700_r8,  0.99773198_r8,  0.99658400_r8,  0.99496001_r8,  0.99266702_r8, &
        0.98943001_r8,  0.98484999_r8,  0.97842997_r8,  0.96945000_r8,  0.95703000_r8, &
        0.93998998_r8,  0.91676998_r8,  0.88540000_r8,  0.84350002_r8,  0.78890002_r8, &
        0.71990001_r8,  0.63639998_r8,  0.54060000_r8,  0.43820000_r8,  0.33639997_r8, &
        0.24080002_r8,  0.15460002_r8,  0.08310002_r8,  0.03310001_r8,  0.00770003_r8, &
        0.00050002_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip=17,17)/ &
       -0.1405E-06_r8, -0.2649E-06_r8, -0.4829E-06_r8, -0.8398E-06_r8, -0.1379E-05_r8, &
       -0.2132E-05_r8, -0.3138E-05_r8, -0.4487E-05_r8, -0.6353E-05_r8, -0.9026E-05_r8, &
       -0.1290E-04_r8, -0.1851E-04_r8, -0.2650E-04_r8, -0.3772E-04_r8, -0.5319E-04_r8, &
       -0.7431E-04_r8, -0.1031E-03_r8, -0.1422E-03_r8, -0.1951E-03_r8, -0.2648E-03_r8, &
       -0.3519E-03_r8, -0.4518E-03_r8, -0.5537E-03_r8, -0.6449E-03_r8, -0.7133E-03_r8, &
       -0.7432E-03_r8, -0.7133E-03_r8, -0.6018E-03_r8, -0.4092E-03_r8, -0.1951E-03_r8, &
       -0.5345E-04_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip=17,17)/  &
        0.6457E-09_r8,  0.1235E-08_r8,  0.2303E-08_r8,  0.4149E-08_r8,  0.7120E-08_r8,  &
        0.1152E-07_r8,  0.1749E-07_r8,  0.2508E-07_r8,  0.3462E-07_r8,  0.4718E-07_r8,  &
        0.6452E-07_r8,  0.8874E-07_r8,  0.1222E-06_r8,  0.1675E-06_r8,  0.2276E-06_r8,  &
        0.3076E-06_r8,  0.4174E-06_r8,  0.5714E-06_r8,  0.7837E-06_r8,  0.1067E-05_r8,  &
        0.1428E-05_r8,  0.1859E-05_r8,  0.2327E-05_r8,  0.2760E-05_r8,  0.3122E-05_r8,  &
        0.3458E-05_r8,  0.3739E-05_r8,  0.3715E-05_r8,  0.3126E-05_r8,  0.1942E-05_r8,  &
        0.6977E-06_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip=18,18)/  &
        0.99990851_r8,  0.99982500_r8,  0.99967349_r8,  0.99941093_r8,  0.99897999_r8,  &
        0.99831200_r8,  0.99732101_r8,  0.99589097_r8,  0.99386197_r8,  0.99099803_r8,  &
        0.98695999_r8,  0.98128998_r8,  0.97333002_r8,  0.96227002_r8,  0.94700998_r8,  &
        0.92614001_r8,  0.89779997_r8,  0.85969996_r8,  0.80949998_r8,  0.74540001_r8,  &
        0.66649997_r8,  0.57420003_r8,  0.47310001_r8,  0.37029999_r8,  0.27200001_r8,  &
        0.18220001_r8,  0.10500002_r8,  0.04710001_r8,  0.01359999_r8,  0.00169998_r8,  &
        0.00000000_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip=18,18)/  &
       -0.1431E-06_r8, -0.2731E-06_r8, -0.5072E-06_r8, -0.9057E-06_r8, -0.1537E-05_r8,  &
       -0.2460E-05_r8, -0.3733E-05_r8, -0.5449E-05_r8, -0.7786E-05_r8, -0.1106E-04_r8,  &
       -0.1574E-04_r8, -0.2249E-04_r8, -0.3212E-04_r8, -0.4564E-04_r8, -0.6438E-04_r8,  &
       -0.9019E-04_r8, -0.1256E-03_r8, -0.1737E-03_r8, -0.2378E-03_r8, -0.3196E-03_r8,  &
       -0.4163E-03_r8, -0.5191E-03_r8, -0.6154E-03_r8, -0.6931E-03_r8, -0.7384E-03_r8,  &
       -0.7313E-03_r8, -0.6492E-03_r8, -0.4805E-03_r8, -0.2629E-03_r8, -0.8897E-04_r8,  &
       -0.1432E-04_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip=18,18)/  &
        0.6607E-09_r8,  0.1282E-08_r8,  0.2441E-08_r8,  0.4522E-08_r8,  0.8027E-08_r8,  &
        0.1348E-07_r8,  0.2122E-07_r8,  0.3139E-07_r8,  0.4435E-07_r8,  0.6095E-07_r8,  &
        0.8319E-07_r8,  0.1139E-06_r8,  0.1557E-06_r8,  0.2107E-06_r8,  0.2819E-06_r8,  &
        0.3773E-06_r8,  0.5107E-06_r8,  0.6982E-06_r8,  0.9542E-06_r8,  0.1290E-05_r8,  &
        0.1703E-05_r8,  0.2170E-05_r8,  0.2628E-05_r8,  0.3013E-05_r8,  0.3352E-05_r8,  &
        0.3669E-05_r8,  0.3780E-05_r8,  0.3397E-05_r8,  0.2386E-05_r8,  0.1062E-05_r8,  &
        0.2216E-06_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip=19,19)/  &
        0.99990678_r8,  0.99981970_r8,  0.99965781_r8,  0.99936831_r8,  0.99887598_r8,  &
        0.99808502_r8,  0.99687898_r8,  0.99510998_r8,  0.99257898_r8,  0.98900002_r8,  &
        0.98398000_r8,  0.97693998_r8,  0.96711999_r8,  0.95353001_r8,  0.93484998_r8,  &
        0.90934002_r8,  0.87479997_r8,  0.82900000_r8,  0.76960003_r8,  0.69550002_r8,  &
        0.60720003_r8,  0.50819999_r8,  0.40490001_r8,  0.30440003_r8,  0.21130002_r8,  &
        0.12910002_r8,  0.06389999_r8,  0.02200001_r8,  0.00389999_r8,  0.00010002_r8,  &
        0.00000000_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip=19,19)/  &
       -0.1454E-06_r8, -0.2805E-06_r8, -0.5296E-06_r8, -0.9685E-06_r8, -0.1695E-05_r8,  &
       -0.2812E-05_r8, -0.4412E-05_r8, -0.6606E-05_r8, -0.9573E-05_r8, -0.1363E-04_r8,  &
       -0.1932E-04_r8, -0.2743E-04_r8, -0.3897E-04_r8, -0.5520E-04_r8, -0.7787E-04_r8,  &
       -0.1094E-03_r8, -0.1529E-03_r8, -0.2117E-03_r8, -0.2880E-03_r8, -0.3809E-03_r8,  &
       -0.4834E-03_r8, -0.5836E-03_r8, -0.6692E-03_r8, -0.7275E-03_r8, -0.7408E-03_r8,  &
       -0.6865E-03_r8, -0.5459E-03_r8, -0.3360E-03_r8, -0.1365E-03_r8, -0.2935E-04_r8,  &
       -0.2173E-05_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip=19,19)/  &
        0.6693E-09_r8,  0.1312E-08_r8,  0.2538E-08_r8,  0.4802E-08_r8,  0.8778E-08_r8,  &
        0.1528E-07_r8,  0.2501E-07_r8,  0.3836E-07_r8,  0.5578E-07_r8,  0.7806E-07_r8,  &
        0.1069E-06_r8,  0.1456E-06_r8,  0.1970E-06_r8,  0.2631E-06_r8,  0.3485E-06_r8,  &
        0.4642E-06_r8,  0.6268E-06_r8,  0.8526E-06_r8,  0.1157E-05_r8,  0.1545E-05_r8,  &
        0.2002E-05_r8,  0.2478E-05_r8,  0.2897E-05_r8,  0.3245E-05_r8,  0.3578E-05_r8,  &
        0.3789E-05_r8,  0.3598E-05_r8,  0.2792E-05_r8,  0.1489E-05_r8,  0.4160E-06_r8,  &
        0.3843E-07_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip=20,20)/  &
        0.99990559_r8,  0.99981570_r8,  0.99964547_r8,  0.99933308_r8,  0.99878299_r8,  &
        0.99786699_r8,  0.99642301_r8,  0.99425799_r8,  0.99111998_r8,  0.98667002_r8,  &
        0.98041999_r8,  0.97170001_r8,  0.95960999_r8,  0.94295001_r8,  0.92012000_r8,  &
        0.88900000_r8,  0.84740001_r8,  0.79280001_r8,  0.72360003_r8,  0.63960004_r8,  &
        0.54330003_r8,  0.44029999_r8,  0.33800000_r8,  0.24190003_r8,  0.15530002_r8,  &
        0.08350003_r8,  0.03320003_r8,  0.00770003_r8,  0.00050002_r8,  0.00000000_r8,  &
        0.00000000_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip=20,20)/  &
       -0.1472E-06_r8, -0.2866E-06_r8, -0.5485E-06_r8, -0.1024E-05_r8, -0.1842E-05_r8,  &
       -0.3160E-05_r8, -0.5136E-05_r8, -0.7922E-05_r8, -0.1171E-04_r8, -0.1682E-04_r8,  &
       -0.2381E-04_r8, -0.3355E-04_r8, -0.4729E-04_r8, -0.6673E-04_r8, -0.9417E-04_r8,  &
       -0.1327E-03_r8, -0.1858E-03_r8, -0.2564E-03_r8, -0.3449E-03_r8, -0.4463E-03_r8,  &
       -0.5495E-03_r8, -0.6420E-03_r8, -0.7116E-03_r8, -0.7427E-03_r8, -0.7139E-03_r8,  &
       -0.6031E-03_r8, -0.4104E-03_r8, -0.1957E-03_r8, -0.5358E-04_r8, -0.6176E-05_r8,  &
       -0.1347E-06_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip=20,20)/ &
        0.6750E-09_r8,  0.1332E-08_r8,  0.2602E-08_r8,  0.5003E-08_r8,  0.9367E-08_r8, &
        0.1684E-07_r8,  0.2863E-07_r8,  0.4566E-07_r8,  0.6865E-07_r8,  0.9861E-07_r8, &
        0.1368E-06_r8,  0.1856E-06_r8,  0.2479E-06_r8,  0.3274E-06_r8,  0.4315E-06_r8, &
        0.5739E-06_r8,  0.7710E-06_r8,  0.1040E-05_r8,  0.1394E-05_r8,  0.1829E-05_r8, &
        0.2309E-05_r8,  0.2759E-05_r8,  0.3131E-05_r8,  0.3472E-05_r8,  0.3755E-05_r8, &
        0.3730E-05_r8,  0.3138E-05_r8,  0.1948E-05_r8,  0.6994E-06_r8,  0.1022E-06_r8, &
        0.2459E-08_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip=21,21)/ &
        0.99990469_r8,  0.99981278_r8,  0.99963617_r8,  0.99930513_r8,  0.99870503_r8, &
        0.99766999_r8,  0.99597800_r8,  0.99336702_r8,  0.98951000_r8,  0.98399001_r8, &
        0.97622001_r8,  0.96543998_r8,  0.95059001_r8,  0.93019998_r8,  0.90235001_r8, &
        0.86470002_r8,  0.81480002_r8,  0.75059998_r8,  0.67129999_r8,  0.57840002_r8, &
        0.47659999_r8,  0.37279999_r8,  0.27389997_r8,  0.18339998_r8,  0.10570002_r8, &
        0.04740000_r8,  0.01370001_r8,  0.00169998_r8,  0.00000000_r8,  0.00000000_r8, &
        0.00000000_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip=21,21)/ &
       -0.1487E-06_r8, -0.2912E-06_r8, -0.5636E-06_r8, -0.1069E-05_r8, -0.1969E-05_r8, &
       -0.3483E-05_r8, -0.5858E-05_r8, -0.9334E-05_r8, -0.1416E-04_r8, -0.2067E-04_r8, &
       -0.2936E-04_r8, -0.4113E-04_r8, -0.5750E-04_r8, -0.8072E-04_r8, -0.1139E-03_r8, &
       -0.1606E-03_r8, -0.2246E-03_r8, -0.3076E-03_r8, -0.4067E-03_r8, -0.5121E-03_r8, &
       -0.6110E-03_r8, -0.6909E-03_r8, -0.7378E-03_r8, -0.7321E-03_r8, -0.6509E-03_r8, &
       -0.4825E-03_r8, -0.2641E-03_r8, -0.8936E-04_r8, -0.1436E-04_r8, -0.5966E-06_r8, &
        0.0000E+00_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip=21,21)/ &
        0.6777E-09_r8,  0.1344E-08_r8,  0.2643E-08_r8,  0.5138E-08_r8,  0.9798E-08_r8, &
        0.1809E-07_r8,  0.3185E-07_r8,  0.5285E-07_r8,  0.8249E-07_r8,  0.1222E-06_r8, &
        0.1730E-06_r8,  0.2351E-06_r8,  0.3111E-06_r8,  0.4078E-06_r8,  0.5366E-06_r8, &
        0.7117E-06_r8,  0.9495E-06_r8,  0.1266E-05_r8,  0.1667E-05_r8,  0.2132E-05_r8, &
        0.2600E-05_r8,  0.3001E-05_r8,  0.3354E-05_r8,  0.3679E-05_r8,  0.3796E-05_r8, &
        0.3414E-05_r8,  0.2399E-05_r8,  0.1067E-05_r8,  0.2222E-06_r8,  0.1075E-07_r8, &
        0.0000E+00_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip=22,22)/ &
        0.99990410_r8,  0.99981070_r8,  0.99962938_r8,  0.99928379_r8,  0.99864298_r8, &
        0.99750000_r8,  0.99556601_r8,  0.99247700_r8,  0.98780000_r8,  0.98100001_r8, &
        0.97140002_r8,  0.95810002_r8,  0.93984997_r8,  0.91491997_r8,  0.88110000_r8, &
        0.83570004_r8,  0.77670002_r8,  0.70239997_r8,  0.61350000_r8,  0.51349998_r8, &
        0.40910000_r8,  0.30750000_r8,  0.21340001_r8,  0.13040000_r8,  0.06449997_r8, &
        0.02219999_r8,  0.00389999_r8,  0.00010002_r8,  0.00000000_r8,  0.00000000_r8, &
        0.00000000_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip=22,22)/ &
       -0.1496E-06_r8, -0.2947E-06_r8, -0.5749E-06_r8, -0.1105E-05_r8, -0.2074E-05_r8, &
       -0.3763E-05_r8, -0.6531E-05_r8, -0.1076E-04_r8, -0.1682E-04_r8, -0.2509E-04_r8, &
       -0.3605E-04_r8, -0.5049E-04_r8, -0.7012E-04_r8, -0.9787E-04_r8, -0.1378E-03_r8, &
       -0.1939E-03_r8, -0.2695E-03_r8, -0.3641E-03_r8, -0.4703E-03_r8, -0.5750E-03_r8, &
       -0.6648E-03_r8, -0.7264E-03_r8, -0.7419E-03_r8, -0.6889E-03_r8, -0.5488E-03_r8, &
       -0.3382E-03_r8, -0.1375E-03_r8, -0.2951E-04_r8, -0.2174E-05_r8,  0.0000E+00_r8, &
        0.0000E+00_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip=22,22)/ &
        0.6798E-09_r8,  0.1350E-08_r8,  0.2667E-08_r8,  0.5226E-08_r8,  0.1010E-07_r8, &
        0.1903E-07_r8,  0.3455E-07_r8,  0.5951E-07_r8,  0.9658E-07_r8,  0.1479E-06_r8, &
        0.2146E-06_r8,  0.2951E-06_r8,  0.3903E-06_r8,  0.5101E-06_r8,  0.6693E-06_r8, &
        0.8830E-06_r8,  0.1168E-05_r8,  0.1532E-05_r8,  0.1968E-05_r8,  0.2435E-05_r8, &
        0.2859E-05_r8,  0.3222E-05_r8,  0.3572E-05_r8,  0.3797E-05_r8,  0.3615E-05_r8, &
        0.2811E-05_r8,  0.1500E-05_r8,  0.4185E-06_r8,  0.3850E-07_r8,  0.0000E+00_r8, &
        0.0000E+00_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip=23,23)/ &
        0.99990374_r8,  0.99980932_r8,  0.99962449_r8,  0.99926788_r8,  0.99859399_r8, &
        0.99736100_r8,  0.99520397_r8,  0.99163198_r8,  0.98606002_r8,  0.97779000_r8, &
        0.96600002_r8,  0.94963002_r8,  0.92720997_r8,  0.89670002_r8,  0.85580003_r8, &
        0.80170000_r8,  0.73269999_r8,  0.64840001_r8,  0.55110002_r8,  0.44669998_r8, &
        0.34280002_r8,  0.24529999_r8,  0.15750003_r8,  0.08469999_r8,  0.03369999_r8, &
        0.00779998_r8,  0.00050002_r8,  0.00000000_r8,  0.00000000_r8,  0.00000000_r8, &
        0.00000000_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip=23,23)/ &
       -0.1503E-06_r8, -0.2971E-06_r8, -0.5832E-06_r8, -0.1131E-05_r8, -0.2154E-05_r8, &
       -0.3992E-05_r8, -0.7122E-05_r8, -0.1211E-04_r8, -0.1954E-04_r8, -0.2995E-04_r8, &
       -0.4380E-04_r8, -0.6183E-04_r8, -0.8577E-04_r8, -0.1191E-03_r8, -0.1668E-03_r8, &
       -0.2333E-03_r8, -0.3203E-03_r8, -0.4237E-03_r8, -0.5324E-03_r8, -0.6318E-03_r8, &
       -0.7075E-03_r8, -0.7429E-03_r8, -0.7168E-03_r8, -0.6071E-03_r8, -0.4139E-03_r8, &
       -0.1976E-03_r8, -0.5410E-04_r8, -0.6215E-05_r8, -0.1343E-06_r8,  0.0000E+00_r8, &
        0.0000E+00_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip=23,23)/ &
        0.6809E-09_r8,  0.1356E-08_r8,  0.2683E-08_r8,  0.5287E-08_r8,  0.1030E-07_r8, &
        0.1971E-07_r8,  0.3665E-07_r8,  0.6528E-07_r8,  0.1100E-06_r8,  0.1744E-06_r8, &
        0.2599E-06_r8,  0.3650E-06_r8,  0.4887E-06_r8,  0.6398E-06_r8,  0.8358E-06_r8, &
        0.1095E-05_r8,  0.1429E-05_r8,  0.1836E-05_r8,  0.2286E-05_r8,  0.2716E-05_r8, &
        0.3088E-05_r8,  0.3444E-05_r8,  0.3748E-05_r8,  0.3740E-05_r8,  0.3157E-05_r8, &
        0.1966E-05_r8,  0.7064E-06_r8,  0.1030E-06_r8,  0.2456E-08_r8,  0.0000E+00_r8, &
        0.0000E+00_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip=24,24)/ &
        0.99990344_r8,  0.99980831_r8,  0.99962109_r8,  0.99925637_r8,  0.99855798_r8, &
        0.99725199_r8,  0.99489999_r8,  0.99087203_r8,  0.98436999_r8,  0.97447002_r8, &
        0.96012998_r8,  0.94006002_r8,  0.91254002_r8,  0.87540001_r8,  0.82609999_r8, &
        0.76240003_r8,  0.68299997_r8,  0.58930004_r8,  0.48589998_r8,  0.38020003_r8, &
        0.27920002_r8,  0.18699998_r8,  0.10769999_r8,  0.04830003_r8,  0.01400000_r8, &
        0.00169998_r8,  0.00000000_r8,  0.00000000_r8,  0.00000000_r8,  0.00000000_r8, &
        0.00000000_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip=24,24)/ &
       -0.1508E-06_r8, -0.2989E-06_r8, -0.5892E-06_r8, -0.1151E-05_r8, -0.2216E-05_r8, &
       -0.4175E-05_r8, -0.7619E-05_r8, -0.1333E-04_r8, -0.2217E-04_r8, -0.3497E-04_r8, &
       -0.5238E-04_r8, -0.7513E-04_r8, -0.1049E-03_r8, -0.1455E-03_r8, -0.2021E-03_r8, &
       -0.2790E-03_r8, -0.3757E-03_r8, -0.4839E-03_r8, -0.5902E-03_r8, -0.6794E-03_r8, &
       -0.7344E-03_r8, -0.7341E-03_r8, -0.6557E-03_r8, -0.4874E-03_r8, -0.2674E-03_r8, &
       -0.9059E-04_r8, -0.1455E-04_r8, -0.5986E-06_r8,  0.0000E+00_r8,  0.0000E+00_r8, &
        0.0000E+00_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip=24,24)/ &
        0.6812E-09_r8,  0.1356E-08_r8,  0.2693E-08_r8,  0.5328E-08_r8,  0.1045E-07_r8, &
        0.2021E-07_r8,  0.3826E-07_r8,  0.6994E-07_r8,  0.1218E-06_r8,  0.1997E-06_r8, &
        0.3069E-06_r8,  0.4428E-06_r8,  0.6064E-06_r8,  0.8015E-06_r8,  0.1043E-05_r8, &
        0.1351E-05_r8,  0.1733E-05_r8,  0.2168E-05_r8,  0.2598E-05_r8,  0.2968E-05_r8, &
        0.3316E-05_r8,  0.3662E-05_r8,  0.3801E-05_r8,  0.3433E-05_r8,  0.2422E-05_r8, &
        0.1081E-05_r8,  0.2256E-06_r8,  0.1082E-07_r8,  0.0000E+00_r8,  0.0000E+00_r8, &
        0.0000E+00_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip=25,25)/ &
        0.99990326_r8,  0.99980772_r8,  0.99961871_r8,  0.99924821_r8,  0.99853098_r8, &
        0.99716800_r8,  0.99465698_r8,  0.99022102_r8,  0.98281002_r8,  0.97118002_r8, &
        0.95393997_r8,  0.92948997_r8,  0.89579999_r8,  0.85070002_r8,  0.79189998_r8, &
        0.71759999_r8,  0.62800002_r8,  0.52639997_r8,  0.41970003_r8,  0.31559998_r8, &
        0.21899998_r8,  0.13370001_r8,  0.06610000_r8,  0.02280003_r8,  0.00400001_r8, &
        0.00010002_r8,  0.00000000_r8,  0.00000000_r8,  0.00000000_r8,  0.00000000_r8, &
        0.00000000_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip=25,25)/ &
       -0.1511E-06_r8, -0.3001E-06_r8, -0.5934E-06_r8, -0.1166E-05_r8, -0.2263E-05_r8, &
       -0.4319E-05_r8, -0.8028E-05_r8, -0.1438E-04_r8, -0.2460E-04_r8, -0.3991E-04_r8, &
       -0.6138E-04_r8, -0.9005E-04_r8, -0.1278E-03_r8, -0.1778E-03_r8, -0.2447E-03_r8, &
       -0.3313E-03_r8, -0.4342E-03_r8, -0.5424E-03_r8, -0.6416E-03_r8, -0.7146E-03_r8, &
       -0.7399E-03_r8, -0.6932E-03_r8, -0.5551E-03_r8, -0.3432E-03_r8, -0.1398E-03_r8, &
       -0.3010E-04_r8, -0.2229E-05_r8,  0.0000E+00_r8,  0.0000E+00_r8,  0.0000E+00_r8, &
        0.0000E+00_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip=25,25)/ &
        0.6815E-09_r8,  0.1358E-08_r8,  0.2698E-08_r8,  0.5355E-08_r8,  0.1054E-07_r8, &
        0.2056E-07_r8,  0.3942E-07_r8,  0.7349E-07_r8,  0.1315E-06_r8,  0.2226E-06_r8, &
        0.3537E-06_r8,  0.5266E-06_r8,  0.7407E-06_r8,  0.9958E-06_r8,  0.1296E-05_r8, &
        0.1657E-05_r8,  0.2077E-05_r8,  0.2512E-05_r8,  0.2893E-05_r8,  0.3216E-05_r8, &
        0.3562E-05_r8,  0.3811E-05_r8,  0.3644E-05_r8,  0.2841E-05_r8,  0.1524E-05_r8, &
        0.4276E-06_r8,  0.3960E-07_r8,  0.0000E+00_r8,  0.0000E+00_r8,  0.0000E+00_r8, &
        0.0000E+00_r8,  &
       ! data ((h11(ip,iw),iw=1,31)_r8, ip=26,26)/ &
        0.99990320_r8,  0.99980718_r8,  0.99961710_r8,  0.99924242_r8,  0.99851102_r8, &
        0.99710602_r8,  0.99446702_r8,  0.98969001_r8,  0.98144001_r8,  0.96805000_r8, &
        0.94762999_r8,  0.91812998_r8,  0.87730002_r8,  0.82290000_r8,  0.75319999_r8, &
        0.66789997_r8,  0.56879997_r8,  0.46160001_r8,  0.35450000_r8,  0.25370002_r8, &
        0.16280001_r8,  0.08740002_r8,  0.03479999_r8,  0.00809997_r8,  0.00059998_r8, &
        0.00000000_r8,  0.00000000_r8,  0.00000000_r8,  0.00000000_r8,  0.00000000_r8, &
        0.00000000_r8,  &
       ! data ((h12(ip,iw),iw=1,31)_r8, ip=26,26)/ &
       -0.1513E-06_r8, -0.3009E-06_r8, -0.5966E-06_r8, -0.1176E-05_r8, -0.2299E-05_r8, &
       -0.4430E-05_r8, -0.8352E-05_r8, -0.1526E-04_r8, -0.2674E-04_r8, -0.4454E-04_r8, &
       -0.7042E-04_r8, -0.1062E-03_r8, -0.1540E-03_r8, -0.2163E-03_r8, -0.2951E-03_r8, &
       -0.3899E-03_r8, -0.4948E-03_r8, -0.5983E-03_r8, -0.6846E-03_r8, -0.7332E-03_r8, &
       -0.7182E-03_r8, -0.6142E-03_r8, -0.4209E-03_r8, -0.2014E-03_r8, -0.5530E-04_r8, &
       -0.6418E-05_r8, -0.1439E-06_r8,  0.0000E+00_r8,  0.0000E+00_r8,  0.0000E+00_r8, &
        0.0000E+00_r8,  &
       ! data ((h13(ip,iw),iw=1,31)_r8, ip=26,26)/ &
        0.6817E-09_r8,  0.1359E-08_r8,  0.2702E-08_r8,  0.5374E-08_r8,  0.1061E-07_r8, &
        0.2079E-07_r8,  0.4022E-07_r8,  0.7610E-07_r8,  0.1392E-06_r8,  0.2428E-06_r8, &
        0.3992E-06_r8,  0.6149E-06_r8,  0.8893E-06_r8,  0.1220E-05_r8,  0.1599E-05_r8, &
        0.2015E-05_r8,  0.2453E-05_r8,  0.2853E-05_r8,  0.3173E-05_r8,  0.3488E-05_r8, &
        0.3792E-05_r8,  0.3800E-05_r8,  0.3210E-05_r8,  0.2002E-05_r8,  0.7234E-06_r8, &
        0.1068E-06_r8,  0.2646E-08_r8,  0.0000E+00_r8,  0.0000E+00_r8,  0.0000E+00_r8, &
        0.0000E+00_r8/),SHAPE=(/nx*3*nh/))

    it=0
    DO k=1,nx
       DO i=1,3
          DO j=1,nh
             it=it+1
              IF(i==1)THEN
                 !WRITE(*,'(a5,2e17.9) ' )'h11   ',data3(it),h11(k,j)
                 h11(k,j)=data3(it)
              END IF 
              IF(i==2) THEN
                 !WRITE(*,'(a5,2e17.9) ' )'h12   ',data3(it),h12(k,j)
                 h12(k,j)=data3(it)
              END IF 
              IF(i==3) THEN
                 !WRITE(*,'(a5,2e17.9) ' )'h13   ',data3(it),h13(k,j)
                 h13(k,j) =data3(it)
              END IF
          END DO
       END DO
    END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      data4(1:nx*3*nh)=RESHAPE(SOURCE=(/&
       ! data ((h21(ip,iw),iw=1,31), ip= 1, 1)/ &
        0.99999607_r8,  0.99999237_r8,  0.99998546_r8,  0.99997294_r8,  0.99995142_r8, &
        0.99991685_r8,  0.99986511_r8,  0.99979371_r8,  0.99970162_r8,  0.99958909_r8, &
        0.99945778_r8,  0.99931037_r8,  0.99914628_r8,  0.99895900_r8,  0.99873799_r8, &
        0.99846601_r8,  0.99813002_r8,  0.99771398_r8,  0.99719697_r8,  0.99655598_r8, &
        0.99575800_r8,  0.99475598_r8,  0.99348903_r8,  0.99186200_r8,  0.98973000_r8, &
        0.98688000_r8,  0.98303002_r8,  0.97777998_r8,  0.97059000_r8,  0.96077001_r8, &
        0.94742000_r8,  & 
       ! data ((h22(ip,iw),iw=1,31)_r8, ip= 1_r8, 1)/ &
       -0.5622E-07_r8, -0.1071E-06_r8, -0.1983E-06_r8, -0.3533E-06_r8, -0.5991E-06_r8, &
       -0.9592E-06_r8, -0.1444E-05_r8, -0.2049E-05_r8, -0.2764E-05_r8, -0.3577E-05_r8, &
       -0.4469E-05_r8, -0.5467E-05_r8, -0.6654E-05_r8, -0.8137E-05_r8, -0.1002E-04_r8, &
       -0.1237E-04_r8, -0.1528E-04_r8, -0.1884E-04_r8, -0.2310E-04_r8, -0.2809E-04_r8, &
       -0.3396E-04_r8, -0.4098E-04_r8, -0.4960E-04_r8, -0.6058E-04_r8, -0.7506E-04_r8, &
       -0.9451E-04_r8, -0.1207E-03_r8, -0.1558E-03_r8, -0.2026E-03_r8, -0.2648E-03_r8, &
       -0.3468E-03_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip= 1_r8, 1)/ &
       -0.2195E-09_r8, -0.4031E-09_r8, -0.7043E-09_r8, -0.1153E-08_r8, -0.1737E-08_r8, &
       -0.2395E-08_r8, -0.3020E-08_r8, -0.3549E-08_r8, -0.4034E-08_r8, -0.4421E-08_r8, &
       -0.4736E-08_r8, -0.5681E-08_r8, -0.8289E-08_r8, -0.1287E-07_r8, -0.1873E-07_r8, &
       -0.2523E-07_r8, -0.3223E-07_r8, -0.3902E-07_r8, -0.4409E-07_r8, -0.4699E-07_r8, &
       -0.4782E-07_r8, -0.4705E-07_r8, -0.4657E-07_r8, -0.4885E-07_r8, -0.5550E-07_r8, &
       -0.6619E-07_r8, -0.7656E-07_r8, -0.8027E-07_r8, -0.7261E-07_r8, -0.4983E-07_r8, &
       -0.1101E-07_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip= 2_r8, 2)/ &
        0.99999607_r8,  0.99999237_r8,  0.99998546_r8,  0.99997294_r8,  0.99995142_r8, &
        0.99991679_r8,  0.99986511_r8,  0.99979353_r8,  0.99970138_r8,  0.99958861_r8, &
        0.99945688_r8,  0.99930882_r8,  0.99914342_r8,  0.99895400_r8,  0.99872798_r8, &
        0.99844801_r8,  0.99809802_r8,  0.99765801_r8,  0.99710101_r8,  0.99639499_r8, &
        0.99549901_r8,  0.99435198_r8,  0.99287099_r8,  0.99093699_r8,  0.98837000_r8, &
        0.98491001_r8,  0.98019999_r8,  0.97373998_r8,  0.96490002_r8,  0.95283002_r8, &
        0.93649000_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip= 2_r8, 2)/  &
       -0.5622E-07_r8, -0.1071E-06_r8, -0.1983E-06_r8, -0.3534E-06_r8, -0.5992E-06_r8,  &
       -0.9594E-06_r8, -0.1445E-05_r8, -0.2050E-05_r8, -0.2766E-05_r8, -0.3580E-05_r8,  &
       -0.4476E-05_r8, -0.5479E-05_r8, -0.6677E-05_r8, -0.8179E-05_r8, -0.1009E-04_r8,  &
       -0.1251E-04_r8, -0.1553E-04_r8, -0.1928E-04_r8, -0.2384E-04_r8, -0.2930E-04_r8,  &
       -0.3588E-04_r8, -0.4393E-04_r8, -0.5403E-04_r8, -0.6714E-04_r8, -0.8458E-04_r8,  &
       -0.1082E-03_r8, -0.1400E-03_r8, -0.1829E-03_r8, -0.2401E-03_r8, -0.3157E-03_r8,  &
       -0.4147E-03_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip= 2_r8, 2)/  &
       -0.2195E-09_r8, -0.4032E-09_r8, -0.7046E-09_r8, -0.1153E-08_r8, -0.1738E-08_r8,  &
       -0.2395E-08_r8, -0.3021E-08_r8, -0.3550E-08_r8, -0.4035E-08_r8, -0.4423E-08_r8,  &
       -0.4740E-08_r8, -0.5692E-08_r8, -0.8314E-08_r8, -0.1292E-07_r8, -0.1882E-07_r8,  &
       -0.2536E-07_r8, -0.3242E-07_r8, -0.3927E-07_r8, -0.4449E-07_r8, -0.4767E-07_r8,  &
       -0.4889E-07_r8, -0.4857E-07_r8, -0.4860E-07_r8, -0.5132E-07_r8, -0.5847E-07_r8,  &
       -0.6968E-07_r8, -0.8037E-07_r8, -0.8400E-07_r8, -0.7521E-07_r8, -0.4830E-07_r8,  &
       -0.7562E-09_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip= 3_r8, 3)/  &
        0.99999607_r8,  0.99999237_r8,  0.99998546_r8,  0.99997294_r8,  0.99995142_r8,  &
        0.99991679_r8,  0.99986500_r8,  0.99979341_r8,  0.99970102_r8,  0.99958777_r8,  &
        0.99945557_r8,  0.99930632_r8,  0.99913889_r8,  0.99894601_r8,  0.99871302_r8,  &
        0.99842101_r8,  0.99805099_r8,  0.99757600_r8,  0.99696302_r8,  0.99617100_r8,  &
        0.99514598_r8,  0.99381000_r8,  0.99205798_r8,  0.98974001_r8,  0.98662001_r8,  &
        0.98238999_r8,  0.97659999_r8,  0.96866000_r8,  0.95776999_r8,  0.94296998_r8,  &
        0.92306000_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip= 3_r8, 3)/  &
       -0.5622E-07_r8, -0.1071E-06_r8, -0.1983E-06_r8, -0.3535E-06_r8, -0.5994E-06_r8,  &
       -0.9599E-06_r8, -0.1446E-05_r8, -0.2052E-05_r8, -0.2769E-05_r8, -0.3586E-05_r8,  &
       -0.4487E-05_r8, -0.5499E-05_r8, -0.6712E-05_r8, -0.8244E-05_r8, -0.1021E-04_r8,  &
       -0.1272E-04_r8, -0.1591E-04_r8, -0.1992E-04_r8, -0.2489E-04_r8, -0.3097E-04_r8,  &
       -0.3845E-04_r8, -0.4782E-04_r8, -0.5982E-04_r8, -0.7558E-04_r8, -0.9674E-04_r8,  &
       -0.1254E-03_r8, -0.1644E-03_r8, -0.2167E-03_r8, -0.2863E-03_r8, -0.3777E-03_r8,  &
       -0.4959E-03_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip= 3_r8, 3)/  &
       -0.2196E-09_r8, -0.4033E-09_r8, -0.7048E-09_r8, -0.1154E-08_r8, -0.1739E-08_r8,  &
       -0.2396E-08_r8, -0.3022E-08_r8, -0.3551E-08_r8, -0.4036E-08_r8, -0.4425E-08_r8,  &
       -0.4746E-08_r8, -0.5710E-08_r8, -0.8354E-08_r8, -0.1300E-07_r8, -0.1894E-07_r8,  &
       -0.2554E-07_r8, -0.3265E-07_r8, -0.3958E-07_r8, -0.4502E-07_r8, -0.4859E-07_r8,  &
       -0.5030E-07_r8, -0.5053E-07_r8, -0.5104E-07_r8, -0.5427E-07_r8, -0.6204E-07_r8,  &
       -0.7388E-07_r8, -0.8477E-07_r8, -0.8760E-07_r8, -0.7545E-07_r8, -0.4099E-07_r8,  &
        0.2046E-07_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip= 4_r8, 4)/  &
        0.99999607_r8,  0.99999237_r8,  0.99998546_r8,  0.99997294_r8,  0.99995142_r8,  &
        0.99991673_r8,  0.99986482_r8,  0.99979299_r8,  0.99970031_r8,  0.99958658_r8,  &
        0.99945343_r8,  0.99930239_r8,  0.99913180_r8,  0.99893302_r8,  0.99869001_r8,  &
        0.99838102_r8,  0.99798000_r8,  0.99745703_r8,  0.99676800_r8,  0.99586397_r8,  &
        0.99467200_r8,  0.99309403_r8,  0.99099600_r8,  0.98817998_r8,  0.98438001_r8,  &
        0.97918999_r8,  0.97206002_r8,  0.96227002_r8,  0.94888997_r8,  0.93080997_r8,  &
        0.90671003_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip= 4_r8, 4)/  &
       -0.5623E-07_r8, -0.1071E-06_r8, -0.1984E-06_r8, -0.3536E-06_r8, -0.5997E-06_r8,  &
       -0.9606E-06_r8, -0.1447E-05_r8, -0.2055E-05_r8, -0.2775E-05_r8, -0.3596E-05_r8,  &
       -0.4504E-05_r8, -0.5529E-05_r8, -0.6768E-05_r8, -0.8345E-05_r8, -0.1039E-04_r8,  &
       -0.1304E-04_r8, -0.1645E-04_r8, -0.2082E-04_r8, -0.2633E-04_r8, -0.3322E-04_r8,  &
       -0.4187E-04_r8, -0.5292E-04_r8, -0.6730E-04_r8, -0.8640E-04_r8, -0.1122E-03_r8,  &
       -0.1472E-03_r8, -0.1948E-03_r8, -0.2585E-03_r8, -0.3428E-03_r8, -0.4523E-03_r8,  &
       -0.5915E-03_r8,  &  
       ! data ((h23(ip,iw),iw=1,31)_r8, ip= 4_r8, 4)/  &
       -0.2196E-09_r8, -0.4034E-09_r8, -0.7050E-09_r8, -0.1154E-08_r8, -0.1740E-08_r8,  &
       -0.2398E-08_r8, -0.3024E-08_r8, -0.3552E-08_r8, -0.4037E-08_r8, -0.4428E-08_r8,  &
       -0.4756E-08_r8, -0.5741E-08_r8, -0.8418E-08_r8, -0.1310E-07_r8, -0.1910E-07_r8,  &
       -0.2575E-07_r8, -0.3293E-07_r8, -0.3998E-07_r8, -0.4572E-07_r8, -0.4980E-07_r8,  &
       -0.5211E-07_r8, -0.5287E-07_r8, -0.5390E-07_r8, -0.5782E-07_r8, -0.6650E-07_r8,  &
       -0.7892E-07_r8, -0.8940E-07_r8, -0.8980E-07_r8, -0.7119E-07_r8, -0.2452E-07_r8,  &
        0.5823E-07_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip= 5_r8, 5)/  &
        0.99999607_r8,  0.99999237_r8,  0.99998546_r8,  0.99997294_r8,  0.99995136_r8,  &
        0.99991661_r8,  0.99986458_r8,  0.99979252_r8,  0.99969929_r8,  0.99958479_r8,  &
        0.99945003_r8,  0.99929619_r8,  0.99912071_r8,  0.99891400_r8,  0.99865502_r8,  &
        0.99831998_r8,  0.99787700_r8,  0.99728799_r8,  0.99650002_r8,  0.99544799_r8,  &
        0.99404198_r8,  0.99215603_r8,  0.98961997_r8,  0.98619002_r8,  0.98153001_r8,  &
        0.97513002_r8,  0.96634001_r8,  0.95428002_r8,  0.93791002_r8,  0.91593999_r8,  &
        0.88700002_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip= 5_r8, 5)/  &
       -0.5623E-07_r8, -0.1071E-06_r8, -0.1985E-06_r8, -0.3538E-06_r8, -0.6002E-06_r8,  &
       -0.9618E-06_r8, -0.1450E-05_r8, -0.2059E-05_r8, -0.2783E-05_r8, -0.3611E-05_r8,  &
       -0.4531E-05_r8, -0.5577E-05_r8, -0.6855E-05_r8, -0.8499E-05_r8, -0.1066E-04_r8,  &
       -0.1351E-04_r8, -0.1723E-04_r8, -0.2207E-04_r8, -0.2829E-04_r8, -0.3621E-04_r8,  &
       -0.4636E-04_r8, -0.5954E-04_r8, -0.7690E-04_r8, -0.1002E-03_r8, -0.1317E-03_r8,  &
       -0.1746E-03_r8, -0.2326E-03_r8, -0.3099E-03_r8, -0.4111E-03_r8, -0.5407E-03_r8,  &
       -0.7020E-03_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip= 5_r8, 5)/  &
       -0.2197E-09_r8, -0.4037E-09_r8, -0.7054E-09_r8, -0.1155E-08_r8, -0.1741E-08_r8,  &
       -0.2401E-08_r8, -0.3027E-08_r8, -0.3553E-08_r8, -0.4039E-08_r8, -0.4431E-08_r8,  &
       -0.4775E-08_r8, -0.5784E-08_r8, -0.8506E-08_r8, -0.1326E-07_r8, -0.1931E-07_r8,  &
       -0.2600E-07_r8, -0.3324E-07_r8, -0.4048E-07_r8, -0.4666E-07_r8, -0.5137E-07_r8,  &
       -0.5428E-07_r8, -0.5558E-07_r8, -0.5730E-07_r8, -0.6228E-07_r8, -0.7197E-07_r8,  &
       -0.8455E-07_r8, -0.9347E-07_r8, -0.8867E-07_r8, -0.5945E-07_r8,  0.5512E-08_r8,  &
        0.1209E-06_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip= 6_r8, 6)/ &
        0.99999607_r8,  0.99999237_r8,  0.99998546_r8,  0.99997288_r8,  0.99995130_r8, &
        0.99991649_r8,  0.99986428_r8,  0.99979180_r8,  0.99969780_r8,  0.99958187_r8, &
        0.99944460_r8,  0.99928659_r8,  0.99910372_r8,  0.99888301_r8,  0.99860299_r8, &
        0.99822998_r8,  0.99773002_r8,  0.99705303_r8,  0.99613500_r8,  0.99489301_r8, &
        0.99321300_r8,  0.99093801_r8,  0.98785001_r8,  0.98365998_r8,  0.97790998_r8, &
        0.97000998_r8,  0.95916998_r8,  0.94437003_r8,  0.92440999_r8,  0.89789999_r8, &
        0.86360002_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip= 6_r8, 6)/ &
       -0.5624E-07_r8, -0.1072E-06_r8, -0.1986E-06_r8, -0.3541E-06_r8, -0.6010E-06_r8, &
       -0.9636E-06_r8, -0.1453E-05_r8, -0.2067E-05_r8, -0.2796E-05_r8, -0.3634E-05_r8, &
       -0.4572E-05_r8, -0.5652E-05_r8, -0.6987E-05_r8, -0.8733E-05_r8, -0.1107E-04_r8, &
       -0.1418E-04_r8, -0.1832E-04_r8, -0.2378E-04_r8, -0.3092E-04_r8, -0.4017E-04_r8, &
       -0.5221E-04_r8, -0.6806E-04_r8, -0.8916E-04_r8, -0.1176E-03_r8, -0.1562E-03_r8, &
       -0.2087E-03_r8, -0.2793E-03_r8, -0.3724E-03_r8, -0.4928E-03_r8, -0.6440E-03_r8, &
       -0.8270E-03_r8,  & 
       ! data ((h23(ip,iw),iw=1,31)_r8, ip= 6_r8, 6)/ &
       -0.2198E-09_r8, -0.4040E-09_r8, -0.7061E-09_r8, -0.1156E-08_r8, -0.1744E-08_r8, &
       -0.2405E-08_r8, -0.3032E-08_r8, -0.3556E-08_r8, -0.4040E-08_r8, -0.4444E-08_r8, &
       -0.4800E-08_r8, -0.5848E-08_r8, -0.8640E-08_r8, -0.1346E-07_r8, -0.1957E-07_r8, &
       -0.2627E-07_r8, -0.3357E-07_r8, -0.4114E-07_r8, -0.4793E-07_r8, -0.5330E-07_r8, &
       -0.5676E-07_r8, -0.5873E-07_r8, -0.6152E-07_r8, -0.6783E-07_r8, -0.7834E-07_r8, &
       -0.9023E-07_r8, -0.9530E-07_r8, -0.8162E-07_r8, -0.3634E-07_r8,  0.5638E-07_r8, &
        0.2189E-06_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip= 7_r8, 7)/ &
        0.99999607_r8,  0.99999237_r8,  0.99998546_r8,  0.99997288_r8,  0.99995124_r8, &
        0.99991626_r8,  0.99986368_r8,  0.99979049_r8,  0.99969530_r8,  0.99957728_r8, &
        0.99943632_r8,  0.99927181_r8,  0.99907762_r8,  0.99883801_r8,  0.99852502_r8, &
        0.99810201_r8,  0.99752498_r8,  0.99673301_r8,  0.99564600_r8,  0.99416101_r8, &
        0.99213398_r8,  0.98936999_r8,  0.98559999_r8,  0.98043001_r8,  0.97333997_r8, &
        0.96359003_r8,  0.95025003_r8,  0.93216002_r8,  0.90798998_r8,  0.87639999_r8, &
        0.83609998_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip= 7_r8, 7)/ &
       -0.5626E-07_r8, -0.1072E-06_r8, -0.1987E-06_r8, -0.3545E-06_r8, -0.6022E-06_r8, &
       -0.9665E-06_r8, -0.1460E-05_r8, -0.2078E-05_r8, -0.2817E-05_r8, -0.3671E-05_r8, &
       -0.4637E-05_r8, -0.5767E-05_r8, -0.7188E-05_r8, -0.9080E-05_r8, -0.1165E-04_r8, &
       -0.1513E-04_r8, -0.1981E-04_r8, -0.2609E-04_r8, -0.3441E-04_r8, -0.4534E-04_r8, &
       -0.5978E-04_r8, -0.7897E-04_r8, -0.1047E-03_r8, -0.1396E-03_r8, -0.1870E-03_r8, &
       -0.2510E-03_r8, -0.3363E-03_r8, -0.4475E-03_r8, -0.5888E-03_r8, -0.7621E-03_r8, &
       -0.9647E-03_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip= 7_r8, 7)/  &
       -0.2200E-09_r8, -0.4045E-09_r8, -0.7071E-09_r8, -0.1159E-08_r8, -0.1748E-08_r8,  &
       -0.2411E-08_r8, -0.3040E-08_r8, -0.3561E-08_r8, -0.4046E-08_r8, -0.4455E-08_r8,  &
       -0.4839E-08_r8, -0.5941E-08_r8, -0.8815E-08_r8, -0.1371E-07_r8, -0.1983E-07_r8,  &
       -0.2652E-07_r8, -0.3400E-07_r8, -0.4207E-07_r8, -0.4955E-07_r8, -0.5554E-07_r8,  &
       -0.5966E-07_r8, -0.6261E-07_r8, -0.6688E-07_r8, -0.7454E-07_r8, -0.8521E-07_r8,  &
       -0.9470E-07_r8, -0.9275E-07_r8, -0.6525E-07_r8,  0.3686E-08_r8,  0.1371E-06_r8,  &
        0.3623E-06_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip= 8_r8, 8)/  &
        0.99999607_r8,  0.99999237_r8,  0.99998540_r8,  0.99997282_r8,  0.99995112_r8,  &
        0.99991590_r8,  0.99986279_r8,  0.99978858_r8,  0.99969149_r8,  0.99957019_r8,  &
        0.99942350_r8,  0.99924922_r8,  0.99903822_r8,  0.99877101_r8,  0.99841398_r8,  &
        0.99792302_r8,  0.99724299_r8,  0.99630302_r8,  0.99500000_r8,  0.99320602_r8,  &
        0.99074000_r8,  0.98736000_r8,  0.98272002_r8,  0.97635001_r8,  0.96758002_r8,  &
        0.95555997_r8,  0.93919998_r8,  0.91722000_r8,  0.88819999_r8,  0.85089999_r8,  &
        0.80439997_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip= 8_r8, 8)/  &
       -0.5628E-07_r8, -0.1073E-06_r8, -0.1990E-06_r8, -0.3553E-06_r8, -0.6042E-06_r8,  &
       -0.9710E-06_r8, -0.1469E-05_r8, -0.2096E-05_r8, -0.2849E-05_r8, -0.3728E-05_r8,  &
       -0.4738E-05_r8, -0.5942E-05_r8, -0.7490E-05_r8, -0.9586E-05_r8, -0.1247E-04_r8,  &
       -0.1644E-04_r8, -0.2184E-04_r8, -0.2916E-04_r8, -0.3898E-04_r8, -0.5205E-04_r8,  &
       -0.6948E-04_r8, -0.9285E-04_r8, -0.1244E-03_r8, -0.1672E-03_r8, -0.2251E-03_r8,  &
       -0.3028E-03_r8, -0.4051E-03_r8, -0.5365E-03_r8, -0.6998E-03_r8, -0.8940E-03_r8,  &
       -0.1112E-02_r8,  &  
       ! data ((h23(ip,iw),iw=1,31)_r8, ip= 8_r8, 8)/  &
       -0.2204E-09_r8, -0.4052E-09_r8, -0.7088E-09_r8, -0.1162E-08_r8, -0.1755E-08_r8,  &
       -0.2422E-08_r8, -0.3053E-08_r8, -0.3572E-08_r8, -0.4052E-08_r8, -0.4474E-08_r8,  &
       -0.4898E-08_r8, -0.6082E-08_r8, -0.9046E-08_r8, -0.1400E-07_r8, -0.2009E-07_r8,  &
       -0.2683E-07_r8, -0.3463E-07_r8, -0.4334E-07_r8, -0.5153E-07_r8, -0.5811E-07_r8,  &
       -0.6305E-07_r8, -0.6749E-07_r8, -0.7346E-07_r8, -0.8208E-07_r8, -0.9173E-07_r8,  &
       -0.9603E-07_r8, -0.8264E-07_r8, -0.3505E-07_r8,  0.6878E-07_r8,  0.2586E-06_r8,  &
        0.5530E-06_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip= 9_r8, 9)/  &
        0.99999607_r8,  0.99999237_r8,  0.99998540_r8,  0.99997276_r8,  0.99995089_r8,  &
        0.99991536_r8,  0.99986148_r8,  0.99978572_r8,  0.99968570_r8,  0.99955928_r8,  &
        0.99940401_r8,  0.99921501_r8,  0.99897999_r8,  0.99867398_r8,  0.99825603_r8,  &
        0.99767601_r8,  0.99686497_r8,  0.99573302_r8,  0.99415499_r8,  0.99196899_r8,  &
        0.98895001_r8,  0.98479998_r8,  0.97907001_r8,  0.97119999_r8,  0.96038002_r8,  &
        0.94559997_r8,  0.92565000_r8,  0.89910001_r8,  0.86470002_r8,  0.82130003_r8,  &
        0.76830000_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip= 9_r8, 9)/  & 
       -0.5630E-07_r8, -0.1074E-06_r8, -0.1994E-06_r8, -0.3564E-06_r8, -0.6072E-06_r8,  & 
       -0.9779E-06_r8, -0.1484E-05_r8, -0.2124E-05_r8, -0.2900E-05_r8, -0.3817E-05_r8,  & 
       -0.4891E-05_r8, -0.6205E-05_r8, -0.7931E-05_r8, -0.1031E-04_r8, -0.1362E-04_r8,  & 
       -0.1821E-04_r8, -0.2454E-04_r8, -0.3320E-04_r8, -0.4493E-04_r8, -0.6068E-04_r8,  & 
       -0.8186E-04_r8, -0.1104E-03_r8, -0.1491E-03_r8, -0.2016E-03_r8, -0.2722E-03_r8,  & 
       -0.3658E-03_r8, -0.4873E-03_r8, -0.6404E-03_r8, -0.8253E-03_r8, -0.1037E-02_r8,  & 
       -0.1267E-02_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip= 9_r8, 9)/  & 
       -0.2207E-09_r8, -0.4061E-09_r8, -0.7117E-09_r8, -0.1169E-08_r8, -0.1767E-08_r8,  & 
       -0.2439E-08_r8, -0.3074E-08_r8, -0.3588E-08_r8, -0.4062E-08_r8, -0.4510E-08_r8,  & 
       -0.4983E-08_r8, -0.6261E-08_r8, -0.9324E-08_r8, -0.1430E-07_r8, -0.2036E-07_r8,  & 
       -0.2725E-07_r8, -0.3561E-07_r8, -0.4505E-07_r8, -0.5384E-07_r8, -0.6111E-07_r8,  & 
       -0.6731E-07_r8, -0.7355E-07_r8, -0.8112E-07_r8, -0.8978E-07_r8, -0.9616E-07_r8,  & 
       -0.9157E-07_r8, -0.6114E-07_r8,  0.1622E-07_r8,  0.1694E-06_r8,  0.4277E-06_r8,  & 
        0.7751E-06_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip=10,10)/  & 
        0.99999607_r8,  0.99999237_r8,  0.99998540_r8,  0.99997264_r8,  0.99995059_r8,  & 
        0.99991453_r8,  0.99985939_r8,  0.99978119_r8,  0.99967682_r8,  0.99954277_r8,  & 
        0.99937469_r8,  0.99916458_r8,  0.99889499_r8,  0.99853599_r8,  0.99804002_r8,  & 
        0.99734300_r8,  0.99636298_r8,  0.99498600_r8,  0.99305803_r8,  0.99037802_r8,  & 
        0.98667002_r8,  0.98154002_r8,  0.97447002_r8,  0.96473998_r8,  0.95141000_r8,  & 
        0.93333000_r8,  0.90916002_r8,  0.87750000_r8,  0.83710003_r8,  0.78729999_r8,  & 
        0.72790003_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip=10,10)/  & 
       -0.5636E-07_r8, -0.1076E-06_r8, -0.2000E-06_r8, -0.3582E-06_r8, -0.6119E-06_r8,  & 
       -0.9888E-06_r8, -0.1507E-05_r8, -0.2168E-05_r8, -0.2978E-05_r8, -0.3952E-05_r8,  & 
       -0.5122E-05_r8, -0.6592E-05_r8, -0.8565E-05_r8, -0.1132E-04_r8, -0.1518E-04_r8,  & 
       -0.2060E-04_r8, -0.2811E-04_r8, -0.3848E-04_r8, -0.5261E-04_r8, -0.7173E-04_r8,  & 
       -0.9758E-04_r8, -0.1326E-03_r8, -0.1801E-03_r8, -0.2442E-03_r8, -0.3296E-03_r8,  & 
       -0.4415E-03_r8, -0.5840E-03_r8, -0.7591E-03_r8, -0.9636E-03_r8, -0.1189E-02_r8,  &
       -0.1427E-02_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip=10,10)/  &
       -0.2214E-09_r8, -0.4080E-09_r8, -0.7156E-09_r8, -0.1178E-08_r8, -0.1784E-08_r8,  &
       -0.2466E-08_r8, -0.3105E-08_r8, -0.3617E-08_r8, -0.4087E-08_r8, -0.4563E-08_r8,  &
       -0.5110E-08_r8, -0.6492E-08_r8, -0.9643E-08_r8, -0.1461E-07_r8, -0.2069E-07_r8,  &
       -0.2796E-07_r8, -0.3702E-07_r8, -0.4717E-07_r8, -0.5662E-07_r8, -0.6484E-07_r8,  &
       -0.7271E-07_r8, -0.8079E-07_r8, -0.8928E-07_r8, -0.9634E-07_r8, -0.9625E-07_r8,  &
       -0.7776E-07_r8, -0.2242E-07_r8,  0.9745E-07_r8,  0.3152E-06_r8,  0.6388E-06_r8,  &
        0.9992E-06_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip=11,11)/  &
        0.99999607_r8,  0.99999237_r8,  0.99998534_r8,  0.99997252_r8,  0.99995011_r8,  &
        0.99991328_r8,  0.99985629_r8,  0.99977452_r8,  0.99966347_r8,  0.99951839_r8,  &
        0.99933177_r8,  0.99909180_r8,  0.99877602_r8,  0.99834698_r8,  0.99774700_r8,  &
        0.99690098_r8,  0.99570400_r8,  0.99401599_r8,  0.99164802_r8,  0.98834997_r8,  &
        0.98377001_r8,  0.97742999_r8,  0.96868002_r8,  0.95666999_r8,  0.94032001_r8,  &
        0.91833997_r8,  0.88929999_r8,  0.85189998_r8,  0.80530000_r8,  0.74909997_r8,  &
        0.68299997_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip=11,11)/  &
       -0.5645E-07_r8, -0.1079E-06_r8, -0.2009E-06_r8, -0.3610E-06_r8, -0.6190E-06_r8,  &
       -0.1005E-05_r8, -0.1541E-05_r8, -0.2235E-05_r8, -0.3096E-05_r8, -0.4155E-05_r8,  &
       -0.5463E-05_r8, -0.7150E-05_r8, -0.9453E-05_r8, -0.1269E-04_r8, -0.1728E-04_r8,  &
       -0.2375E-04_r8, -0.3279E-04_r8, -0.4531E-04_r8, -0.6246E-04_r8, -0.8579E-04_r8,  &
       -0.1175E-03_r8, -0.1604E-03_r8, -0.2186E-03_r8, -0.2964E-03_r8, -0.3990E-03_r8,  &
       -0.5311E-03_r8, -0.6957E-03_r8, -0.8916E-03_r8, -0.1112E-02_r8, -0.1346E-02_r8,  &
       -0.1590E-02_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip=11,11)/  &
       -0.2225E-09_r8, -0.4104E-09_r8, -0.7217E-09_r8, -0.1192E-08_r8, -0.1811E-08_r8,  &
       -0.2509E-08_r8, -0.3155E-08_r8, -0.3668E-08_r8, -0.4138E-08_r8, -0.4650E-08_r8,  &
       -0.5296E-08_r8, -0.6785E-08_r8, -0.9991E-08_r8, -0.1494E-07_r8, -0.2122E-07_r8,  &
       -0.2911E-07_r8, -0.3895E-07_r8, -0.4979E-07_r8, -0.6002E-07_r8, -0.6964E-07_r8,  &
       -0.7935E-07_r8, -0.8887E-07_r8, -0.9699E-07_r8, -0.9967E-07_r8, -0.8883E-07_r8,  &
       -0.4988E-07_r8,  0.4156E-07_r8,  0.2197E-06_r8,  0.5081E-06_r8,  0.8667E-06_r8,  &
        0.1212E-05_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip=12,12)/  &
        0.99999607_r8,  0.99999237_r8,  0.99998528_r8,  0.99997234_r8,  0.99994951_r8,  &
        0.99991143_r8,  0.99985188_r8,  0.99976480_r8,  0.99964428_r8,  0.99948311_r8,  &
        0.99927050_r8,  0.99899000_r8,  0.99861199_r8,  0.99809098_r8,  0.99735999_r8,  &
        0.99632198_r8,  0.99484903_r8,  0.99276900_r8,  0.98984998_r8,  0.98576999_r8,  &
        0.98009998_r8,  0.97224998_r8,  0.96144003_r8,  0.94667000_r8,  0.92672002_r8,  &
        0.90020001_r8,  0.86570001_r8,  0.82220000_r8,  0.76919997_r8,  0.70640004_r8,  &
        0.63330001_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip=12,12)/  &
       -0.5658E-07_r8, -0.1083E-06_r8, -0.2023E-06_r8, -0.3650E-06_r8, -0.6295E-06_r8,  &
       -0.1030E-05_r8, -0.1593E-05_r8, -0.2334E-05_r8, -0.3273E-05_r8, -0.4455E-05_r8,  &
       -0.5955E-05_r8, -0.7935E-05_r8, -0.1067E-04_r8, -0.1455E-04_r8, -0.2007E-04_r8,  &
       -0.2788E-04_r8, -0.3885E-04_r8, -0.5409E-04_r8, -0.7503E-04_r8, -0.1036E-03_r8,  &
       -0.1425E-03_r8, -0.1951E-03_r8, -0.2660E-03_r8, -0.3598E-03_r8, -0.4817E-03_r8,  &
       -0.6355E-03_r8, -0.8218E-03_r8, -0.1035E-02_r8, -0.1267E-02_r8, -0.1508E-02_r8,  &
       -0.1755E-02_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip=12,12)/  &
       -0.2241E-09_r8, -0.4142E-09_r8, -0.7312E-09_r8, -0.1214E-08_r8, -0.1854E-08_r8,  &
       -0.2578E-08_r8, -0.3250E-08_r8, -0.3765E-08_r8, -0.4238E-08_r8, -0.4809E-08_r8,  &
       -0.5553E-08_r8, -0.7132E-08_r8, -0.1035E-07_r8, -0.1538E-07_r8, -0.2211E-07_r8,  &
       -0.3079E-07_r8, -0.4142E-07_r8, -0.5303E-07_r8, -0.6437E-07_r8, -0.7566E-07_r8,  &
       -0.8703E-07_r8, -0.9700E-07_r8, -0.1025E-06_r8, -0.9718E-07_r8, -0.6973E-07_r8,  &
       -0.5265E-09_r8,  0.1413E-06_r8,  0.3895E-06_r8,  0.7321E-06_r8,  0.1085E-05_r8,  &
        0.1449E-05_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip=13,13)/  &
        0.99999607_r8,  0.99999231_r8,  0.99998522_r8,  0.99997205_r8,  0.99994856_r8,  &
        0.99990892_r8,  0.99984580_r8,  0.99975121_r8,  0.99961728_r8,  0.99943388_r8,  &
        0.99918568_r8,  0.99884999_r8,  0.99839199_r8,  0.99775398_r8,  0.99685299_r8,  &
        0.99557197_r8,  0.99375200_r8,  0.99117702_r8,  0.98755997_r8,  0.98250002_r8,  &
        0.97548002_r8,  0.96575999_r8,  0.95244002_r8,  0.93436003_r8,  0.91018999_r8,  &
        0.87849998_r8,  0.83810002_r8,  0.78820002_r8,  0.72870004_r8,  0.65910000_r8,  &
        0.57850003_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip=13,13)/  &
       -0.5677E-07_r8, -0.1090E-06_r8, -0.2043E-06_r8, -0.3709E-06_r8, -0.6448E-06_r8,  &
       -0.1066E-05_r8, -0.1669E-05_r8, -0.2480E-05_r8, -0.3532E-05_r8, -0.4887E-05_r8,  &
       -0.6650E-05_r8, -0.9017E-05_r8, -0.1232E-04_r8, -0.1701E-04_r8, -0.2372E-04_r8,  &
       -0.3325E-04_r8, -0.4664E-04_r8, -0.6528E-04_r8, -0.9095E-04_r8, -0.1260E-03_r8,  &
       -0.1737E-03_r8, -0.2381E-03_r8, -0.3238E-03_r8, -0.4359E-03_r8, -0.5789E-03_r8,  &
       -0.7549E-03_r8, -0.9607E-03_r8, -0.1188E-02_r8, -0.1427E-02_r8, -0.1674E-02_r8,  &
       -0.1914E-02_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip=13,13)/  &
       -0.2262E-09_r8, -0.4191E-09_r8, -0.7441E-09_r8, -0.1244E-08_r8, -0.1916E-08_r8,  &
       -0.2687E-08_r8, -0.3407E-08_r8, -0.3947E-08_r8, -0.4432E-08_r8, -0.5059E-08_r8,  &
       -0.5896E-08_r8, -0.7538E-08_r8, -0.1079E-07_r8, -0.1606E-07_r8, -0.2346E-07_r8,  &
       -0.3305E-07_r8, -0.4459E-07_r8, -0.5720E-07_r8, -0.6999E-07_r8, -0.8300E-07_r8,  &
       -0.9536E-07_r8, -0.1039E-06_r8, -0.1036E-06_r8, -0.8526E-07_r8, -0.3316E-07_r8,  &
        0.7909E-07_r8,  0.2865E-06_r8,  0.6013E-06_r8,  0.9580E-06_r8,  0.1303E-05_r8,  &
        0.1792E-05_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip=14,14)/  &
        0.99999607_r8,  0.99999231_r8,  0.99998510_r8,  0.99997163_r8,  0.99994737_r8,  &
        0.99990571_r8,  0.99983770_r8,  0.99973333_r8,  0.99958128_r8,  0.99936771_r8,  &
        0.99907219_r8,  0.99866599_r8,  0.99810302_r8,  0.99731499_r8,  0.99620003_r8,  &
        0.99461198_r8,  0.99235398_r8,  0.98916000_r8,  0.98466998_r8,  0.97839999_r8,  &
        0.96969002_r8,  0.95769000_r8,  0.94133997_r8,  0.91935998_r8,  0.89029998_r8,  &
        0.85290003_r8,  0.80620003_r8,  0.74989998_r8,  0.68369997_r8,  0.60679996_r8,  &
        0.51899999_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip=14,14)/  &
       -0.5703E-07_r8, -0.1098E-06_r8, -0.2071E-06_r8, -0.3788E-06_r8, -0.6657E-06_r8,  &
       -0.1116E-05_r8, -0.1776E-05_r8, -0.2687E-05_r8, -0.3898E-05_r8, -0.5493E-05_r8,  &
       -0.7607E-05_r8, -0.1048E-04_r8, -0.1450E-04_r8, -0.2024E-04_r8, -0.2845E-04_r8,  &
       -0.4014E-04_r8, -0.5658E-04_r8, -0.7947E-04_r8, -0.1110E-03_r8, -0.1541E-03_r8,  &
       -0.2125E-03_r8, -0.2907E-03_r8, -0.3936E-03_r8, -0.5261E-03_r8, -0.6912E-03_r8,  &
       -0.8880E-03_r8, -0.1109E-02_r8, -0.1346E-02_r8, -0.1591E-02_r8, -0.1837E-02_r8,  &
       -0.2054E-02_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip=14,14)/  &
       -0.2288E-09_r8, -0.4265E-09_r8, -0.7627E-09_r8, -0.1289E-08_r8, -0.2011E-08_r8,  &
       -0.2861E-08_r8, -0.3673E-08_r8, -0.4288E-08_r8, -0.4812E-08_r8, -0.5475E-08_r8,  &
       -0.6365E-08_r8, -0.8052E-08_r8, -0.1142E-07_r8, -0.1711E-07_r8, -0.2533E-07_r8,  &
       -0.3597E-07_r8, -0.4862E-07_r8, -0.6259E-07_r8, -0.7708E-07_r8, -0.9150E-07_r8,  &
       -0.1035E-06_r8, -0.1079E-06_r8, -0.9742E-07_r8, -0.5928E-07_r8,  0.2892E-07_r8,  &
        0.1998E-06_r8,  0.4789E-06_r8,  0.8298E-06_r8,  0.1172E-05_r8,  0.1583E-05_r8,  &
        0.2329E-05_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip=15,15)/  &
        0.99999607_r8,  0.99999225_r8,  0.99998498_r8,  0.99997115_r8,  0.99994600_r8,  &
        0.99990171_r8,  0.99982780_r8,  0.99971092_r8,  0.99953562_r8,  0.99928278_r8,  &
        0.99892598_r8,  0.99842799_r8,  0.99773300_r8,  0.99675500_r8,  0.99536800_r8,  &
        0.99339402_r8,  0.99058902_r8,  0.98662001_r8,  0.98104000_r8,  0.97325999_r8,  &
        0.96249002_r8,  0.94773000_r8,  0.92778003_r8,  0.90125000_r8,  0.86680001_r8,  &
        0.82319999_r8,  0.77010000_r8,  0.70730001_r8,  0.63400000_r8,  0.54960001_r8,  &
        0.45560002_r8,  &  
       ! data ((h22(ip,iw),iw=1,31)_r8, ip=15,15)/  &
       -0.5736E-07_r8, -0.1109E-06_r8, -0.2106E-06_r8, -0.3890E-06_r8, -0.6928E-06_r8,  &
       -0.1181E-05_r8, -0.1917E-05_r8, -0.2965E-05_r8, -0.4396E-05_r8, -0.6315E-05_r8,  &
       -0.8891E-05_r8, -0.1242E-04_r8, -0.1736E-04_r8, -0.2442E-04_r8, -0.3454E-04_r8,  &
       -0.4892E-04_r8, -0.6916E-04_r8, -0.9735E-04_r8, -0.1362E-03_r8, -0.1891E-03_r8,  &
       -0.2602E-03_r8, -0.3545E-03_r8, -0.4768E-03_r8, -0.6310E-03_r8, -0.8179E-03_r8,  &
       -0.1032E-02_r8, -0.1265E-02_r8, -0.1508E-02_r8, -0.1757E-02_r8, -0.1989E-02_r8,  &
       -0.2159E-02_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip=15,15)/  &
       -0.2321E-09_r8, -0.4350E-09_r8, -0.7861E-09_r8, -0.1347E-08_r8, -0.2144E-08_r8,  &
       -0.3120E-08_r8, -0.4107E-08_r8, -0.4892E-08_r8, -0.5511E-08_r8, -0.6164E-08_r8,  &
       -0.7054E-08_r8, -0.8811E-08_r8, -0.1240E-07_r8, -0.1862E-07_r8, -0.2773E-07_r8,  &
       -0.3957E-07_r8, -0.5371E-07_r8, -0.6939E-07_r8, -0.8564E-07_r8, -0.1006E-06_r8,  &
       -0.1100E-06_r8, -0.1066E-06_r8, -0.8018E-07_r8, -0.1228E-07_r8,  0.1263E-06_r8,  &
        0.3678E-06_r8,  0.7022E-06_r8,  0.1049E-05_r8,  0.1411E-05_r8,  0.2015E-05_r8,  &
        0.3099E-05_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip=16,16)/  &
        0.99999791_r8,  0.99999589_r8,  0.99999183_r8,  0.99998391_r8,  0.99996853_r8,  &
        0.99993920_r8,  0.99988472_r8,  0.99978709_r8,  0.99961978_r8,  0.99934620_r8,  &
        0.99892199_r8,  0.99830103_r8,  0.99742401_r8,  0.99619502_r8,  0.99445999_r8,  &
        0.99199599_r8,  0.98851001_r8,  0.98360002_r8,  0.97671002_r8,  0.96710002_r8,  &
        0.95378000_r8,  0.93559003_r8,  0.91112000_r8,  0.87900001_r8,  0.83780003_r8,  &
        0.78680003_r8,  0.72549999_r8,  0.65300000_r8,  0.56830001_r8,  0.47140002_r8,  &
        0.36650002_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip=16,16)/  &
       -0.3122E-07_r8, -0.6175E-07_r8, -0.1214E-06_r8, -0.2361E-06_r8, -0.4518E-06_r8,  &
       -0.8438E-06_r8, -0.1524E-05_r8, -0.2643E-05_r8, -0.4380E-05_r8, -0.6922E-05_r8,  &
       -0.1042E-04_r8, -0.1504E-04_r8, -0.2125E-04_r8, -0.2987E-04_r8, -0.4200E-04_r8,  &
       -0.5923E-04_r8, -0.8383E-04_r8, -0.1186E-03_r8, -0.1670E-03_r8, -0.2328E-03_r8,  &
       -0.3204E-03_r8, -0.4347E-03_r8, -0.5802E-03_r8, -0.7595E-03_r8, -0.9703E-03_r8,  &
       -0.1205E-02_r8, -0.1457E-02_r8, -0.1720E-02_r8, -0.1980E-02_r8, -0.2191E-02_r8,  &
       -0.2290E-02_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip=16,16)/  &
       -0.1376E-09_r8, -0.2699E-09_r8, -0.5220E-09_r8, -0.9897E-09_r8, -0.1819E-08_r8,  &
       -0.3186E-08_r8, -0.5224E-08_r8, -0.7896E-08_r8, -0.1090E-07_r8, -0.1349E-07_r8,  &
       -0.1443E-07_r8, -0.1374E-07_r8, -0.1386E-07_r8, -0.1673E-07_r8, -0.2237E-07_r8,  &
       -0.3248E-07_r8, -0.5050E-07_r8, -0.7743E-07_r8, -0.1097E-06_r8, -0.1369E-06_r8,  &
       -0.1463E-06_r8, -0.1268E-06_r8, -0.6424E-07_r8,  0.5941E-07_r8,  0.2742E-06_r8,  &
        0.5924E-06_r8,  0.9445E-06_r8,  0.1286E-05_r8,  0.1819E-05_r8,  0.2867E-05_r8,  &
        0.4527E-05_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip=17,17)/  &
        0.99999756_r8,  0.99999511_r8,  0.99999028_r8,  0.99998081_r8,  0.99996233_r8,  &
        0.99992681_r8,  0.99986011_r8,  0.99973929_r8,  0.99953061_r8,  0.99918979_r8,  &
        0.99866599_r8,  0.99790198_r8,  0.99681997_r8,  0.99528998_r8,  0.99312103_r8,  &
        0.99004799_r8,  0.98571002_r8,  0.97961998_r8,  0.97105998_r8,  0.95915002_r8,  &
        0.94278002_r8,  0.92061001_r8,  0.89120001_r8,  0.85320002_r8,  0.80550003_r8,  &
        0.74759996_r8,  0.67900002_r8,  0.59829998_r8,  0.50520003_r8,  0.40219998_r8,  &
        0.29600000_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip=17,17)/  &
       -0.3547E-07_r8, -0.7029E-07_r8, -0.1386E-06_r8, -0.2709E-06_r8, -0.5218E-06_r8,  &
       -0.9840E-06_r8, -0.1799E-05_r8, -0.3156E-05_r8, -0.5272E-05_r8, -0.8357E-05_r8,  &
       -0.1260E-04_r8, -0.1827E-04_r8, -0.2598E-04_r8, -0.3667E-04_r8, -0.5169E-04_r8,  &
       -0.7312E-04_r8, -0.1037E-03_r8, -0.1467E-03_r8, -0.2060E-03_r8, -0.2857E-03_r8,  &
       -0.3907E-03_r8, -0.5257E-03_r8, -0.6940E-03_r8, -0.8954E-03_r8, -0.1124E-02_r8,  &
       -0.1371E-02_r8, -0.1632E-02_r8, -0.1897E-02_r8, -0.2131E-02_r8, -0.2275E-02_r8,  &
       -0.2265E-02_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip=17,17)/  &
       -0.1482E-09_r8, -0.2910E-09_r8, -0.5667E-09_r8, -0.1081E-08_r8, -0.2005E-08_r8,  &
       -0.3554E-08_r8, -0.5902E-08_r8, -0.8925E-08_r8, -0.1209E-07_r8, -0.1448E-07_r8,  &
       -0.1536E-07_r8, -0.1565E-07_r8, -0.1763E-07_r8, -0.2088E-07_r8, -0.2564E-07_r8,  &
       -0.3635E-07_r8, -0.5791E-07_r8, -0.8907E-07_r8, -0.1213E-06_r8, -0.1418E-06_r8,  &
       -0.1397E-06_r8, -0.1000E-06_r8, -0.4427E-08_r8,  0.1713E-06_r8,  0.4536E-06_r8,  &
        0.8086E-06_r8,  0.1153E-05_r8,  0.1588E-05_r8,  0.2437E-05_r8,  0.3905E-05_r8,  &
        0.5874E-05_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip=18,18)/  &
        0.99999714_r8,  0.99999428_r8,  0.99998862_r8,  0.99997741_r8,  0.99995553_r8,  &
        0.99991333_r8,  0.99983358_r8,  0.99968803_r8,  0.99943441_r8,  0.99901879_r8,  &
        0.99837899_r8,  0.99744099_r8,  0.99609798_r8,  0.99418801_r8,  0.99147803_r8,  &
        0.98764998_r8,  0.98227000_r8,  0.97469997_r8,  0.96410000_r8,  0.94941998_r8,  &
        0.92940998_r8,  0.90263999_r8,  0.86769998_r8,  0.82330000_r8,  0.76880002_r8,  &
        0.70379996_r8,  0.62720001_r8,  0.53810000_r8,  0.43769997_r8,  0.33149999_r8,  &
        0.22839999_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip=18,18)/  &
       -0.4064E-07_r8, -0.8066E-07_r8, -0.1593E-06_r8, -0.3124E-06_r8, -0.6049E-06_r8,  &
       -0.1148E-05_r8, -0.2118E-05_r8, -0.3751E-05_r8, -0.6314E-05_r8, -0.1006E-04_r8,  &
       -0.1526E-04_r8, -0.2232E-04_r8, -0.3196E-04_r8, -0.4525E-04_r8, -0.6394E-04_r8,  &
       -0.9058E-04_r8, -0.1284E-03_r8, -0.1812E-03_r8, -0.2533E-03_r8, -0.3493E-03_r8,  &
       -0.4740E-03_r8, -0.6315E-03_r8, -0.8228E-03_r8, -0.1044E-02_r8, -0.1286E-02_r8,  &
       -0.1544E-02_r8, -0.1810E-02_r8, -0.2061E-02_r8, -0.2243E-02_r8, -0.2291E-02_r8,  &
       -0.2152E-02_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip=18,18)/  &
       -0.1630E-09_r8, -0.3213E-09_r8, -0.6266E-09_r8, -0.1201E-08_r8, -0.2248E-08_r8,  &
       -0.4030E-08_r8, -0.6770E-08_r8, -0.1033E-07_r8, -0.1392E-07_r8, -0.1640E-07_r8,  &
       -0.1768E-07_r8, -0.1932E-07_r8, -0.2229E-07_r8, -0.2508E-07_r8, -0.2940E-07_r8,  &
       -0.4200E-07_r8, -0.6717E-07_r8, -0.1002E-06_r8, -0.1286E-06_r8, -0.1402E-06_r8,  &
       -0.1216E-06_r8, -0.5487E-07_r8,  0.8418E-07_r8,  0.3246E-06_r8,  0.6610E-06_r8,  &
        0.1013E-05_r8,  0.1394E-05_r8,  0.2073E-05_r8,  0.3337E-05_r8,  0.5175E-05_r8,  &
        0.7255E-05_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip=19,19)/  &
        0.99999672_r8,  0.99999344_r8,  0.99998701_r8,  0.99997419_r8,  0.99994916_r8,  &
        0.99990064_r8,  0.99980861_r8,  0.99963921_r8,  0.99934143_r8,  0.99884701_r8,  &
        0.99807602_r8,  0.99692601_r8,  0.99525797_r8,  0.99287403_r8,  0.98948997_r8,  &
        0.98474002_r8,  0.97804999_r8,  0.96867001_r8,  0.95559001_r8,  0.93761998_r8,  &
        0.91336000_r8,  0.88139999_r8,  0.84029996_r8,  0.78920001_r8,  0.72770000_r8,  &
        0.65499997_r8,  0.56999999_r8,  0.47280002_r8,  0.36760002_r8,  0.26220000_r8,  &
        0.16700000_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip=19,19)/  &
       -0.4629E-07_r8, -0.9195E-07_r8, -0.1819E-06_r8, -0.3572E-06_r8, -0.6936E-06_r8,  &
       -0.1323E-05_r8, -0.2456E-05_r8, -0.4385E-05_r8, -0.7453E-05_r8, -0.1200E-04_r8,  &
       -0.1843E-04_r8, -0.2731E-04_r8, -0.3943E-04_r8, -0.5606E-04_r8, -0.7936E-04_r8,  &
       -0.1123E-03_r8, -0.1588E-03_r8, -0.2231E-03_r8, -0.3101E-03_r8, -0.4247E-03_r8,  &
       -0.5713E-03_r8, -0.7522E-03_r8, -0.9651E-03_r8, -0.1202E-02_r8, -0.1456E-02_r8,  &
       -0.1721E-02_r8, -0.1983E-02_r8, -0.2196E-02_r8, -0.2296E-02_r8, -0.2224E-02_r8,  &
       -0.1952E-02_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip=19,19)/  &
       -0.1827E-09_r8, -0.3607E-09_r8, -0.7057E-09_r8, -0.1359E-08_r8, -0.2552E-08_r8,  &
       -0.4615E-08_r8, -0.7854E-08_r8, -0.1218E-07_r8, -0.1670E-07_r8, -0.2008E-07_r8,  &
       -0.2241E-07_r8, -0.2516E-07_r8, -0.2796E-07_r8, -0.3015E-07_r8, -0.3506E-07_r8,  &
       -0.4958E-07_r8, -0.7627E-07_r8, -0.1070E-06_r8, -0.1289E-06_r8, -0.1286E-06_r8,  &
       -0.8843E-07_r8,  0.1492E-07_r8,  0.2118E-06_r8,  0.5155E-06_r8,  0.8665E-06_r8,  &
        0.1220E-05_r8,  0.1765E-05_r8,  0.2825E-05_r8,  0.4498E-05_r8,  0.6563E-05_r8,  &
        0.8422E-05_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip=20,20)/  &
        0.99999636_r8,  0.99999279_r8,  0.99998569_r8,  0.99997163_r8,  0.99994397_r8,  &
        0.99989033_r8,  0.99978799_r8,  0.99959832_r8,  0.99926043_r8,  0.99868900_r8,  &
        0.99777400_r8,  0.99637598_r8,  0.99431503_r8,  0.99134803_r8,  0.98714000_r8,  &
        0.98122001_r8,  0.97290999_r8,  0.96131998_r8,  0.94528997_r8,  0.92346001_r8,  &
        0.89429998_r8,  0.85650003_r8,  0.80879998_r8,  0.75080001_r8,  0.68190002_r8,  &
        0.60100001_r8,  0.50740004_r8,  0.40399998_r8,  0.29729998_r8,  0.19730002_r8,  &
        0.11479998_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip=20,20)/  &
       -0.5164E-07_r8, -0.1026E-06_r8, -0.2031E-06_r8, -0.3994E-06_r8, -0.7771E-06_r8,  &
       -0.1488E-05_r8, -0.2776E-05_r8, -0.5001E-05_r8, -0.8610E-05_r8, -0.1411E-04_r8,  &
       -0.2209E-04_r8, -0.3328E-04_r8, -0.4860E-04_r8, -0.6954E-04_r8, -0.9861E-04_r8,  &
       -0.1393E-03_r8, -0.1961E-03_r8, -0.2738E-03_r8, -0.3778E-03_r8, -0.5132E-03_r8,  &
       -0.6831E-03_r8, -0.8868E-03_r8, -0.1118E-02_r8, -0.1368E-02_r8, -0.1632E-02_r8,  &
       -0.1899E-02_r8, -0.2136E-02_r8, -0.2282E-02_r8, -0.2273E-02_r8, -0.2067E-02_r8,  &
       -0.1679E-02_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip=20,20)/  &
       -0.2058E-09_r8, -0.4066E-09_r8, -0.7967E-09_r8, -0.1539E-08_r8, -0.2904E-08_r8,  &
       -0.5293E-08_r8, -0.9116E-08_r8, -0.1447E-07_r8, -0.2058E-07_r8, -0.2608E-07_r8,  &
       -0.3053E-07_r8, -0.3418E-07_r8, -0.3619E-07_r8, -0.3766E-07_r8, -0.4313E-07_r8,  &
       -0.5817E-07_r8, -0.8299E-07_r8, -0.1072E-06_r8, -0.1195E-06_r8, -0.1031E-06_r8,  &
       -0.3275E-07_r8,  0.1215E-06_r8,  0.3835E-06_r8,  0.7220E-06_r8,  0.1062E-05_r8,  &
        0.1504E-05_r8,  0.2367E-05_r8,  0.3854E-05_r8,  0.5842E-05_r8,  0.7875E-05_r8,  &
        0.9082E-05_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip=21,21)/  &
        0.99999619_r8,  0.99999237_r8,  0.99998480_r8,  0.99996990_r8,  0.99994045_r8,  &
        0.99988312_r8,  0.99977320_r8,  0.99956751_r8,  0.99919540_r8,  0.99855101_r8,  &
        0.99748802_r8,  0.99581498_r8,  0.99329299_r8,  0.98961997_r8,  0.98439002_r8,  &
        0.97702003_r8,  0.96671999_r8,  0.95245999_r8,  0.93291998_r8,  0.90660000_r8,  &
        0.87199998_r8,  0.82780004_r8,  0.77329999_r8,  0.70819998_r8,  0.63119996_r8,  &
        0.54159999_r8,  0.44059998_r8,  0.33380002_r8,  0.23000002_r8,  0.14029998_r8,  &
        0.07340002_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip=21,21)/  &
       -0.5584E-07_r8, -0.1110E-06_r8, -0.2198E-06_r8, -0.4329E-06_r8, -0.8444E-06_r8,  &
       -0.1623E-05_r8, -0.3049E-05_r8, -0.5551E-05_r8, -0.9714E-05_r8, -0.1627E-04_r8,  &
       -0.2609E-04_r8, -0.4015E-04_r8, -0.5955E-04_r8, -0.8603E-04_r8, -0.1223E-03_r8,  &
       -0.1724E-03_r8, -0.2413E-03_r8, -0.3346E-03_r8, -0.4578E-03_r8, -0.6155E-03_r8,  &
       -0.8087E-03_r8, -0.1033E-02_r8, -0.1279E-02_r8, -0.1540E-02_r8, -0.1811E-02_r8,  &
       -0.2065E-02_r8, -0.2251E-02_r8, -0.2301E-02_r8, -0.2163E-02_r8, -0.1828E-02_r8,  &
       -0.1365E-02_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip=21,21)/  &
       -0.2274E-09_r8, -0.4498E-09_r8, -0.8814E-09_r8, -0.1708E-08_r8, -0.3247E-08_r8,  &
       -0.5972E-08_r8, -0.1045E-07_r8, -0.1707E-07_r8, -0.2545E-07_r8, -0.3440E-07_r8,  &
       -0.4259E-07_r8, -0.4822E-07_r8, -0.5004E-07_r8, -0.5061E-07_r8, -0.5485E-07_r8,  &
       -0.6687E-07_r8, -0.8483E-07_r8, -0.9896E-07_r8, -0.9646E-07_r8, -0.5557E-07_r8,  &
        0.5765E-07_r8,  0.2752E-06_r8,  0.5870E-06_r8,  0.9188E-06_r8,  0.1291E-05_r8,  &
        0.1971E-05_r8,  0.3251E-05_r8,  0.5115E-05_r8,  0.7221E-05_r8,  0.8825E-05_r8,  &
        0.9032E-05_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip=22,22)/  &
        0.99999607_r8,  0.99999213_r8,  0.99998438_r8,  0.99996895_r8,  0.99993849_r8,  &
        0.99987900_r8,  0.99976391_r8,  0.99954629_r8,  0.99914569_r8,  0.99843502_r8,  &
        0.99722600_r8,  0.99526101_r8,  0.99221897_r8,  0.98771000_r8,  0.98123002_r8,  &
        0.97209001_r8,  0.95936000_r8,  0.94187999_r8,  0.91820002_r8,  0.88679999_r8,  &
        0.84609997_r8,  0.79530001_r8,  0.73379999_r8,  0.66090000_r8,  0.57529998_r8,  &
        0.47740000_r8,  0.37129998_r8,  0.26490003_r8,  0.16890001_r8,  0.09340000_r8,  &
        0.04310000_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip=22,22)/  &
       -0.5833E-07_r8, -0.1160E-06_r8, -0.2300E-06_r8, -0.4540E-06_r8, -0.8885E-06_r8,  &
       -0.1718E-05_r8, -0.3256E-05_r8, -0.6010E-05_r8, -0.1072E-04_r8, -0.1838E-04_r8,  &
       -0.3026E-04_r8, -0.4772E-04_r8, -0.7223E-04_r8, -0.1057E-03_r8, -0.1512E-03_r8,  &
       -0.2128E-03_r8, -0.2961E-03_r8, -0.4070E-03_r8, -0.5514E-03_r8, -0.7320E-03_r8,  &
       -0.9467E-03_r8, -0.1187E-02_r8, -0.1446E-02_r8, -0.1717E-02_r8, -0.1985E-02_r8,  &
       -0.2203E-02_r8, -0.2307E-02_r8, -0.2237E-02_r8, -0.1965E-02_r8, -0.1532E-02_r8,  &
       -0.1044E-02_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip=22,22)/  &
       -0.2426E-09_r8, -0.4805E-09_r8, -0.9447E-09_r8, -0.1841E-08_r8, -0.3519E-08_r8,  &
       -0.6565E-08_r8, -0.1172E-07_r8, -0.1979E-07_r8, -0.3095E-07_r8, -0.4443E-07_r8,  &
       -0.5821E-07_r8, -0.6868E-07_r8, -0.7282E-07_r8, -0.7208E-07_r8, -0.7176E-07_r8,  &
       -0.7562E-07_r8, -0.8110E-07_r8, -0.7934E-07_r8, -0.5365E-07_r8,  0.2483E-07_r8,  &
        0.1959E-06_r8,  0.4731E-06_r8,  0.7954E-06_r8,  0.1123E-05_r8,  0.1652E-05_r8,  &
        0.2711E-05_r8,  0.4402E-05_r8,  0.6498E-05_r8,  0.8392E-05_r8,  0.9154E-05_r8,  &
        0.8261E-05_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip=23,23)/  &
        0.99999601_r8,  0.99999207_r8,  0.99998420_r8,  0.99996859_r8,  0.99993771_r8,  &
        0.99987692_r8,  0.99975860_r8,  0.99953198_r8,  0.99910772_r8,  0.99833697_r8,  &
        0.99698901_r8,  0.99473202_r8,  0.99113101_r8,  0.98566002_r8,  0.97770000_r8,  &
        0.96640998_r8,  0.95076001_r8,  0.92943001_r8,  0.90092003_r8,  0.86370003_r8,  &
        0.81659997_r8,  0.75889999_r8,  0.69000000_r8,  0.60870004_r8,  0.51440001_r8,  &
        0.40990001_r8,  0.30190003_r8,  0.20060003_r8,  0.11680001_r8,  0.05760002_r8,  &
        0.02270001_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip=23,23)/  &
       -0.5929E-07_r8, -0.1180E-06_r8, -0.2344E-06_r8, -0.4638E-06_r8, -0.9118E-06_r8,  &
       -0.1775E-05_r8, -0.3401E-05_r8, -0.6375E-05_r8, -0.1160E-04_r8, -0.2039E-04_r8,  &
       -0.3444E-04_r8, -0.5575E-04_r8, -0.8641E-04_r8, -0.1287E-03_r8, -0.1856E-03_r8,  &
       -0.2615E-03_r8, -0.3618E-03_r8, -0.4928E-03_r8, -0.6594E-03_r8, -0.8621E-03_r8,  &
       -0.1095E-02_r8, -0.1349E-02_r8, -0.1618E-02_r8, -0.1894E-02_r8, -0.2140E-02_r8,  &
       -0.2293E-02_r8, -0.2290E-02_r8, -0.2085E-02_r8, -0.1696E-02_r8, -0.1212E-02_r8,  &
       -0.7506E-03_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip=23,23)/  &
       -0.2496E-09_r8, -0.4954E-09_r8, -0.9780E-09_r8, -0.1915E-08_r8, -0.3697E-08_r8,  &
       -0.6991E-08_r8, -0.1279E-07_r8, -0.2231E-07_r8, -0.3653E-07_r8, -0.5541E-07_r8,  &
       -0.7688E-07_r8, -0.9614E-07_r8, -0.1069E-06_r8, -0.1065E-06_r8, -0.9866E-07_r8,  &
       -0.8740E-07_r8, -0.7192E-07_r8, -0.4304E-07_r8,  0.1982E-07_r8,  0.1525E-06_r8,  &
        0.3873E-06_r8,  0.6947E-06_r8,  0.1000E-05_r8,  0.1409E-05_r8,  0.2253E-05_r8,  &
        0.3739E-05_r8,  0.5744E-05_r8,  0.7812E-05_r8,  0.9067E-05_r8,  0.8746E-05_r8,  &
        0.6940E-05_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip=24,24)/  &
        0.99999601_r8,  0.99999207_r8,  0.99998420_r8,  0.99996853_r8,  0.99993742_r8,  &
        0.99987602_r8,  0.99975550_r8,  0.99952233_r8,  0.99907869_r8,  0.99825698_r8,  &
        0.99678302_r8,  0.99424398_r8,  0.99007100_r8,  0.98356003_r8,  0.97387999_r8,  &
        0.96004999_r8,  0.94090003_r8,  0.91503000_r8,  0.88099998_r8,  0.83740002_r8,  &
        0.78350002_r8,  0.71869999_r8,  0.64170003_r8,  0.55149996_r8,  0.44950002_r8,  &
        0.34100002_r8,  0.23540002_r8,  0.14380002_r8,  0.07550001_r8,  0.03240001_r8,  &
        0.01029998_r8,  &  
       ! data ((h22(ip,iw),iw=1,31)_r8, ip=24,24)/  &
       -0.5950E-07_r8, -0.1185E-06_r8, -0.2358E-06_r8, -0.4678E-06_r8, -0.9235E-06_r8,  &
       -0.1810E-05_r8, -0.3503E-05_r8, -0.6664E-05_r8, -0.1236E-04_r8, -0.2223E-04_r8,  &
       -0.3849E-04_r8, -0.6396E-04_r8, -0.1017E-03_r8, -0.1545E-03_r8, -0.2257E-03_r8,  &
       -0.3192E-03_r8, -0.4399E-03_r8, -0.5933E-03_r8, -0.7824E-03_r8, -0.1005E-02_r8,  &
       -0.1251E-02_r8, -0.1516E-02_r8, -0.1794E-02_r8, -0.2060E-02_r8, -0.2257E-02_r8,  &
       -0.2318E-02_r8, -0.2186E-02_r8, -0.1853E-02_r8, -0.1386E-02_r8, -0.9021E-03_r8,  &
       -0.5050E-03_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip=24,24)/  &
       -0.2515E-09_r8, -0.5001E-09_r8, -0.9904E-09_r8, -0.1951E-08_r8, -0.3800E-08_r8,  &
       -0.7288E-08_r8, -0.1362E-07_r8, -0.2452E-07_r8, -0.4184E-07_r8, -0.6663E-07_r8,  &
       -0.9770E-07_r8, -0.1299E-06_r8, -0.1533E-06_r8, -0.1584E-06_r8, -0.1425E-06_r8,  &
       -0.1093E-06_r8, -0.5972E-07_r8,  0.1426E-07_r8,  0.1347E-06_r8,  0.3364E-06_r8,  &
        0.6209E-06_r8,  0.9169E-06_r8,  0.1243E-05_r8,  0.1883E-05_r8,  0.3145E-05_r8,  &
        0.5011E-05_r8,  0.7136E-05_r8,  0.8785E-05_r8,  0.9048E-05_r8,  0.7686E-05_r8,  &
        0.5368E-05_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip=25,25)/  &
        0.99999601_r8,  0.99999207_r8,  0.99998420_r8,  0.99996847_r8,  0.99993724_r8,  &
        0.99987543_r8,  0.99975342_r8,  0.99951530_r8,  0.99905682_r8,  0.99819201_r8,  &
        0.99660802_r8,  0.99381101_r8,  0.98908001_r8,  0.98148000_r8,  0.96991003_r8,  &
        0.95317000_r8,  0.92992002_r8,  0.89880002_r8,  0.85839999_r8,  0.80799997_r8,  &
        0.74689996_r8,  0.67420000_r8,  0.58850002_r8,  0.48979998_r8,  0.38200003_r8,  &
        0.27329999_r8,  0.17479998_r8,  0.09700000_r8,  0.04500002_r8,  0.01620001_r8,  &
        0.00349998_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip=25,25)/  &
       -0.5953E-07_r8, -0.1187E-06_r8, -0.2363E-06_r8, -0.4697E-06_r8, -0.9304E-06_r8,  &
       -0.1833E-05_r8, -0.3578E-05_r8, -0.6889E-05_r8, -0.1299E-04_r8, -0.2384E-04_r8,  &
       -0.4227E-04_r8, -0.7201E-04_r8, -0.1174E-03_r8, -0.1825E-03_r8, -0.2710E-03_r8,  &
       -0.3861E-03_r8, -0.5313E-03_r8, -0.7095E-03_r8, -0.9203E-03_r8, -0.1158E-02_r8,  &
       -0.1417E-02_r8, -0.1691E-02_r8, -0.1968E-02_r8, -0.2200E-02_r8, -0.2320E-02_r8,  &
       -0.2263E-02_r8, -0.1998E-02_r8, -0.1563E-02_r8, -0.1068E-02_r8, -0.6300E-03_r8,  &
       -0.3174E-03_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip=25,25)/  &
       -0.2520E-09_r8, -0.5016E-09_r8, -0.9963E-09_r8, -0.1971E-08_r8, -0.3867E-08_r8,  &
       -0.7500E-08_r8, -0.1427E-07_r8, -0.2634E-07_r8, -0.4656E-07_r8, -0.7753E-07_r8,  &
       -0.1196E-06_r8, -0.1683E-06_r8, -0.2113E-06_r8, -0.2309E-06_r8, -0.2119E-06_r8,  &
       -0.1518E-06_r8, -0.5200E-07_r8,  0.9274E-07_r8,  0.2973E-06_r8,  0.5714E-06_r8,  &
        0.8687E-06_r8,  0.1152E-05_r8,  0.1621E-05_r8,  0.2634E-05_r8,  0.4310E-05_r8,  &
        0.6418E-05_r8,  0.8347E-05_r8,  0.9162E-05_r8,  0.8319E-05_r8,  0.6209E-05_r8,  &
        0.3844E-05_r8,  &
       ! data ((h21(ip,iw),iw=1,31)_r8, ip=26,26)/  &
        0.99999601_r8,  0.99999207_r8,  0.99998420_r8,  0.99996847_r8,  0.99993718_r8,  &
        0.99987501_r8,  0.99975210_r8,  0.99951041_r8,  0.99904078_r8,  0.99814302_r8,  &
        0.99646801_r8,  0.99344200_r8,  0.98819000_r8,  0.97952002_r8,  0.96597999_r8,  &
        0.94600999_r8,  0.91812003_r8,  0.88099998_r8,  0.83359998_r8,  0.77569997_r8,  &
        0.70669997_r8,  0.62529999_r8,  0.53049999_r8,  0.42449999_r8,  0.31419998_r8,  &
        0.20969999_r8,  0.12269998_r8,  0.06089997_r8,  0.02420002_r8,  0.00660002_r8,  &
        0.00040001_r8,  &
       ! data ((h22(ip,iw),iw=1,31)_r8, ip=26,26)/  &
       -0.5954E-07_r8, -0.1187E-06_r8, -0.2366E-06_r8, -0.4709E-06_r8, -0.9349E-06_r8,  &
       -0.1849E-05_r8, -0.3632E-05_r8, -0.7058E-05_r8, -0.1350E-04_r8, -0.2521E-04_r8,  &
       -0.4564E-04_r8, -0.7958E-04_r8, -0.1329E-03_r8, -0.2114E-03_r8, -0.3200E-03_r8,  &
       -0.4611E-03_r8, -0.6353E-03_r8, -0.8409E-03_r8, -0.1072E-02_r8, -0.1324E-02_r8,  &
       -0.1592E-02_r8, -0.1871E-02_r8, -0.2127E-02_r8, -0.2297E-02_r8, -0.2312E-02_r8,  &
       -0.2123E-02_r8, -0.1738E-02_r8, -0.1247E-02_r8, -0.7744E-03_r8, -0.4117E-03_r8,  &
       -0.1850E-03_r8,  &
       ! data ((h23(ip,iw),iw=1,31)_r8, ip=26,26)/  &
       -0.2522E-09_r8, -0.5025E-09_r8, -0.9997E-09_r8, -0.1983E-08_r8, -0.3912E-08_r8,  &
       -0.7650E-08_r8, -0.1474E-07_r8, -0.2777E-07_r8, -0.5055E-07_r8, -0.8745E-07_r8,  &
       -0.1414E-06_r8, -0.2095E-06_r8, -0.2790E-06_r8, -0.3241E-06_r8, -0.3135E-06_r8,  &
       -0.2269E-06_r8, -0.5896E-07_r8,  0.1875E-06_r8,  0.4996E-06_r8,  0.8299E-06_r8,  &
        0.1115E-05_r8,  0.1467E-05_r8,  0.2236E-05_r8,  0.3672E-05_r8,  0.5668E-05_r8,  &
        0.7772E-05_r8,  0.9094E-05_r8,  0.8827E-05_r8,  0.7041E-05_r8,  0.4638E-05_r8,  &
        0.2539E-05_r8/),SHAPE=(/nx*3*nh/))

   it=0
   DO k=1,nx
      DO i=1,3
         DO j=1,nh
            it=it+1              
             IF(i==1) THEN
                !WRITE(*,'(a5,2e17.9) ' )'h21   ',data4(it),h21(k,j)
                h21(k,j) = data4(it)
             END IF
             IF(i==2) THEN
                !WRITE(*,'(a5,2e17.9) ' )'h22   ',data4(it),h22(k,j)
                h22(k,j) =data4(it)
             END IF
             IF(i==3) THEN
                !WRITE(*,'(a5,2e17.9) ' )'h23   ',data4(it),h23(k,j)
                h23(k,j)=data4(it)
             END IF
         END DO
      END DO
   END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      data5(1:nx*3*nh)=RESHAPE(SOURCE=(/&
       ! data ((h81(ip,iw),iw=1,31), ip= 1, 1)/  &
        0.99998659_r8,  0.99997360_r8,  0.99994862_r8,  0.99990171_r8,  0.99981678_r8,  &
        0.99967158_r8,  0.99944150_r8,  0.99910933_r8,  0.99867302_r8,  0.99814397_r8,  &
        0.99753898_r8,  0.99686199_r8,  0.99610198_r8,  0.99523401_r8,  0.99421698_r8,  &
        0.99299300_r8,  0.99147898_r8,  0.98958999_r8,  0.98721999_r8,  0.98430002_r8,  &
        0.98071998_r8,  0.97639000_r8,  0.97115999_r8,  0.96480000_r8,  0.95695001_r8,  &
        0.94713998_r8,  0.93469000_r8,  0.91873002_r8,  0.89810002_r8,  0.87129998_r8,  &
        0.83679998_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip= 1_r8, 1)/  &
       -0.5685E-08_r8, -0.1331E-07_r8, -0.3249E-07_r8, -0.8137E-07_r8, -0.2048E-06_r8,  &
       -0.4973E-06_r8, -0.1118E-05_r8, -0.2246E-05_r8, -0.3982E-05_r8, -0.6290E-05_r8,  &
       -0.9040E-05_r8, -0.1215E-04_r8, -0.1567E-04_r8, -0.1970E-04_r8, -0.2449E-04_r8,  &
       -0.3046E-04_r8, -0.3798E-04_r8, -0.4725E-04_r8, -0.5831E-04_r8, -0.7123E-04_r8,  &
       -0.8605E-04_r8, -0.1028E-03_r8, -0.1212E-03_r8, -0.1413E-03_r8, -0.1635E-03_r8,  &
       -0.1884E-03_r8, -0.2160E-03_r8, -0.2461E-03_r8, -0.2778E-03_r8, -0.3098E-03_r8,  &
       -0.3411E-03_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip= 1_r8, 1)/  &
        0.2169E-10_r8,  0.5237E-10_r8,  0.1296E-09_r8,  0.3204E-09_r8,  0.7665E-09_r8,  &
        0.1691E-08_r8,  0.3222E-08_r8,  0.5110E-08_r8,  0.6779E-08_r8,  0.7681E-08_r8,  &
        0.7378E-08_r8,  0.5836E-08_r8,  0.3191E-08_r8, -0.1491E-08_r8, -0.1022E-07_r8,  &
       -0.2359E-07_r8, -0.3957E-07_r8, -0.5553E-07_r8, -0.6927E-07_r8, -0.7849E-07_r8,  &
       -0.8139E-07_r8, -0.7853E-07_r8, -0.7368E-07_r8, -0.7220E-07_r8, -0.7780E-07_r8,  &
       -0.9091E-07_r8, -0.1038E-06_r8, -0.9929E-07_r8, -0.5422E-07_r8,  0.5379E-07_r8,  &
        0.2350E-06_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip= 2_r8, 2)/  &
        0.99998659_r8,  0.99997360_r8,  0.99994862_r8,  0.99990171_r8,  0.99981678_r8,  &
        0.99967158_r8,  0.99944139_r8,  0.99910921_r8,  0.99867302_r8,  0.99814397_r8,  &
        0.99753797_r8,  0.99686003_r8,  0.99609798_r8,  0.99522603_r8,  0.99420297_r8,  &
        0.99296701_r8,  0.99142998_r8,  0.98949999_r8,  0.98706001_r8,  0.98400998_r8,  &
        0.98021001_r8,  0.97552001_r8,  0.96976000_r8,  0.96262002_r8,  0.95367002_r8,  &
        0.94234002_r8,  0.92781997_r8,  0.90903997_r8,  0.88459998_r8,  0.85290003_r8,  &
        0.81200004_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip= 2_r8, 2)/  &
       -0.5684E-08_r8, -0.1331E-07_r8, -0.3248E-07_r8, -0.8133E-07_r8, -0.2047E-06_r8,  &
       -0.4971E-06_r8, -0.1117E-05_r8, -0.2245E-05_r8, -0.3981E-05_r8, -0.6287E-05_r8,  &
       -0.9035E-05_r8, -0.1215E-04_r8, -0.1565E-04_r8, -0.1967E-04_r8, -0.2444E-04_r8,  &
       -0.3036E-04_r8, -0.3780E-04_r8, -0.4694E-04_r8, -0.5779E-04_r8, -0.7042E-04_r8,  &
       -0.8491E-04_r8, -0.1013E-03_r8, -0.1196E-03_r8, -0.1399E-03_r8, -0.1625E-03_r8,  &
       -0.1879E-03_r8, -0.2163E-03_r8, -0.2474E-03_r8, -0.2803E-03_r8, -0.3140E-03_r8,  &
       -0.3478E-03_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip= 2_r8, 2)/  &
        0.2168E-10_r8,  0.5242E-10_r8,  0.1295E-09_r8,  0.3201E-09_r8,  0.7662E-09_r8,  &
        0.1690E-08_r8,  0.3220E-08_r8,  0.5106E-08_r8,  0.6776E-08_r8,  0.7673E-08_r8,  &
        0.7362E-08_r8,  0.5808E-08_r8,  0.3138E-08_r8, -0.1595E-08_r8, -0.1041E-07_r8,  &
       -0.2390E-07_r8, -0.4010E-07_r8, -0.5636E-07_r8, -0.7045E-07_r8, -0.7972E-07_r8,  &
       -0.8178E-07_r8, -0.7677E-07_r8, -0.6876E-07_r8, -0.6381E-07_r8, -0.6583E-07_r8,  &
       -0.7486E-07_r8, -0.8229E-07_r8, -0.7017E-07_r8, -0.1497E-07_r8,  0.1051E-06_r8,  &
        0.2990E-06_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip= 3_r8, 3)/  &
        0.99998659_r8,  0.99997360_r8,  0.99994862_r8,  0.99990171_r8,  0.99981678_r8,  &
        0.99967152_r8,  0.99944133_r8,  0.99910891_r8,  0.99867201_r8,  0.99814302_r8,  &
        0.99753499_r8,  0.99685597_r8,  0.99609101_r8,  0.99521297_r8,  0.99418002_r8,  &
        0.99292499_r8,  0.99135399_r8,  0.98935997_r8,  0.98681003_r8,  0.98356998_r8,  &
        0.97947001_r8,  0.97430998_r8,  0.96784997_r8,  0.95972002_r8,  0.94941002_r8,  &
        0.93620002_r8,  0.91912001_r8,  0.89690000_r8,  0.86790001_r8,  0.83020002_r8,  &
        0.78210002_r8,  &  
       ! data ((h82(ip,iw),iw=1,31)_r8, ip= 3_r8, 3)/  &
       -0.5682E-08_r8, -0.1330E-07_r8, -0.3247E-07_r8, -0.8129E-07_r8, -0.2046E-06_r8,  &
       -0.4968E-06_r8, -0.1117E-05_r8, -0.2244E-05_r8, -0.3978E-05_r8, -0.6283E-05_r8,  &
       -0.9027E-05_r8, -0.1213E-04_r8, -0.1563E-04_r8, -0.1963E-04_r8, -0.2436E-04_r8,  &
       -0.3021E-04_r8, -0.3754E-04_r8, -0.4649E-04_r8, -0.5709E-04_r8, -0.6940E-04_r8,  &
       -0.8359E-04_r8, -0.9986E-04_r8, -0.1182E-03_r8, -0.1388E-03_r8, -0.1620E-03_r8,  &
       -0.1882E-03_r8, -0.2175E-03_r8, -0.2498E-03_r8, -0.2843E-03_r8, -0.3203E-03_r8,  &
       -0.3573E-03_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip= 3_r8, 3)/  &
        0.2167E-10_r8,  0.5238E-10_r8,  0.1294E-09_r8,  0.3198E-09_r8,  0.7656E-09_r8,  &
        0.1688E-08_r8,  0.3217E-08_r8,  0.5104E-08_r8,  0.6767E-08_r8,  0.7661E-08_r8,  &
        0.7337E-08_r8,  0.5764E-08_r8,  0.3051E-08_r8, -0.1752E-08_r8, -0.1068E-07_r8,  &
       -0.2436E-07_r8, -0.4081E-07_r8, -0.5740E-07_r8, -0.7165E-07_r8, -0.8046E-07_r8,  &
       -0.8082E-07_r8, -0.7289E-07_r8, -0.6141E-07_r8, -0.5294E-07_r8, -0.5134E-07_r8,  &
       -0.5552E-07_r8, -0.5609E-07_r8, -0.3464E-07_r8,  0.3275E-07_r8,  0.1669E-06_r8,  &
        0.3745E-06_r8,  &  
       ! data ((h81(ip,iw),iw=1,31)_r8, ip= 4_r8, 4)/  &
        0.99998659_r8,  0.99997360_r8,  0.99994862_r8,  0.99990171_r8,  0.99981678_r8,  &
        0.99967140_r8,  0.99944109_r8,  0.99910849_r8,  0.99867100_r8,  0.99814099_r8,  &
        0.99753201_r8,  0.99685001_r8,  0.99607998_r8,  0.99519402_r8,  0.99414498_r8,  &
        0.99286002_r8,  0.99123698_r8,  0.98914999_r8,  0.98644000_r8,  0.98293000_r8,  &
        0.97842002_r8,  0.97263998_r8,  0.96529001_r8,  0.95592999_r8,  0.94392002_r8,  &
        0.92839003_r8,  0.90815997_r8,  0.88169998_r8,  0.84720004_r8,  0.80269998_r8,  &
        0.74629998_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip= 4_r8, 4)/  &
       -0.5680E-08_r8, -0.1329E-07_r8, -0.3243E-07_r8, -0.8121E-07_r8, -0.2044E-06_r8,  &
       -0.4963E-06_r8, -0.1115E-05_r8, -0.2242E-05_r8, -0.3974E-05_r8, -0.6276E-05_r8,  &
       -0.9015E-05_r8, -0.1211E-04_r8, -0.1559E-04_r8, -0.1956E-04_r8, -0.2423E-04_r8,  &
       -0.2999E-04_r8, -0.3716E-04_r8, -0.4588E-04_r8, -0.5618E-04_r8, -0.6818E-04_r8,  &
       -0.8218E-04_r8, -0.9847E-04_r8, -0.1171E-03_r8, -0.1382E-03_r8, -0.1621E-03_r8,  &
       -0.1892E-03_r8, -0.2197E-03_r8, -0.2535E-03_r8, -0.2902E-03_r8, -0.3293E-03_r8,  &
       -0.3700E-03_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip= 4_r8, 4)/  &
        0.2166E-10_r8,  0.5229E-10_r8,  0.1294E-09_r8,  0.3193E-09_r8,  0.7644E-09_r8,  &
        0.1686E-08_r8,  0.3213E-08_r8,  0.5092E-08_r8,  0.6753E-08_r8,  0.7640E-08_r8,  &
        0.7302E-08_r8,  0.5696E-08_r8,  0.2917E-08_r8, -0.1984E-08_r8, -0.1108E-07_r8,  &
       -0.2497E-07_r8, -0.4171E-07_r8, -0.5849E-07_r8, -0.7254E-07_r8, -0.8017E-07_r8,  &
       -0.7802E-07_r8, -0.6662E-07_r8, -0.5153E-07_r8, -0.3961E-07_r8, -0.3387E-07_r8,  &
       -0.3219E-07_r8, -0.2426E-07_r8,  0.8700E-08_r8,  0.9027E-07_r8,  0.2400E-06_r8,  &
        0.4623E-06_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip= 5_r8, 5)/  &
        0.99998659_r8,  0.99997360_r8,  0.99994862_r8,  0.99990165_r8,  0.99981672_r8,  &
        0.99967128_r8,  0.99944091_r8,  0.99910778_r8,  0.99866998_r8,  0.99813801_r8,  &
        0.99752700_r8,  0.99684101_r8,  0.99606299_r8,  0.99516302_r8,  0.99408901_r8,  &
        0.99276000_r8,  0.99105698_r8,  0.98882997_r8,  0.98588997_r8,  0.98202002_r8,  &
        0.97696000_r8,  0.97039002_r8,  0.96192998_r8,  0.95104003_r8,  0.93691999_r8,  &
        0.91851997_r8,  0.89440000_r8,  0.86290002_r8,  0.82200003_r8,  0.76969999_r8,  &
        0.70420003_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip= 5_r8, 5)/  &
       -0.5675E-08_r8, -0.1328E-07_r8, -0.3239E-07_r8, -0.8110E-07_r8, -0.2040E-06_r8,  &
       -0.4954E-06_r8, -0.1114E-05_r8, -0.2238E-05_r8, -0.3968E-05_r8, -0.6265E-05_r8,  &
       -0.8996E-05_r8, -0.1208E-04_r8, -0.1553E-04_r8, -0.1945E-04_r8, -0.2404E-04_r8,  &
       -0.2966E-04_r8, -0.3663E-04_r8, -0.4508E-04_r8, -0.5508E-04_r8, -0.6686E-04_r8,  &
       -0.8082E-04_r8, -0.9732E-04_r8, -0.1165E-03_r8, -0.1382E-03_r8, -0.1630E-03_r8,  &
       -0.1913E-03_r8, -0.2234E-03_r8, -0.2593E-03_r8, -0.2989E-03_r8, -0.3417E-03_r8,  &
       -0.3857E-03_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip= 5_r8, 5)/  &
        0.2163E-10_r8,  0.5209E-10_r8,  0.1291E-09_r8,  0.3186E-09_r8,  0.7626E-09_r8,  &
        0.1682E-08_r8,  0.3203E-08_r8,  0.5078E-08_r8,  0.6730E-08_r8,  0.7606E-08_r8,  &
        0.7246E-08_r8,  0.5592E-08_r8,  0.2735E-08_r8, -0.2325E-08_r8, -0.1162E-07_r8,  &
       -0.2576E-07_r8, -0.4268E-07_r8, -0.5938E-07_r8, -0.7262E-07_r8, -0.7827E-07_r8,  &
       -0.7297E-07_r8, -0.5786E-07_r8, -0.3930E-07_r8, -0.2373E-07_r8, -0.1295E-07_r8,  &
       -0.3728E-08_r8,  0.1465E-07_r8,  0.6114E-07_r8,  0.1590E-06_r8,  0.3257E-06_r8,  &
        0.5622E-06_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip= 6_r8, 6)/ &
        0.99998659_r8,  0.99997360_r8,  0.99994862_r8,  0.99990165_r8,  0.99981672_r8, &
        0.99967122_r8,  0.99944037_r8,  0.99910682_r8,  0.99866802_r8,  0.99813402_r8, &
        0.99751902_r8,  0.99682599_r8,  0.99603701_r8,  0.99511498_r8,  0.99400300_r8, &
        0.99260598_r8,  0.99078500_r8,  0.98835999_r8,  0.98510998_r8,  0.98075998_r8, &
        0.97499001_r8,  0.96741998_r8,  0.95757002_r8,  0.94476998_r8,  0.92804998_r8, &
        0.90613002_r8,  0.87739998_r8,  0.83990002_r8,  0.79159999_r8,  0.73049998_r8, &
        0.65540004_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip= 6_r8, 6)/ &
       -0.5671E-08_r8, -0.1326E-07_r8, -0.3234E-07_r8, -0.8091E-07_r8, -0.2035E-06_r8, &
       -0.4941E-06_r8, -0.1111E-05_r8, -0.2232E-05_r8, -0.3958E-05_r8, -0.6247E-05_r8, &
       -0.8966E-05_r8, -0.1202E-04_r8, -0.1544E-04_r8, -0.1929E-04_r8, -0.2377E-04_r8, &
       -0.2921E-04_r8, -0.3593E-04_r8, -0.4409E-04_r8, -0.5385E-04_r8, -0.6555E-04_r8, &
       -0.7965E-04_r8, -0.9656E-04_r8, -0.1163E-03_r8, -0.1390E-03_r8, -0.1649E-03_r8, &
       -0.1947E-03_r8, -0.2288E-03_r8, -0.2675E-03_r8, -0.3109E-03_r8, -0.3575E-03_r8, &
       -0.4039E-03_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip= 6_r8, 6)/ &
        0.2155E-10_r8,  0.5188E-10_r8,  0.1288E-09_r8,  0.3175E-09_r8,  0.7599E-09_r8, &
        0.1675E-08_r8,  0.3190E-08_r8,  0.5059E-08_r8,  0.6699E-08_r8,  0.7551E-08_r8, &
        0.7154E-08_r8,  0.5435E-08_r8,  0.2452E-08_r8, -0.2802E-08_r8, -0.1235E-07_r8, &
       -0.2668E-07_r8, -0.4353E-07_r8, -0.5962E-07_r8, -0.7134E-07_r8, -0.7435E-07_r8, &
       -0.6551E-07_r8, -0.4676E-07_r8, -0.2475E-07_r8, -0.4876E-08_r8,  0.1235E-07_r8, &
        0.3092E-07_r8,  0.6192E-07_r8,  0.1243E-06_r8,  0.2400E-06_r8,  0.4247E-06_r8, &
        0.6755E-06_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip= 7_r8, 7)/ &
        0.99998659_r8,  0.99997360_r8,  0.99994862_r8,  0.99990165_r8,  0.99981660_r8, &
        0.99967092_r8,  0.99943972_r8,  0.99910510_r8,  0.99866402_r8,  0.99812698_r8, &
        0.99750602_r8,  0.99680197_r8,  0.99599499_r8,  0.99504000_r8,  0.99387002_r8, &
        0.99237198_r8,  0.99038202_r8,  0.98768002_r8,  0.98400998_r8,  0.97903001_r8, &
        0.97236001_r8,  0.96354002_r8,  0.95196998_r8,  0.93681002_r8,  0.91688001_r8, &
        0.89069998_r8,  0.85640001_r8,  0.81190002_r8,  0.75520003_r8,  0.68470001_r8, &
        0.59990001_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip= 7_r8, 7)/ &
       -0.5665E-08_r8, -0.1322E-07_r8, -0.3224E-07_r8, -0.8063E-07_r8, -0.2027E-06_r8, &
       -0.4921E-06_r8, -0.1106E-05_r8, -0.2223E-05_r8, -0.3942E-05_r8, -0.6220E-05_r8, &
       -0.8920E-05_r8, -0.1194E-04_r8, -0.1530E-04_r8, -0.1905E-04_r8, -0.2337E-04_r8, &
       -0.2860E-04_r8, -0.3505E-04_r8, -0.4296E-04_r8, -0.5259E-04_r8, -0.6439E-04_r8, &
       -0.7884E-04_r8, -0.9635E-04_r8, -0.1170E-03_r8, -0.1407E-03_r8, -0.1681E-03_r8, &
       -0.1998E-03_r8, -0.2366E-03_r8, -0.2790E-03_r8, -0.3265E-03_r8, -0.3763E-03_r8, &
       -0.4235E-03_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip= 7_r8, 7)/  &
        0.2157E-10_r8,  0.5178E-10_r8,  0.1283E-09_r8,  0.3162E-09_r8,  0.7558E-09_r8,  &
        0.1665E-08_r8,  0.3169E-08_r8,  0.5027E-08_r8,  0.6645E-08_r8,  0.7472E-08_r8,  &
        0.7017E-08_r8,  0.5212E-08_r8,  0.2059E-08_r8, -0.3443E-08_r8, -0.1321E-07_r8,  &
       -0.2754E-07_r8, -0.4389E-07_r8, -0.5869E-07_r8, -0.6825E-07_r8, -0.6819E-07_r8,  &
       -0.5570E-07_r8, -0.3343E-07_r8, -0.7592E-08_r8,  0.1778E-07_r8,  0.4322E-07_r8,  &
        0.7330E-07_r8,  0.1190E-06_r8,  0.1993E-06_r8,  0.3348E-06_r8,  0.5376E-06_r8,  &
        0.8030E-06_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip= 8_r8, 8)/  &
        0.99998659_r8,  0.99997360_r8,  0.99994856_r8,  0.99990159_r8,  0.99981642_r8,  &
        0.99967051_r8,  0.99943858_r8,  0.99910247_r8,  0.99865901_r8,  0.99811602_r8,  &
        0.99748600_r8,  0.99676597_r8,  0.99592900_r8,  0.99492502_r8,  0.99366802_r8,  &
        0.99202400_r8,  0.98979002_r8,  0.98672003_r8,  0.98249000_r8,  0.97671002_r8,  &
        0.96890998_r8,  0.95854002_r8,  0.94483000_r8,  0.92675000_r8,  0.90289998_r8,  &
        0.87160003_r8,  0.83069998_r8,  0.77829999_r8,  0.71239996_r8,  0.63209999_r8,  &
        0.53839999_r8,  &  
       ! data ((h82(ip,iw),iw=1,31)_r8, ip= 8_r8, 8)/  &
       -0.5652E-08_r8, -0.1318E-07_r8, -0.3210E-07_r8, -0.8018E-07_r8, -0.2014E-06_r8,  &
       -0.4888E-06_r8, -0.1099E-05_r8, -0.2210E-05_r8, -0.3918E-05_r8, -0.6179E-05_r8,  &
       -0.8849E-05_r8, -0.1182E-04_r8, -0.1509E-04_r8, -0.1871E-04_r8, -0.2284E-04_r8,  &
       -0.2782E-04_r8, -0.3403E-04_r8, -0.4177E-04_r8, -0.5145E-04_r8, -0.6354E-04_r8,  &
       -0.7853E-04_r8, -0.9681E-04_r8, -0.1185E-03_r8, -0.1437E-03_r8, -0.1729E-03_r8,  &
       -0.2072E-03_r8, -0.2475E-03_r8, -0.2942E-03_r8, -0.3457E-03_r8, -0.3973E-03_r8,  &
       -0.4434E-03_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip= 8_r8, 8)/  &
        0.2153E-10_r8,  0.5151E-10_r8,  0.1273E-09_r8,  0.3136E-09_r8,  0.7488E-09_r8,  &
        0.1649E-08_r8,  0.3142E-08_r8,  0.4980E-08_r8,  0.6559E-08_r8,  0.7346E-08_r8,  &
        0.6813E-08_r8,  0.4884E-08_r8,  0.1533E-08_r8, -0.4209E-08_r8, -0.1409E-07_r8,  &
       -0.2801E-07_r8, -0.4320E-07_r8, -0.5614E-07_r8, -0.6312E-07_r8, -0.5976E-07_r8,  &
       -0.4369E-07_r8, -0.1775E-07_r8,  0.1280E-07_r8,  0.4534E-07_r8,  0.8106E-07_r8,  &
        0.1246E-06_r8,  0.1874E-06_r8,  0.2873E-06_r8,  0.4433E-06_r8,  0.6651E-06_r8,  &
        0.9477E-06_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip= 9_r8, 9)/  &
        0.99998659_r8,  0.99997360_r8,  0.99994856_r8,  0.99990153_r8,  0.99981618_r8,  &
        0.99966979_r8,  0.99943691_r8,  0.99909842_r8,  0.99865001_r8,  0.99809903_r8,  &
        0.99745399_r8,  0.99670798_r8,  0.99582899_r8,  0.99475002_r8,  0.99336600_r8,  &
        0.99151403_r8,  0.98896003_r8,  0.98540002_r8,  0.98044997_r8,  0.97364998_r8,  &
        0.96445000_r8,  0.95213997_r8,  0.93579000_r8,  0.91412002_r8,  0.88550001_r8,  &
        0.84810001_r8,  0.79970002_r8,  0.73850000_r8,  0.66280001_r8,  0.57309997_r8,  &
        0.47189999_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip= 9_r8, 9)/    &
       -0.5629E-08_r8, -0.1310E-07_r8, -0.3186E-07_r8, -0.7948E-07_r8, -0.1995E-06_r8,  &
       -0.4837E-06_r8, -0.1088E-05_r8, -0.2188E-05_r8, -0.3880E-05_r8, -0.6115E-05_r8,  &
       -0.8743E-05_r8, -0.1165E-04_r8, -0.1480E-04_r8, -0.1824E-04_r8, -0.2216E-04_r8,  &
       -0.2691E-04_r8, -0.3293E-04_r8, -0.4067E-04_r8, -0.5057E-04_r8, -0.6314E-04_r8,  &
       -0.7885E-04_r8, -0.9813E-04_r8, -0.1212E-03_r8, -0.1482E-03_r8, -0.1799E-03_r8,  &
       -0.2175E-03_r8, -0.2622E-03_r8, -0.3135E-03_r8, -0.3678E-03_r8, -0.4193E-03_r8,  &
       -0.4627E-03_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip= 9_r8, 9)/  &
        0.2121E-10_r8,  0.5076E-10_r8,  0.1257E-09_r8,  0.3091E-09_r8,  0.7379E-09_r8,  &
        0.1623E-08_r8,  0.3097E-08_r8,  0.4904E-08_r8,  0.6453E-08_r8,  0.7168E-08_r8,  &
        0.6534E-08_r8,  0.4458E-08_r8,  0.8932E-09_r8, -0.5026E-08_r8, -0.1469E-07_r8,  &
       -0.2765E-07_r8, -0.4103E-07_r8, -0.5169E-07_r8, -0.5585E-07_r8, -0.4913E-07_r8,  &
       -0.2954E-07_r8,  0.6372E-09_r8,  0.3738E-07_r8,  0.7896E-07_r8,  0.1272E-06_r8,  &
        0.1867E-06_r8,  0.2682E-06_r8,  0.3895E-06_r8,  0.5672E-06_r8,  0.8091E-06_r8,  &
        0.1114E-05_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip=10,10)/  &
        0.99998659_r8,  0.99997360_r8,  0.99994850_r8,  0.99990141_r8,  0.99981582_r8,  &
        0.99966878_r8,  0.99943417_r8,  0.99909198_r8,  0.99863601_r8,  0.99807203_r8,  &
        0.99740499_r8,  0.99662101_r8,  0.99567503_r8,  0.99448699_r8,  0.99292302_r8,  &
        0.99078500_r8,  0.98780000_r8,  0.98360002_r8,  0.97773999_r8,  0.96968001_r8,  &
        0.95872998_r8,  0.94401997_r8,  0.92440999_r8,  0.89840001_r8,  0.86409998_r8,  &
        0.81959999_r8,  0.76279998_r8,  0.69190001_r8,  0.60650003_r8,  0.50839996_r8,  &
        0.40249997_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip=10,10)/  &
       -0.5597E-08_r8, -0.1300E-07_r8, -0.3148E-07_r8, -0.7838E-07_r8, -0.1964E-06_r8,  &
       -0.4759E-06_r8, -0.1071E-05_r8, -0.2155E-05_r8, -0.3822E-05_r8, -0.6019E-05_r8,  &
       -0.8586E-05_r8, -0.1139E-04_r8, -0.1439E-04_r8, -0.1764E-04_r8, -0.2134E-04_r8,  &
       -0.2591E-04_r8, -0.3188E-04_r8, -0.3978E-04_r8, -0.5011E-04_r8, -0.6334E-04_r8,  &
       -0.7998E-04_r8, -0.1006E-03_r8, -0.1253E-03_r8, -0.1547E-03_r8, -0.1895E-03_r8,  &
       -0.2315E-03_r8, -0.2811E-03_r8, -0.3363E-03_r8, -0.3917E-03_r8, -0.4413E-03_r8,  &
       -0.4809E-03_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip=10,10)/  &
        0.2109E-10_r8,  0.5017E-10_r8,  0.1235E-09_r8,  0.3021E-09_r8,  0.7217E-09_r8,  &
        0.1585E-08_r8,  0.3028E-08_r8,  0.4796E-08_r8,  0.6285E-08_r8,  0.6910E-08_r8,  &
        0.6178E-08_r8,  0.3945E-08_r8,  0.2436E-09_r8, -0.5632E-08_r8, -0.1464E-07_r8,  &
       -0.2596E-07_r8, -0.3707E-07_r8, -0.4527E-07_r8, -0.4651E-07_r8, -0.3644E-07_r8,  &
       -0.1296E-07_r8,  0.2250E-07_r8,  0.6722E-07_r8,  0.1202E-06_r8,  0.1831E-06_r8,  &
        0.2605E-06_r8,  0.3627E-06_r8,  0.5062E-06_r8,  0.7064E-06_r8,  0.9725E-06_r8,  &
        0.1304E-05_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip=11,11)/ &
        0.99998659_r8,  0.99997354_r8,  0.99994850_r8,  0.99990124_r8,  0.99981529_r8, &
        0.99966723_r8,  0.99942988_r8,  0.99908209_r8,  0.99861503_r8,  0.99803102_r8, &
        0.99732900_r8,  0.99648702_r8,  0.99544603_r8,  0.99409997_r8,  0.99228698_r8, &
        0.98977000_r8,  0.98620999_r8,  0.98120999_r8,  0.97421998_r8,  0.96458000_r8, &
        0.95146000_r8,  0.93378001_r8,  0.91017002_r8,  0.87889999_r8,  0.83810002_r8, &
        0.78549999_r8,  0.71930003_r8,  0.63859999_r8,  0.54409999_r8,  0.43970001_r8, &
        0.33249998_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip=11,11)/ &
       -0.5538E-08_r8, -0.1280E-07_r8, -0.3089E-07_r8, -0.7667E-07_r8, -0.1917E-06_r8, &
       -0.4642E-06_r8, -0.1045E-05_r8, -0.2106E-05_r8, -0.3736E-05_r8, -0.5878E-05_r8, &
       -0.8363E-05_r8, -0.1104E-04_r8, -0.1387E-04_r8, -0.1692E-04_r8, -0.2044E-04_r8, &
       -0.2493E-04_r8, -0.3101E-04_r8, -0.3926E-04_r8, -0.5020E-04_r8, -0.6429E-04_r8, &
       -0.8213E-04_r8, -0.1044E-03_r8, -0.1314E-03_r8, -0.1637E-03_r8, -0.2027E-03_r8, &
       -0.2498E-03_r8, -0.3042E-03_r8, -0.3617E-03_r8, -0.4163E-03_r8, -0.4625E-03_r8, &
       -0.4969E-03_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip=11,11)/ &
        0.2067E-10_r8,  0.4903E-10_r8,  0.1200E-09_r8,  0.2917E-09_r8,  0.6965E-09_r8, &
        0.1532E-08_r8,  0.2925E-08_r8,  0.4632E-08_r8,  0.6054E-08_r8,  0.6590E-08_r8, &
        0.5746E-08_r8,  0.3436E-08_r8, -0.2251E-09_r8, -0.5703E-08_r8, -0.1344E-07_r8, &
       -0.2256E-07_r8, -0.3120E-07_r8, -0.3690E-07_r8, -0.3520E-07_r8, -0.2164E-07_r8, &
        0.6510E-08_r8,  0.4895E-07_r8,  0.1037E-06_r8,  0.1702E-06_r8,  0.2502E-06_r8, &
        0.3472E-06_r8,  0.4710E-06_r8,  0.6379E-06_r8,  0.8633E-06_r8,  0.1159E-05_r8, &
        0.1514E-05_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip=12,12)/ &
        0.99998659_r8,  0.99997354_r8,  0.99994838_r8,  0.99990094_r8,  0.99981439_r8, &
        0.99966472_r8,  0.99942350_r8,  0.99906689_r8,  0.99858302_r8,  0.99796802_r8, &
        0.99721497_r8,  0.99628800_r8,  0.99510801_r8,  0.99354398_r8,  0.99139601_r8, &
        0.98838001_r8,  0.98409998_r8,  0.97807997_r8,  0.96967000_r8,  0.95806998_r8, &
        0.94226003_r8,  0.92093998_r8,  0.89249998_r8,  0.85510004_r8,  0.80659997_r8, &
        0.74510002_r8,  0.66909999_r8,  0.57870001_r8,  0.47680002_r8,  0.36919999_r8, &
        0.26520002_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip=12,12)/ &
       -0.5476E-08_r8, -0.1257E-07_r8, -0.3008E-07_r8, -0.7418E-07_r8, -0.1848E-06_r8, &
       -0.4468E-06_r8, -0.1006E-05_r8, -0.2032E-05_r8, -0.3611E-05_r8, -0.5679E-05_r8, &
       -0.8058E-05_r8, -0.1059E-04_r8, -0.1324E-04_r8, -0.1612E-04_r8, -0.1956E-04_r8, &
       -0.2411E-04_r8, -0.3046E-04_r8, -0.3925E-04_r8, -0.5098E-04_r8, -0.6619E-04_r8, &
       -0.8562E-04_r8, -0.1100E-03_r8, -0.1399E-03_r8, -0.1761E-03_r8, -0.2202E-03_r8, &
       -0.2726E-03_r8, -0.3306E-03_r8, -0.3885E-03_r8, -0.4404E-03_r8, -0.4820E-03_r8, &
       -0.5082E-03_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip=12,12)/ &
        0.2041E-10_r8,  0.4771E-10_r8,  0.1149E-09_r8,  0.2782E-09_r8,  0.6614E-09_r8, &
        0.1451E-08_r8,  0.2778E-08_r8,  0.4401E-08_r8,  0.5736E-08_r8,  0.6189E-08_r8, &
        0.5315E-08_r8,  0.3087E-08_r8, -0.2518E-09_r8, -0.4806E-08_r8, -0.1071E-07_r8, &
       -0.1731E-07_r8, -0.2346E-07_r8, -0.2659E-07_r8, -0.2184E-07_r8, -0.4261E-08_r8, &
        0.2975E-07_r8,  0.8112E-07_r8,  0.1484E-06_r8,  0.2308E-06_r8,  0.3296E-06_r8, &
        0.4475E-06_r8,  0.5942E-06_r8,  0.7859E-06_r8,  0.1041E-05_r8,  0.1369E-05_r8, &
        0.1726E-05_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip=13,13)/ &
        0.99998653_r8,  0.99997348_r8,  0.99994826_r8,  0.99990052_r8,  0.99981320_r8, &
        0.99966109_r8,  0.99941391_r8,  0.99904412_r8,  0.99853402_r8,  0.99787498_r8, &
        0.99704498_r8,  0.99599600_r8,  0.99462402_r8,  0.99276501_r8,  0.99017602_r8, &
        0.98651999_r8,  0.98133999_r8,  0.97403997_r8,  0.96386999_r8,  0.94984001_r8, &
        0.93071002_r8,  0.90495998_r8,  0.87080002_r8,  0.82620001_r8,  0.76910001_r8, &
        0.69790000_r8,  0.61199999_r8,  0.51320004_r8,  0.40640002_r8,  0.29970002_r8, &
        0.20359999_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip=13,13)/ &
       -0.5362E-08_r8, -0.1223E-07_r8, -0.2895E-07_r8, -0.7071E-07_r8, -0.1748E-06_r8, &
       -0.4219E-06_r8, -0.9516E-06_r8, -0.1928E-05_r8, -0.3436E-05_r8, -0.5409E-05_r8, &
       -0.7666E-05_r8, -0.1005E-04_r8, -0.1254E-04_r8, -0.1533E-04_r8, -0.1880E-04_r8, &
       -0.2358E-04_r8, -0.3038E-04_r8, -0.3988E-04_r8, -0.5264E-04_r8, -0.6934E-04_r8, &
       -0.9083E-04_r8, -0.1179E-03_r8, -0.1515E-03_r8, -0.1927E-03_r8, -0.2424E-03_r8, &
       -0.2994E-03_r8, -0.3591E-03_r8, -0.4155E-03_r8, -0.4634E-03_r8, -0.4982E-03_r8, &
       -0.5096E-03_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip=13,13)/ &
        0.1976E-10_r8,  0.4551E-10_r8,  0.1086E-09_r8,  0.2601E-09_r8,  0.6126E-09_r8, &
        0.1345E-08_r8,  0.2583E-08_r8,  0.4112E-08_r8,  0.5365E-08_r8,  0.5796E-08_r8, &
        0.5031E-08_r8,  0.3182E-08_r8,  0.5970E-09_r8, -0.2547E-08_r8, -0.6172E-08_r8, &
       -0.1017E-07_r8, -0.1388E-07_r8, -0.1430E-07_r8, -0.6118E-08_r8,  0.1624E-07_r8, &
        0.5791E-07_r8,  0.1205E-06_r8,  0.2025E-06_r8,  0.3032E-06_r8,  0.4225E-06_r8, &
        0.5619E-06_r8,  0.7322E-06_r8,  0.9528E-06_r8,  0.1243E-05_r8,  0.1592E-05_r8, &
        0.1904E-05_r8,  & 
       ! data ((h81(ip,iw),iw=1,31)_r8, ip=14,14)/ &
        0.99998653_r8,  0.99997348_r8,  0.99994808_r8,  0.99989992_r8,  0.99981129_r8, &
        0.99965578_r8,  0.99939990_r8,  0.99901080_r8,  0.99846399_r8,  0.99773800_r8, &
        0.99680001_r8,  0.99558002_r8,  0.99394703_r8,  0.99169999_r8,  0.98853999_r8, &
        0.98408002_r8,  0.97776002_r8,  0.96888000_r8,  0.95652002_r8,  0.93949002_r8, &
        0.91631001_r8,  0.88529998_r8,  0.84439999_r8,  0.79159999_r8,  0.72510004_r8, &
        0.64390004_r8,  0.54890001_r8,  0.44379997_r8,  0.33560002_r8,  0.23449999_r8, &
        0.15009999_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip=14,14)/ &
       -0.5210E-08_r8, -0.1172E-07_r8, -0.2731E-07_r8, -0.6598E-07_r8, -0.1615E-06_r8, &
       -0.3880E-06_r8, -0.8769E-06_r8, -0.1787E-05_r8, -0.3204E-05_r8, -0.5066E-05_r8, &
       -0.7197E-05_r8, -0.9451E-05_r8, -0.1185E-04_r8, -0.1465E-04_r8, -0.1831E-04_r8, &
       -0.2346E-04_r8, -0.3088E-04_r8, -0.4132E-04_r8, -0.5545E-04_r8, -0.7410E-04_r8, &
       -0.9820E-04_r8, -0.1288E-03_r8, -0.1670E-03_r8, -0.2140E-03_r8, -0.2692E-03_r8, &
       -0.3293E-03_r8, -0.3886E-03_r8, -0.4417E-03_r8, -0.4840E-03_r8, -0.5073E-03_r8, &
       -0.4944E-03_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip=14,14)/  &
        0.1880E-10_r8,  0.4271E-10_r8,  0.9966E-10_r8,  0.2352E-09_r8,  0.5497E-09_r8,  &
        0.1205E-08_r8,  0.2334E-08_r8,  0.3765E-08_r8,  0.4993E-08_r8,  0.5532E-08_r8,  &
        0.5148E-08_r8,  0.4055E-08_r8,  0.2650E-08_r8,  0.1326E-08_r8,  0.2019E-09_r8,  &
       -0.1124E-08_r8, -0.2234E-08_r8,  0.2827E-09_r8,  0.1247E-07_r8,  0.4102E-07_r8,  &
        0.9228E-07_r8,  0.1682E-06_r8,  0.2676E-06_r8,  0.3885E-06_r8,  0.5286E-06_r8,  &
        0.6904E-06_r8,  0.8871E-06_r8,  0.1142E-05_r8,  0.1466E-05_r8,  0.1800E-05_r8,  &
        0.2004E-05_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip=15,15)/  &
        0.99998653_r8,  0.99997336_r8,  0.99994785_r8,  0.99989909_r8,  0.99980879_r8,  &
        0.99964851_r8,  0.99938041_r8,  0.99896401_r8,  0.99836302_r8,  0.99754399_r8,  &
        0.99645603_r8,  0.99500400_r8,  0.99302697_r8,  0.99027801_r8,  0.98640001_r8,  &
        0.98092002_r8,  0.97319001_r8,  0.96234000_r8,  0.94727999_r8,  0.92657000_r8,  &
        0.89850003_r8,  0.86119998_r8,  0.81260002_r8,  0.75080001_r8,  0.67429996_r8,  & 
        0.58350003_r8,  0.48089999_r8,  0.37250000_r8,  0.26760000_r8,  0.17650002_r8,  &
        0.10610002_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip=15,15)/  &
       -0.5045E-08_r8, -0.1113E-07_r8, -0.2540E-07_r8, -0.6008E-07_r8, -0.1449E-06_r8,  &
       -0.3457E-06_r8, -0.7826E-06_r8, -0.1609E-05_r8, -0.2920E-05_r8, -0.4665E-05_r8,  &
       -0.6691E-05_r8, -0.8868E-05_r8, -0.1127E-04_r8, -0.1422E-04_r8, -0.1820E-04_r8,  &
       -0.2389E-04_r8, -0.3213E-04_r8, -0.4380E-04_r8, -0.5975E-04_r8, -0.8092E-04_r8,  &
       -0.1083E-03_r8, -0.1433E-03_r8, -0.1873E-03_r8, -0.2402E-03_r8, -0.2997E-03_r8,  &
       -0.3607E-03_r8, -0.4178E-03_r8, -0.4662E-03_r8, -0.4994E-03_r8, -0.5028E-03_r8,  &
       -0.4563E-03_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip=15,15)/  &
        0.1804E-10_r8,  0.3983E-10_r8,  0.9045E-10_r8,  0.2080E-09_r8,  0.4786E-09_r8,  &
        0.1046E-08_r8,  0.2052E-08_r8,  0.3413E-08_r8,  0.4704E-08_r8,  0.5565E-08_r8,  &
        0.5887E-08_r8,  0.5981E-08_r8,  0.6202E-08_r8,  0.6998E-08_r8,  0.8493E-08_r8,  &
        0.1002E-07_r8,  0.1184E-07_r8,  0.1780E-07_r8,  0.3483E-07_r8,  0.7122E-07_r8,  &
        0.1341E-06_r8,  0.2259E-06_r8,  0.3446E-06_r8,  0.4866E-06_r8,  0.6486E-06_r8,  &
        0.8343E-06_r8,  0.1063E-05_r8,  0.1356E-05_r8,  0.1690E-05_r8,  0.1951E-05_r8,  &
        0.2005E-05_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip=16,16)/  &
        0.99998647_r8,  0.99997330_r8,  0.99994755_r8,  0.99989808_r8,  0.99980563_r8,  &
        0.99963909_r8,  0.99935490_r8,  0.99890202_r8,  0.99822801_r8,  0.99728203_r8,  &
        0.99599099_r8,  0.99423301_r8,  0.99181002_r8,  0.98842001_r8,  0.98364002_r8,  &
        0.97689003_r8,  0.96740001_r8,  0.95414001_r8,  0.93575001_r8,  0.91060001_r8,  &
        0.87680000_r8,  0.83219999_r8,  0.77490002_r8,  0.70330000_r8,  0.61689997_r8,  &
        0.51750004_r8,  0.40990001_r8,  0.30239999_r8,  0.20539999_r8,  0.12750000_r8,  &
        0.07150000_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip=16,16)/  &
       -0.4850E-08_r8, -0.1045E-07_r8, -0.2334E-07_r8, -0.5367E-07_r8, -0.1265E-06_r8,  &
       -0.2980E-06_r8, -0.6750E-06_r8, -0.1406E-05_r8, -0.2601E-05_r8, -0.4239E-05_r8,  &
       -0.6201E-05_r8, -0.8389E-05_r8, -0.1091E-04_r8, -0.1413E-04_r8, -0.1859E-04_r8,  &
       -0.2500E-04_r8, -0.3432E-04_r8, -0.4761E-04_r8, -0.6595E-04_r8, -0.9030E-04_r8,  &
       -0.1219E-03_r8, -0.1624E-03_r8, -0.2126E-03_r8, -0.2708E-03_r8, -0.3327E-03_r8,  &
       -0.3926E-03_r8, -0.4458E-03_r8, -0.4871E-03_r8, -0.5045E-03_r8, -0.4777E-03_r8,  &
       -0.3954E-03_r8,  &  
       ! data ((h83(ip,iw),iw=1,31)_r8, ip=16,16)/  &
        0.1717E-10_r8,  0.3723E-10_r8,  0.8093E-10_r8,  0.1817E-09_r8,  0.4100E-09_r8,  &
        0.8932E-09_r8,  0.1791E-08_r8,  0.3126E-08_r8,  0.4634E-08_r8,  0.6095E-08_r8,  &
        0.7497E-08_r8,  0.9170E-08_r8,  0.1136E-07_r8,  0.1453E-07_r8,  0.1892E-07_r8,  &
        0.2369E-07_r8,  0.2909E-07_r8,  0.3922E-07_r8,  0.6232E-07_r8,  0.1083E-06_r8,  &
        0.1847E-06_r8,  0.2943E-06_r8,  0.4336E-06_r8,  0.5970E-06_r8,  0.7815E-06_r8,  &
        0.9959E-06_r8,  0.1263E-05_r8,  0.1583E-05_r8,  0.1880E-05_r8,  0.2009E-05_r8,  &
        0.1914E-05_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip=17,17)/  &
        0.99998647_r8,  0.99997318_r8,  0.99994719_r8,  0.99989688_r8,  0.99980187_r8,  &
        0.99962789_r8,  0.99932390_r8,  0.99882400_r8,  0.99805701_r8,  0.99694502_r8,  &
        0.99538797_r8,  0.99323398_r8,  0.99023998_r8,  0.98604000_r8,  0.98013997_r8,  &
        0.97182000_r8,  0.96016002_r8,  0.94391000_r8,  0.92149997_r8,  0.89100003_r8,  &
        0.85049999_r8,  0.79769999_r8,  0.73089999_r8,  0.64919996_r8,  0.55350000_r8,  &
        0.44760001_r8,  0.33870000_r8,  0.23670000_r8,  0.15149999_r8,  0.08810002_r8,  &
        0.04570001_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip=17,17)/  &
       -0.4673E-08_r8, -0.9862E-08_r8, -0.2135E-07_r8, -0.4753E-07_r8, -0.1087E-06_r8,  &
       -0.2512E-06_r8, -0.5671E-06_r8, -0.1199E-05_r8, -0.2281E-05_r8, -0.3842E-05_r8,  &
       -0.5804E-05_r8, -0.8110E-05_r8, -0.1088E-04_r8, -0.1452E-04_r8, -0.1961E-04_r8,  &
       -0.2696E-04_r8, -0.3768E-04_r8, -0.5311E-04_r8, -0.7444E-04_r8, -0.1028E-03_r8,  &
       -0.1397E-03_r8, -0.1865E-03_r8, -0.2427E-03_r8, -0.3047E-03_r8, -0.3667E-03_r8,  &
       -0.4237E-03_r8, -0.4712E-03_r8, -0.5003E-03_r8, -0.4921E-03_r8, -0.4286E-03_r8,  &
       -0.3188E-03_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip=17,17)/  &
        0.1653E-10_r8,  0.3436E-10_r8,  0.7431E-10_r8,  0.1605E-09_r8,  0.3548E-09_r8,  &
        0.7723E-09_r8,  0.1595E-08_r8,  0.2966E-08_r8,  0.4849E-08_r8,  0.7169E-08_r8,  &
        0.1003E-07_r8,  0.1366E-07_r8,  0.1825E-07_r8,  0.2419E-07_r8,  0.3186E-07_r8,  &
        0.4068E-07_r8,  0.5064E-07_r8,  0.6618E-07_r8,  0.9684E-07_r8,  0.1536E-06_r8,  &
        0.2450E-06_r8,  0.3730E-06_r8,  0.5328E-06_r8,  0.7184E-06_r8,  0.9291E-06_r8,  &
        0.1180E-05_r8,  0.1484E-05_r8,  0.1798E-05_r8,  0.1992E-05_r8,  0.1968E-05_r8,  &
        0.1736E-05_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip=18,18)/  &
        0.99998647_r8,  0.99997312_r8,  0.99994683_r8,  0.99989569_r8,  0.99979800_r8,  &
        0.99961591_r8,  0.99928999_r8,  0.99873698_r8,  0.99785602_r8,  0.99653602_r8,  &
        0.99464101_r8,  0.99198103_r8,  0.98825997_r8,  0.98306000_r8,  0.97574002_r8,  &
        0.96548998_r8,  0.95117003_r8,  0.93129998_r8,  0.90407002_r8,  0.86739999_r8,  &
        0.81910002_r8,  0.75720000_r8,  0.68040001_r8,  0.58880001_r8,  0.48530000_r8,  &
        0.37610000_r8,  0.27029997_r8,  0.17830002_r8,  0.10720003_r8,  0.05790001_r8,  &
        0.02740002_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip=18,18)/  &
       -0.4532E-08_r8, -0.9395E-08_r8, -0.1978E-07_r8, -0.4272E-07_r8, -0.9442E-07_r8,  &
       -0.2124E-06_r8, -0.4747E-06_r8, -0.1017E-05_r8, -0.2003E-05_r8, -0.3524E-05_r8,  &
       -0.5567E-05_r8, -0.8108E-05_r8, -0.1127E-04_r8, -0.1547E-04_r8, -0.2138E-04_r8,  &
       -0.2996E-04_r8, -0.4251E-04_r8, -0.6059E-04_r8, -0.8563E-04_r8, -0.1190E-03_r8,  &
       -0.1623E-03_r8, -0.2156E-03_r8, -0.2767E-03_r8, -0.3403E-03_r8, -0.4006E-03_r8,  &
       -0.4530E-03_r8, -0.4912E-03_r8, -0.4995E-03_r8, -0.4563E-03_r8, -0.3592E-03_r8,  &
       -0.2383E-03_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip=18,18)/  &
        0.1593E-10_r8,  0.3276E-10_r8,  0.6896E-10_r8,  0.1476E-09_r8,  0.3190E-09_r8,  &
        0.6944E-09_r8,  0.1474E-08_r8,  0.2935E-08_r8,  0.5300E-08_r8,  0.8697E-08_r8,  &
        0.1336E-07_r8,  0.1946E-07_r8,  0.2707E-07_r8,  0.3637E-07_r8,  0.4800E-07_r8,  &
        0.6187E-07_r8,  0.7806E-07_r8,  0.1008E-06_r8,  0.1404E-06_r8,  0.2089E-06_r8,  &
        0.3153E-06_r8,  0.4613E-06_r8,  0.6416E-06_r8,  0.8506E-06_r8,  0.1095E-05_r8,  &
        0.1387E-05_r8,  0.1708E-05_r8,  0.1956E-05_r8,  0.2003E-05_r8,  0.1836E-05_r8,  &
        0.1483E-05_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip=19,19)/  &
        0.99998641_r8,  0.99997300_r8,  0.99994648_r8,  0.99989462_r8,  0.99979430_r8,  &
        0.99960452_r8,  0.99925661_r8,  0.99864697_r8,  0.99763900_r8,  0.99607199_r8,  &
        0.99376297_r8,  0.99046898_r8,  0.98584002_r8,  0.97937000_r8,  0.97031999_r8,  &
        0.95766997_r8,  0.94010001_r8,  0.91588002_r8,  0.88300002_r8,  0.83920002_r8,  &
        0.78230000_r8,  0.71060002_r8,  0.62360001_r8,  0.52320004_r8,  0.41450000_r8,  &
        0.30589998_r8,  0.20789999_r8,  0.12900001_r8,  0.07239997_r8,  0.03590000_r8,  &
        0.01539999_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip=19,19)/  &
       -0.4448E-08_r8, -0.9085E-08_r8, -0.1877E-07_r8, -0.3946E-07_r8, -0.8472E-07_r8,  &
       -0.1852E-06_r8, -0.4074E-06_r8, -0.8791E-06_r8, -0.1789E-05_r8, -0.3314E-05_r8,  &
       -0.5521E-05_r8, -0.8425E-05_r8, -0.1215E-04_r8, -0.1711E-04_r8, -0.2407E-04_r8,  &
       -0.3421E-04_r8, -0.4905E-04_r8, -0.7032E-04_r8, -0.9985E-04_r8, -0.1394E-03_r8,  &
       -0.1897E-03_r8, -0.2491E-03_r8, -0.3132E-03_r8, -0.3763E-03_r8, -0.4332E-03_r8,  &
       -0.4786E-03_r8, -0.5005E-03_r8, -0.4775E-03_r8, -0.3970E-03_r8, -0.2794E-03_r8,  &
       -0.1652E-03_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip=19,19)/  &
        0.1566E-10_r8,  0.3219E-10_r8,  0.6635E-10_r8,  0.1400E-09_r8,  0.2999E-09_r8,  &
        0.6513E-09_r8,  0.1406E-08_r8,  0.2953E-08_r8,  0.5789E-08_r8,  0.1037E-07_r8,  &
        0.1709E-07_r8,  0.2623E-07_r8,  0.3777E-07_r8,  0.5159E-07_r8,  0.6823E-07_r8,  &
        0.8864E-07_r8,  0.1134E-06_r8,  0.1461E-06_r8,  0.1960E-06_r8,  0.2761E-06_r8,  &
        0.3962E-06_r8,  0.5583E-06_r8,  0.7580E-06_r8,  0.9957E-06_r8,  0.1282E-05_r8,  &
        0.1607E-05_r8,  0.1898E-05_r8,  0.2020E-05_r8,  0.1919E-05_r8,  0.1623E-05_r8,  &
        0.1171E-05_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip=20,20)/  &
        0.99998641_r8,  0.99997294_r8,  0.99994624_r8,  0.99989372_r8,  0.99979132_r8,  &
        0.99959481_r8,  0.99922693_r8,  0.99856299_r8,  0.99742502_r8,  0.99558598_r8,  &
        0.99278802_r8,  0.98872000_r8,  0.98295999_r8,  0.97491002_r8,  0.96368998_r8,  &
        0.94812000_r8,  0.92662001_r8,  0.89719999_r8,  0.85780001_r8,  0.80599999_r8,  &
        0.73969996_r8,  0.65779996_r8,  0.56130004_r8,  0.45410001_r8,  0.34369999_r8,  &
        0.24030000_r8,  0.15390003_r8,  0.08950001_r8,  0.04640001_r8,  0.02090001_r8,  &
        0.00800002_r8,  &  
       ! data ((h82(ip,iw),iw=1,31)_r8, ip=20,20)/  &
       -0.4403E-08_r8, -0.8896E-08_r8, -0.1818E-07_r8, -0.3751E-07_r8, -0.7880E-07_r8,  &
       -0.1683E-06_r8, -0.3640E-06_r8, -0.7852E-06_r8, -0.1640E-05_r8, -0.3191E-05_r8,  &
       -0.5634E-05_r8, -0.9046E-05_r8, -0.1355E-04_r8, -0.1953E-04_r8, -0.2786E-04_r8,  &
       -0.3995E-04_r8, -0.5752E-04_r8, -0.8256E-04_r8, -0.1174E-03_r8, -0.1638E-03_r8,  &
       -0.2211E-03_r8, -0.2854E-03_r8, -0.3507E-03_r8, -0.4116E-03_r8, -0.4633E-03_r8,  &
       -0.4966E-03_r8, -0.4921E-03_r8, -0.4309E-03_r8, -0.3215E-03_r8, -0.2016E-03_r8,  &
       -0.1061E-03_r8,  &  
       ! data ((h83(ip,iw),iw=1,31)_r8, ip=20,20)/  &
        0.1551E-10_r8,  0.3147E-10_r8,  0.6419E-10_r8,  0.1356E-09_r8,  0.2860E-09_r8,  &
        0.6178E-09_r8,  0.1353E-08_r8,  0.2934E-08_r8,  0.6095E-08_r8,  0.1174E-07_r8,  &
        0.2067E-07_r8,  0.3346E-07_r8,  0.5014E-07_r8,  0.7024E-07_r8,  0.9377E-07_r8,  &
        0.1226E-06_r8,  0.1592E-06_r8,  0.2056E-06_r8,  0.2678E-06_r8,  0.3584E-06_r8,  &
        0.4892E-06_r8,  0.6651E-06_r8,  0.8859E-06_r8,  0.1160E-05_r8,  0.1488E-05_r8,  &
        0.1814E-05_r8,  0.2010E-05_r8,  0.1984E-05_r8,  0.1748E-05_r8,  0.1338E-05_r8,  &
        0.8445E-06_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip=21,21)/  &
        0.99998641_r8,  0.99997288_r8,  0.99994606_r8,  0.99989301_r8,  0.99978900_r8,  &
        0.99958712_r8,  0.99920273_r8,  0.99849200_r8,  0.99723101_r8,  0.99511403_r8,  &
        0.99177098_r8,  0.98677999_r8,  0.97962999_r8,  0.96961999_r8,  0.95573002_r8,  &
        0.93658000_r8,  0.91036999_r8,  0.87500000_r8,  0.82800001_r8,  0.76730001_r8,  &
        0.69110000_r8,  0.59930003_r8,  0.49479997_r8,  0.38370001_r8,  0.27590001_r8,  &
        0.18210000_r8,  0.10949999_r8,  0.05919999_r8,  0.02800000_r8,  0.01130003_r8,  &
        0.00389999_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip=21,21)/  &
       -0.4379E-08_r8, -0.8801E-08_r8, -0.1782E-07_r8, -0.3642E-07_r8, -0.7536E-07_r8,  &
       -0.1581E-06_r8, -0.3366E-06_r8, -0.7227E-06_r8, -0.1532E-05_r8, -0.3106E-05_r8,  &
       -0.5810E-05_r8, -0.9862E-05_r8, -0.1540E-04_r8, -0.2279E-04_r8, -0.3292E-04_r8,  &
       -0.4738E-04_r8, -0.6817E-04_r8, -0.9765E-04_r8, -0.1384E-03_r8, -0.1918E-03_r8,  &
       -0.2551E-03_r8, -0.3226E-03_r8, -0.3876E-03_r8, -0.4452E-03_r8, -0.4883E-03_r8,  &
       -0.5005E-03_r8, -0.4598E-03_r8, -0.3633E-03_r8, -0.2416E-03_r8, -0.1349E-03_r8,  &
       -0.6278E-04_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip=21,21)/  &
        0.1542E-10_r8,  0.3111E-10_r8,  0.6345E-10_r8,  0.1310E-09_r8,  0.2742E-09_r8,  &
        0.5902E-09_r8,  0.1289E-08_r8,  0.2826E-08_r8,  0.6103E-08_r8,  0.1250E-07_r8,  &
        0.2355E-07_r8,  0.4041E-07_r8,  0.6347E-07_r8,  0.9217E-07_r8,  0.1256E-06_r8,  &
        0.1658E-06_r8,  0.2175E-06_r8,  0.2824E-06_r8,  0.3607E-06_r8,  0.4614E-06_r8,  &
        0.6004E-06_r8,  0.7880E-06_r8,  0.1034E-05_r8,  0.1349E-05_r8,  0.1698E-05_r8,  &
        0.1965E-05_r8,  0.2021E-05_r8,  0.1857E-05_r8,  0.1500E-05_r8,  0.1015E-05_r8,  &
        0.5467E-06_r8,  &  
       ! data ((h81(ip,iw),iw=1,31)_r8, ip=22,22)/  &
        0.99998635_r8,  0.99997288_r8,  0.99994594_r8,  0.99989259_r8,  0.99978727_r8,  &
        0.99958128_r8,  0.99918407_r8,  0.99843502_r8,  0.99706697_r8,  0.99468601_r8,  &
        0.99077803_r8,  0.98474997_r8,  0.97593999_r8,  0.96350998_r8,  0.94633001_r8,  &
        0.92282999_r8,  0.89100003_r8,  0.84860003_r8,  0.79330003_r8,  0.72299999_r8,  &
        0.63670003_r8,  0.53600001_r8,  0.42580003_r8,  0.31470001_r8,  0.21410000_r8,  &  
        0.13300002_r8,  0.07470000_r8,  0.03710002_r8,  0.01580000_r8,  0.00580001_r8,  &
        0.00169998_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip=22,22)/  &
       -0.4366E-08_r8, -0.8749E-08_r8, -0.1761E-07_r8, -0.3578E-07_r8, -0.7322E-07_r8,  &
       -0.1517E-06_r8, -0.3189E-06_r8, -0.6785E-06_r8, -0.1446E-05_r8, -0.3014E-05_r8,  &
       -0.5933E-05_r8, -0.1069E-04_r8, -0.1755E-04_r8, -0.2683E-04_r8, -0.3936E-04_r8,  &
       -0.5675E-04_r8, -0.8137E-04_r8, -0.1160E-03_r8, -0.1630E-03_r8, -0.2223E-03_r8,  &
       -0.2899E-03_r8, -0.3589E-03_r8, -0.4230E-03_r8, -0.4755E-03_r8, -0.5031E-03_r8,  &
       -0.4834E-03_r8, -0.4036E-03_r8, -0.2849E-03_r8, -0.1687E-03_r8, -0.8356E-04_r8,  &
       -0.3388E-04_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip=22,22)/  &
        0.1536E-10_r8,  0.3086E-10_r8,  0.6248E-10_r8,  0.1288E-09_r8,  0.2664E-09_r8,  &
        0.5637E-09_r8,  0.1222E-08_r8,  0.2680E-08_r8,  0.5899E-08_r8,  0.1262E-07_r8,  &
        0.2527E-07_r8,  0.4621E-07_r8,  0.7678E-07_r8,  0.1165E-06_r8,  0.1640E-06_r8,  &
        0.2199E-06_r8,  0.2904E-06_r8,  0.3783E-06_r8,  0.4787E-06_r8,  0.5925E-06_r8,  &
        0.7377E-06_r8,  0.9389E-06_r8,  0.1216E-05_r8,  0.1560E-05_r8,  0.1879E-05_r8,  &
        0.2025E-05_r8,  0.1940E-05_r8,  0.1650E-05_r8,  0.1194E-05_r8,  0.6981E-06_r8,  &
        0.3103E-06_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip=23,23)/  &
        0.99998635_r8,  0.99997282_r8,  0.99994588_r8,  0.99989229_r8,  0.99978608_r8,  &
        0.99957722_r8,  0.99917048_r8,  0.99839097_r8,  0.99693698_r8,  0.99432403_r8,  &
        0.98987001_r8,  0.98273998_r8,  0.97201002_r8,  0.95668000_r8,  0.93548000_r8,  &
        0.90671998_r8,  0.86830002_r8,  0.81800002_r8,  0.75330001_r8,  0.67299998_r8,  &
        0.57720000_r8,  0.46920002_r8,  0.35659999_r8,  0.25010002_r8,  0.16049999_r8,  &
        0.09350002_r8,  0.04850000_r8,  0.02179998_r8,  0.00840002_r8,  0.00269997_r8,  &
        0.00070000_r8,  &  
       ! data ((h82(ip,iw),iw=1,31)_r8, ip=23,23)/  &
       -0.4359E-08_r8, -0.8720E-08_r8, -0.1749E-07_r8, -0.3527E-07_r8, -0.7175E-07_r8,  &
       -0.1473E-06_r8, -0.3062E-06_r8, -0.6451E-06_r8, -0.1372E-05_r8, -0.2902E-05_r8,  &
       -0.5936E-05_r8, -0.1133E-04_r8, -0.1971E-04_r8, -0.3143E-04_r8, -0.4715E-04_r8,  &
       -0.6833E-04_r8, -0.9759E-04_r8, -0.1379E-03_r8, -0.1907E-03_r8, -0.2542E-03_r8,  &
       -0.3239E-03_r8, -0.3935E-03_r8, -0.4559E-03_r8, -0.4991E-03_r8, -0.5009E-03_r8,  &
       -0.4414E-03_r8, -0.3306E-03_r8, -0.2077E-03_r8, -0.1093E-03_r8, -0.4754E-04_r8,  &
       -0.1642E-04_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip=23,23)/  &
        0.1531E-10_r8,  0.3070E-10_r8,  0.6184E-10_r8,  0.1257E-09_r8,  0.2578E-09_r8,  &
        0.5451E-09_r8,  0.1159E-08_r8,  0.2526E-08_r8,  0.5585E-08_r8,  0.1225E-07_r8,  &
        0.2576E-07_r8,  0.5017E-07_r8,  0.8855E-07_r8,  0.1417E-06_r8,  0.2078E-06_r8,  &
        0.2858E-06_r8,  0.3802E-06_r8,  0.4946E-06_r8,  0.6226E-06_r8,  0.7572E-06_r8,  &
        0.9137E-06_r8,  0.1133E-05_r8,  0.1438E-05_r8,  0.1772E-05_r8,  0.1994E-05_r8,  &
        0.1994E-05_r8,  0.1779E-05_r8,  0.1375E-05_r8,  0.8711E-06_r8,  0.4273E-06_r8,  &
        0.1539E-06_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip=24,24)/  &
        0.99998635_r8,  0.99997282_r8,  0.99994582_r8,  0.99989212_r8,  0.99978542_r8,  &
        0.99957442_r8,  0.99916071_r8,  0.99835902_r8,  0.99683702_r8,  0.99403203_r8,  &
        0.98908001_r8,  0.98084998_r8,  0.96805000_r8,  0.94933999_r8,  0.92330998_r8,  &
        0.88830000_r8,  0.84219998_r8,  0.78270000_r8,  0.70809996_r8,  0.61759996_r8,  &
        0.51330000_r8,  0.40079999_r8,  0.29000002_r8,  0.19239998_r8,  0.11619997_r8,  &
        0.06300002_r8,  0.02980000_r8,  0.01200002_r8,  0.00410002_r8,  0.00120002_r8,  &
        0.00019997_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip=24,24)/  &
       -0.4354E-08_r8, -0.8703E-08_r8, -0.1742E-07_r8, -0.3499E-07_r8, -0.7074E-07_r8,  &
       -0.1441E-06_r8, -0.2971E-06_r8, -0.6195E-06_r8, -0.1309E-05_r8, -0.2780E-05_r8,  &
       -0.5823E-05_r8, -0.1165E-04_r8, -0.2152E-04_r8, -0.3616E-04_r8, -0.5604E-04_r8,  &
       -0.8230E-04_r8, -0.1173E-03_r8, -0.1635E-03_r8, -0.2211E-03_r8, -0.2868E-03_r8,  &
       -0.3567E-03_r8, -0.4260E-03_r8, -0.4844E-03_r8, -0.5097E-03_r8, -0.4750E-03_r8,  &
       -0.3779E-03_r8, -0.2522E-03_r8, -0.1409E-03_r8, -0.6540E-04_r8, -0.2449E-04_r8,  &
       -0.6948E-05_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip=24,24)/  &
        0.1529E-10_r8,  0.3060E-10_r8,  0.6142E-10_r8,  0.1241E-09_r8,  0.2535E-09_r8,  &
        0.5259E-09_r8,  0.1107E-08_r8,  0.2383E-08_r8,  0.5243E-08_r8,  0.1161E-07_r8,  &
        0.2523E-07_r8,  0.5188E-07_r8,  0.9757E-07_r8,  0.1657E-06_r8,  0.2553E-06_r8,  &
        0.3629E-06_r8,  0.4878E-06_r8,  0.6323E-06_r8,  0.7923E-06_r8,  0.9575E-06_r8,  &
        0.1139E-05_r8,  0.1381E-05_r8,  0.1687E-05_r8,  0.1952E-05_r8,  0.2029E-05_r8,  &
        0.1890E-05_r8,  0.1552E-05_r8,  0.1062E-05_r8,  0.5728E-06_r8,  0.2280E-06_r8,  &
        0.6762E-07_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip=25,25)/  &
        0.99998635_r8,  0.99997282_r8,  0.99994582_r8,  0.99989200_r8,  0.99978489_r8,  &
        0.99957252_r8,  0.99915391_r8,  0.99833602_r8,  0.99676299_r8,  0.99380499_r8,  &
        0.98843998_r8,  0.97920001_r8,  0.96427000_r8,  0.94182003_r8,  0.91018999_r8,  &
        0.86769998_r8,  0.81260002_r8,  0.74300003_r8,  0.65770000_r8,  0.55750000_r8,  &
        0.44660002_r8,  0.33310002_r8,  0.22860003_r8,  0.14319998_r8,  0.08090001_r8,  &
        0.04030001_r8,  0.01719999_r8,  0.00620002_r8,  0.00190002_r8,  0.00040001_r8,  &
        0.00000000_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip=25,25)/  &
       -0.4352E-08_r8, -0.8693E-08_r8, -0.1738E-07_r8, -0.3483E-07_r8, -0.7006E-07_r8,  &
       -0.1423E-06_r8, -0.2905E-06_r8, -0.6008E-06_r8, -0.1258E-05_r8, -0.2663E-05_r8,  &
       -0.5638E-05_r8, -0.1165E-04_r8, -0.2270E-04_r8, -0.4044E-04_r8, -0.6554E-04_r8,  &
       -0.9855E-04_r8, -0.1407E-03_r8, -0.1928E-03_r8, -0.2534E-03_r8, -0.3197E-03_r8,  &
       -0.3890E-03_r8, -0.4563E-03_r8, -0.5040E-03_r8, -0.4998E-03_r8, -0.4249E-03_r8,  &
       -0.3025E-03_r8, -0.1794E-03_r8, -0.8860E-04_r8, -0.3575E-04_r8, -0.1122E-04_r8,  &
       -0.2506E-05_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip=25,25)/  &
        0.1527E-10_r8,  0.3053E-10_r8,  0.6115E-10_r8,  0.1230E-09_r8,  0.2492E-09_r8,  &
        0.5149E-09_r8,  0.1068E-08_r8,  0.2268E-08_r8,  0.4932E-08_r8,  0.1089E-07_r8,  &
        0.2408E-07_r8,  0.5156E-07_r8,  0.1028E-06_r8,  0.1859E-06_r8,  0.3028E-06_r8,  &
        0.4476E-06_r8,  0.6124E-06_r8,  0.7932E-06_r8,  0.9879E-06_r8,  0.1194E-05_r8,  &
        0.1417E-05_r8,  0.1673E-05_r8,  0.1929E-05_r8,  0.2064E-05_r8,  0.1997E-05_r8,  &
        0.1725E-05_r8,  0.1267E-05_r8,  0.7464E-06_r8,  0.3312E-06_r8,  0.1066E-06_r8,  &
        0.2718E-07_r8,  &
       ! data ((h81(ip,iw),iw=1,31)_r8, ip=26,26)/ &
        0.99998635_r8,  0.99997282_r8,  0.99994576_r8,  0.99989188_r8,  0.99978459_r8, &
        0.99957132_r8,  0.99914938_r8,  0.99831998_r8,  0.99670899_r8,  0.99363601_r8, &
        0.98794001_r8,  0.97781998_r8,  0.96087998_r8,  0.93456000_r8,  0.89670002_r8, &
        0.84560001_r8,  0.78020000_r8,  0.69920003_r8,  0.60299999_r8,  0.49400002_r8, &
        0.37910002_r8,  0.26889998_r8,  0.17460001_r8,  0.10280001_r8,  0.05379999_r8, &
        0.02429998_r8,  0.00929999_r8,  0.00300002_r8,  0.00080001_r8,  0.00010002_r8, &
        0.00000000_r8,  &
       ! data ((h82(ip,iw),iw=1,31)_r8, ip=26,26)/ &
       -0.4351E-08_r8, -0.8688E-08_r8, -0.1736E-07_r8, -0.3473E-07_r8, -0.6966E-07_r8, &
       -0.1405E-06_r8, -0.2857E-06_r8, -0.5867E-06_r8, -0.1218E-05_r8, -0.2563E-05_r8, &
       -0.5435E-05_r8, -0.1144E-04_r8, -0.2321E-04_r8, -0.4379E-04_r8, -0.7487E-04_r8, &
       -0.1163E-03_r8, -0.1670E-03_r8, -0.2250E-03_r8, -0.2876E-03_r8, -0.3535E-03_r8, &
       -0.4215E-03_r8, -0.4826E-03_r8, -0.5082E-03_r8, -0.4649E-03_r8, -0.3564E-03_r8, &
       -0.2264E-03_r8, -0.1188E-03_r8, -0.5128E-04_r8, -0.1758E-04_r8, -0.4431E-05_r8, &
       -0.7275E-06_r8,  &
       ! data ((h83(ip,iw),iw=1,31)_r8, ip=26,26)/ &
        0.1525E-10_r8,  0.3048E-10_r8,  0.6097E-10_r8,  0.1223E-09_r8,  0.2466E-09_r8, &
        0.5021E-09_r8,  0.1032E-08_r8,  0.2195E-08_r8,  0.4688E-08_r8,  0.1027E-07_r8, &
        0.2279E-07_r8,  0.4999E-07_r8,  0.1046E-06_r8,  0.2009E-06_r8,  0.3460E-06_r8, &
        0.5335E-06_r8,  0.7478E-06_r8,  0.9767E-06_r8,  0.1216E-05_r8,  0.1469E-05_r8, &
        0.1735E-05_r8,  0.1977E-05_r8,  0.2121E-05_r8,  0.2103E-05_r8,  0.1902E-05_r8, &
        0.1495E-05_r8,  0.9541E-06_r8,  0.4681E-06_r8,  0.1672E-06_r8,  0.4496E-07_r8, &
        0.9859E-08_r8/),SHAPE=(/nx*3*nh/))

    it=0
    DO k=1,nx
       DO i=1,3
          DO j=1,nh
             it=it+1              
              IF(i==1)THEN
                ! WRITE(*,'(a5,2e17.9) ' )'h81   ',data5(it),h81(k,j)
                 h81(k,j)=data5(it)
              END IF
              IF(i==2)THEN
                ! WRITE(*,'(a5,2e17.9) ' )'h82   ',data5(it),h82(k,j)
                 h82(k,j)=data5(it)
              END IF
              IF(i==3)THEN
                 !WRITE(*,'(a5,2e17.9) ' )'h83   ',data5(it),h83(k,j)
                 h83(k,j)=data5(it)
              END IF
          END DO
       END DO
    END DO

   END SUBROUTINE read_table

END MODULE Rad_Cliradlw


!PROGRAM Main
! USE Rad_Cliradlw, Only:InitCliRadLW
! CALL InitCliRadLW()
!END PROGRAM Main
