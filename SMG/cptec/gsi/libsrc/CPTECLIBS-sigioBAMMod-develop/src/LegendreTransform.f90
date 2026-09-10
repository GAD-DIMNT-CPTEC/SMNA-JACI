! -------------------------------------------------------------------------------
! INPE/CPTEC, DMD, Modelling and Development Division
! -------------------------------------------------------------------------------
!
! modULE: Mod_LegendreTransform
!
! REVISION HISTORY:
!            - Jose P. Bonatti         - version: 1
! 01-08-2007 - Tomita                  - version: 1.1.1.1
! 01-04-2018 - Daniel M. Lamosa        - version: 2.0
!              Barbara A. G. P. Yamada - version: 2.0
! 08-04-2019 - Eduardo Khamis
!> @author
!> Jose P. Bonatti \n
!> Daniel M. Lamosa (last revision) \n
!> Barbara A. G. P. Yamada (last revision) \n
!> Eduardo Khamis \n
!!
!> @brief 
!!
!!
!! @version 2.0
!! @date 01-04-2018
!! 
!! @copyright Under GLP-3.0
!! @link: https://opensource.org/licenses/GPL-3.0
! -------------------------------------------------------------------------------
module LegendreTransform
  use ModConstants, only: p_r8 => r8, p_i4 => i4
  use ModConstants, only: emRad
  use ModConstants, only: emRad1
  use ModConstants, only: emRad2
  use ModConstants, only: emRad12
  implicit none
  private

  integer, parameter :: p_nferr = 0

  type, public :: legendre
!     private

     logical :: created

     integer :: Mend
     integer :: Mend1
     integer :: Mend2
     integer :: Mend3
     integer :: Mends
     integer :: MnWv0
     integer :: MnWv1
     integer :: MnWv2
     integer :: MnWv3
     integer :: xMax
     integer :: yMax
     integer :: yMaxHf
     integer :: xMx
     
     integer(p_i4), pointer, dimension(:,:)   :: la0       => null()
     integer(p_i4), pointer, dimension(:,:)   :: la1       => null()
     real(kind=p_r8), pointer, dimension(:)   :: e0        => null()
     real(kind=p_r8), pointer, dimension(:)   :: e1        => null()
     real(kind=p_r8), pointer, dimension(:)   :: eps       => null()
     real(kind=p_r8), pointer, dimension(:)   :: colRad    => null()
     real(kind=p_r8), pointer, dimension(:)   :: rCs2      => null()
     real(kind=p_r8), pointer, dimension(:)   :: wgt       => null()
     real(kind=p_r8), pointer, dimension(:)   :: gLats     => null()
     real(kind=p_r8), pointer, dimension(:)   :: snnp      => null()
     real(kind=p_r8), pointer, dimension(:,:) :: legS2F    => null()
     real(kind=p_r8), pointer, dimension(:,:) :: legExtS2F => null()
     real(kind=p_r8), pointer, dimension(:,:) :: legDerS2F => null()
     real(kind=p_r8), pointer, dimension(:,:) :: legDerNS  => null()
     real(kind=p_r8), pointer, dimension(:,:) :: legDerEW  => null()
     contains
        procedure, public  :: initLegendre => initialization_
        procedure, public  :: destroyLegendre => destroylegendre
        procedure, public  :: transp => specTransp_
        procedure, public  :: SymAsy => SymAsy_

        procedure, public  :: Four2Spec => Fourier2SpecCoef_
        procedure, public  :: Spec2Four => SpecCoef2Fourier_

        procedure, public  :: GetSize

        procedure, private :: createSpectralRep
        procedure, private :: createGaussRep
        procedure, private :: createLegPols

        procedure, private :: DestroySpectralRep
        procedure, private :: DestroyGaussRep
        procedure, private :: DestroyLegPols

  end type


  contains
  function initialization_(self, mend) result(iret)

     class(legendre) :: self
     integer         :: Mend
     integer         :: iret

     iret = 0 ! precisa implementar isso
     
     self%Mend  = Mend
     self%Mend1 = Mend + 1
     self%Mend2 = Mend + 2
     self%Mend3 = Mend + 3
     self%Mends = 2*self%Mend1

     self%MnWv2 = (Mend + 1) * (Mend + 2)
     self%MnWv0 = self%MnWv2 / 2
     self%MnWv3 = self%MnWv2 + 2 * (self%Mend1)
     self%MnWv1 = self%MnWv3 / 2

     call getxMaxyMax(Mend, self%xMax, self%yMax)
     self%yMaxHf = self%yMax/2
     self%xmx    = self%xMax + 2

     self%created = .false.
     call self%createSpectralRep( )
     call self%createGaussRep( )
     call self%createLegPols( )

  end function

  function destroylegendre(self) result(iret)
     class(legendre) :: self
     integer         :: iret

     iret = 0 ! precisa implementar isso

     self%Mend  = -1 
     self%Mend1 = -1
     self%Mend2 = -1

     self%MnWv0 = -1
     self%MnWv1 = -1
     self%MnWv2 = -1
     self%MnWv3 = -1

     self%xMax  = -1
     self%yMax  = -1
     self%yMaxHf= -1

     self%created = .false.
     call self%destroySpectralRep( )
     call self%destroyGaussRep( )
     call self%destroyLegPols( )

  end function 

  ! ---------------------------------------------------------------------------
  !> @brief createSpectralRep
  !!
  !! @details 
  !!
  !! @author Jose P. Bonatti \n
  !> Daniel M. Lamosa (last revision) \n
  !> Barbara A. G. P. Yamada (last revision)
  !> Eduardo Khamis 08-04-2019
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  subroutine createSpectralRep(self)

    class(legendre) :: self
    integer :: k
    integer :: l          !<
    integer :: m          !<
    integer :: n          !<
    real(kind=p_r8) :: am !<
    real(kind=p_r8) :: an !<
    real(kind=p_r8) :: sn !<


    allocate(self%la0(self%Mend1,self%Mend1))
    allocate(self%la1(self%Mend1,self%Mend2))
    allocate(self%eps(self%mnwv1))
    allocate(self%e0(self%mnwv1))
    allocate(self%e1(self%mnwv0))
    allocate(self%snnp(self%mnwv2))

    l = 0
    k = 0
    do n=1, self%mend1
      l = l + 1
      self%la1(n,1) = l

      self%eps(l) = 0.0_p_r8
      self%e0(l)  = 0.0_p_r8
      self%e1(l)  = emRad/REAL(l,p_r8)

      do m=1, self%mend2 - n
        k = k + 1
        self%la0(m,n) = k

        sn=REAL(m+n-2,p_r8)
        if(sn.ne.0.0_p_r8)then
           sn = -EMRad2/(sn*(sn+1.0_p_r8))
           self%snnp(2*l-1) = sn
           self%snnp(2*l  ) = sn
        else
           self%snnp(2*l-1) = 0.0_p_r8
           self%snnp(2*l  ) = 0.0_p_r8
        endif

        if(n.ge.2)then
           an=n+m-2
           am=m-1
           self%e1(k) = emRad*am/(an+an*an)
        endif
       
      end do
    end do
    self%e1(1) = 0.0_p_r8

    l = self%mend1
    do n=2, self%mend2
      do m=1, self%mend1 + 2 - n
        l = l + 1
        self%la1(m,n) = l

        am = m - 1
        an = m + n - 2
        self%eps(l) = sqrt((an * an - am * am) / (4.0_p_r8 * an * an - 1.0_p_r8))
        self%e0(l)  = emRad*self%eps(l)/REAL(n+m-2,p_r8)

      end do
    end do

  end subroutine createSpectralRep

  ! ---------------------------------------------------------------------------
  !> @brief createGaussRep
  !!
  !! @details 
  !!
  !! @author Jose P. Bonatti \n
  !> Daniel M. Lamosa (last revision) \n
  !> Barbara A. G. P. Yamada (last revision)
  !> Eduardo Khamis 09-04-2019
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  subroutine createGaussRep(self)
    class(legendre) :: self
  
    allocate(self%colRad(self%yMaxHf))
    allocate(self%rCs2(self%yMaxHf))
    allocate(self%wgt(self%yMaxHf))
    allocate(self%gLats(self%yMax))

    call gaussianLatitudes(self)
    
  end subroutine createGaussRep

  subroutine CreateLegPols(self)
     class(legendre) :: self

     allocate(self%legS2F   (self%mnwv0, self%yMaxHf))
     allocate(self%legExtS2F(self%mnwv1, self%yMaxHf))
     allocate(self%legDerS2F(self%mnwv0, self%yMaxHf))
     allocate(self%legDerNS (self%mnwv0, self%yMaxHf))
     allocate(self%legDerEW (self%mnwv0, self%yMaxHf))    

     call legPols(self)

  endsubroutine

  ! ---------------------------------------------------------------------------
  !> @brief gLats
  !!
  !! @details Calculates Gaussian Latitudes and Gaussian Weights 
  !! for Use in Grid-Spectral and Spectral-Grid Transforms
  !!
  !! @author Jose P. Bonatti \n
  !> Daniel M. Lamosa (last revision) \n
  !> Barbara A. G. P. Yamada (last revision)
  !> Eduardo Khamis     09-04-2019
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  !
  !
  subroutine gaussianLatitudes(self)
  !     glats: calculates gaussian latitudes and 
  !            gaussian weights for use in grid-spectral 
  !            and spectral-grid transforms.
  ! 
  !      glats calls the subroutine poly
  ! 
  !     argument(dimensions)        description
  ! 
  !     colrad(Jmaxhf)      output: co-latitudes for gaussian
  !                                 latitudes in one hemisphere.
  !     wgt(Jmaxhf)         output: gaussian weights.
  !     rcs2(Jmaxhf)        output: 1.0/cos(gaussian latitude)**2
  !     Jmaxhf               input: number of gaussian latitudes
  !                                 in one hemisphere.

    type(legendre), intent(inout) :: self
  
    real(kind=p_r8) :: epsil, scal, dgColIn
  
    integer         :: j       !< Loop iterator
    real(kind=p_r8) :: gCol    !<
    real(kind=p_r8) :: dgCol   !<
    real(kind=p_r8) :: p2      !<
    real(kind=p_r8) :: p1      !<
    real(kind=p_r8) :: rad     !<


    epsil   = epsilon(1.0_p_r8) * 100.0_p_r8   
    rad     = 45.0_p_r8 / atan(1.0_p_r8) 
    scal    = 2.0_p_r8 / (real(self%yMax, p_r8) * real(self%yMax, p_r8))
    dgColIn = atan(1.0_p_r8) / real(self%yMax, p_r8) 

    gCol = 0.0_p_r8
    do j=1, self%yMaxHf

      dgCol = dgColIn

      do
        call legendrePolynomial(self%yMax, gCol, p2)
        do
          p1   = p2
          gCol = gCol + dgCol
          call legendrePolynomial(self%yMax, gCol, p2)
          if(sign(1.0_p_r8, p1) /= sign(1.0_p_r8, p2)) exit
        end do
        if(dgCol <= epsil) exit
          gCol  = gCol  - dgCol
          dgCol = dgCol * 0.25_p_r8
       end do

       self%colRad(j) = gCol
       self%gLats(j)  = 90.0_p_r8 - rad * gCol
       self%gLats(self%yMax - j + 1) = -self%gLats(j)

       call legendrePolynomial(self%yMax-1, gCol, p1)

       self%wgt(j)  = scal * (1.0_p_r8 - cos(gCol) * cos(gCol)) / (p1 * p1)
       self%rCs2(j) = 1.0_p_r8 / (sin(gCol) * sin(gCol))
    end do

  end subroutine gaussianLatitudes

  ! ---------------------------------------------------------------------------
  !> @brief legendrePolynomial
  !!
  !! @details Calculates the Value of the Ordinary Legendre 
  !! Function of Given Order at a Specified Colatitude.  
  !! Used to Determine Gaussian Latitudes.
  !!
  !! @author Jose P. Bonatti \n
  !> Daniel M. Lamosa (last revision) \n
  !> Barbara A. G. P. Yamada (last revision)
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  subroutine legendrePolynomial(n, rad, pln)
  ! legendrePolynomial:
  !        calculates the value of the ordinary legendre function
  !        of given order at a specified latitude.  used to
  !        determine gaussian latitudes.
  !
  !***********************************************************************
  !
  ! This routine is called by the subroutine glats.
  !
  ! this routine calls no subroutines.
  !
  !***********************************************************************
  !
  ! argument(dimensions)                       description
  !
  !             n                   input : order of the ordinary legendre
  !                                         function whose value is to be
  !                                         calculated. set in routine
  !                                         "glats".
  !             rad                 input : colatitude (in radians) at
  !                                         which the value of the ordinar
  !                                         legendre function is to be
  !                                         calculated. set in routine
  !                                         "glats".
  !             pln                output : value of the ordinary legendre
  !                                         function of order "n" at
  !                                         colatitude "rad".
  !
  !***********************************************************************
  !

    integer,         intent(in)  :: n   !<
    real(kind=p_r8), intent(in)  :: rad !<
    real(kind=p_r8), intent(out) :: pln !<

    integer         :: i  !< Loop iterator
    real(kind=p_r8) :: x  !<
    real(kind=p_r8) :: y1 !<
    real(kind=p_r8) :: y2 !<
    real(kind=p_r8) :: y3 !<
    real(kind=p_r8) :: g  !<

    x  = cos(rad)
    y1 = 1.0_p_r8
    y2 = x
    do i=2, n
      g  = x * y2
      y3 = g - y1 + g - (g - y1) / real(i, p_r8)
      y1 = y2
      y2 = y3
    end do
    pln = y3

  end subroutine legendrePolynomial


  ! ---------------------------------------------------------------------------
  !> @brief legPols
  !!
  !! @details
  !!
  !!
  !! @author Jose P. Bonatti \n
  !> Daniel M. Lamosa (last revision) \n
  !> Barbara A. G. P. Yamada (last revision)
  !> Eduardo Khamis 19-04-2019 
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  subroutine legPols(self)
    type(legendre), intent(inout) :: self

    integer :: j  !<
    integer :: l  !<
    integer :: nn !<
    integer :: mm !<
    integer :: mn !<
    integer :: lx !<
    
    real(kind=p_r8) :: pln(self%mnwv1)

    do j=1, self%yMaxHf
      
      call pln2(self, pln, j)
      l = 0
      do nn=1, self%mend1
        do mm=1, self%mend2 - nn
          l  = l + 1
          lx = self%la1(mm,nn)
          self%legS2F(l,j) = pln(lx) ! just half part
        end do
      end do

      do mn=1, self%mnwv1
        self%legExtS2F(mn,j) = pln(mn)
      end do

      call plnder(self, j, pln,        &
                  self%legDerS2F(:,j), &
                  self%legDerNS (:,j), &
                  self%legDerEW (:,j)  &
                  )

    end do

  end subroutine legPols


  subroutine pln2 (self, sln, plat)

    !      pln2: calculates the associated legendre functions
    !            at one specified latitude.
    !
    !      plat                input: gaussian latitude index. set
    !                                 by calling routine.
    !      sln (Mnwv1)         output: values of associated legendre
    !                                  functions at one gaussian
    !                                  latitude specified by the
    !                                  argument "lat".
    ! uses from self:
    !      colrad(Jmaxhf)       colatitudes of gaussian grid
    !                           (in radians). calculated
    !                           in routine "glats".
    !      eps                  factor that appears in recusio
    !                           formula of a.l.f.
    !      Mend                 triangular truncation wave number
    !      Mnwv1                number of elements
    !      la1(Mend1,Mend1+1)   numbering array of pln1
    !-------------------------------------------------------!

    type(legendre), intent(in) :: self

    real(kind=p_r8), intent(out) :: sln(self%mnwv1)

    integer, intent(in) :: plat

    integer :: mm, nn, lx, ly, lz

    real(kind=p_r8) :: colr, sinlat, coslat, prod

    logical, SAVE :: first = .true.

    real(kind=p_r8), allocatable, dimension(:), SAVE :: x, y
    real(kind=p_r8), SAVE :: rthf

    if (first) then
       allocate (x(self%mend1))
       allocate (y(self%mend1))
       first=.false.
       do mm=1,self%mend1
          x(mm)=sqrt(2.0_p_r8*mm+1.0_p_r8)
          y(mm)=sqrt(1.0_p_r8+0.5_p_r8/real(mm,p_r8))
       end do
       rthf=sqrt(0.5_p_r8)
   endif
    colr=self%colrad(plat)
    sinlat=cos(colr)
    coslat=sin(colr)
    prod=1.0_p_r8
    do mm=1,self%mend1
       sln(mm)=rthf*prod
       !     line below should only be used where exponent range is limted
       !     if(prod < flim) prod=0.0_p_r8
       prod=prod*coslat*y(mm)
    end do

    do mm=1,self%mend1
       sln(mm+self%mend1)=x(mm)*sinlat*sln(mm)
    end do
    do nn=3,self%mend2
       do mm=1,self%mend1+2-nn
          lx=self%la1(mm,nn)
          ly=self%la1(mm,nn-1)
          lz=self%la1(mm,nn-2)
          sln(lx)=(sinlat*sln(ly)-self%eps(ly)*sln(lz))/self%eps(lx)
       end do
    end do

  end subroutine pln2


  subroutine plnder(self, plat, pln, dpln, der, plnwcs)

    !     calculates zonal and meridional pseudo-derivatives as
    !     well as laplacians of the associated legendre functions.
    !
    !     argument(dimensions)             description
    ! 
    !     pln   (Mnwv1)            input : associated legendre function
    !                                      values at gaussian latitude
    !                             output : pln(l,n)=
    !                                      pln(l,n)*(l+n-2)*(l+n-1)/a**2.
    !                                      ("a" denotes earth's radius)
    !                                      used in calculating the
    !                                      laplacian of global fields
    !                                      in spectral form.
    !     dpln  (Mnwv0)           output : derivatives of
    !                                      associated legendre functions
    !                                      at gaussian latitude
    !     der   (Mnwv0)           output : cosine-weighted derivatives of
    !                                      associated legendre functions
    !                                      at gaussian latitude
    !     plnwcs(Mnwv0)           output : plnwcs(l,n)=
    !                                      pln(l,n)*(l-1)/sin(latitude)**2.
    !
    !uses from self:
    !
    !     rcs2l                    input : 1.0/sin(latitude)**2 at
    !                                      gaussian latitude
    !     wgtl                     input : gaussian weight, at gaussian
    !                                      latitude
    !     eps   (Mnwv1)            input : array of constants used to
    !                                      calculate "der" from "pln".
    !                                      computed in routine "epslon".
    !     la1(Mend1,Mend2)         input : numbering array for pln    
    type(legendre),  intent(in   ) :: self
    integer,         intent(in   ) :: plat
    real(kind=p_r8), intent(inout) :: pln(self%mnwv1)
    real(kind=p_r8), intent(  out) :: dpln(self%mnwv0)
    real(kind=p_r8), intent(  out) :: der(self%mnwv0)
    real(kind=p_r8), intent(  out) :: plnwcs(self%mnwv0)

    integer :: n, l, nn, mm, mn, lm, l0, lp

    real(kind=p_r8) :: raa, wcsa
    real(kind=p_r8) :: x(self%mnwv1)

    logical, SAVE :: first=.true.

    real(kind=p_r8), allocatable, SAVE :: an(:)
    ! 
    !     compute pln derivatives
    ! 
    if (first) then
       allocate (an(self%mend2))
       do n=1,self%mend2
          an(n)=real(n-1,p_r8)
       end do
       first=.false.
    end if
    raa=self%wgt(plat)*emRad12
    wcsa=self%rCs2(plat)*self%wgt(plat)*emRad1
    l=0
    do mm=1,self%mend1
       l=l+1
       x(l)=an(mm)
    end do
    do nn=2,self%mend2
       do mm=1,self%mend1+2-nn
          l=l+1
          x(l)=an(mm+nn-1)
       end do
    end do
    l=self%mend1
    do nn=2,self%mend1
       do mm=1,self%mend2-nn
          l=l+1
          lm=self%la1(mm,nn-1)
          l0=self%la1(mm,nn)
          lp=self%la1(mm,nn+1)
          der(l)=x(lp)*self%eps(l0)*pln(lm)-x(l0)*self%eps(lp)*pln(lp)
       end do
    end do
    do mm=1,self%mend1
       der(mm)=-x(mm)*self%eps(mm+self%Mend1)*pln(mm+self%Mend1)
    end do
    do mn=1,self%mnwv0
       dpln(mn) = der(mn)
       der(mn)  = wcsa*der(mn)
    end do
    l=0
    do nn=1,self%mend1
       do mm=1,self%mend2-nn
          l=l+1
          l0=self%la1(mm,nn)
          plnwcs(l)=an(mm)*pln(l0)
       end do
    end do
    do mn=1,self%mnwv0
       plnwcs(mn)=wcsa*plnwcs(mn)
    end do
    do nn=1,self%mend1
       do mm=1,self%mend2-nn
          l0=self%la1(mm,nn)
          lp=self%la1(mm,nn+1)
          pln(l0)=x(l0)*x(lp)*raa*pln(l0)
       end do
    end do


  end subroutine plnder

   subroutine SpecTransp_ (self, itype, a, MnWv)
    !     transp: after input, transposes scalar arrays
    !             of spectral coefficients by swapping
    !             the order of the subscripts
    !             representing the degree and order
    !             of the associated legendre functions.
    ! 
    !     argument(dimensions)        description
    ! 
    !     a(:)                input: spectral representation of a
    !                                 global field. Size should be
    !                                 mnwv2 or mnwv3
    !                                 signal=+1 diagonalwise storage
    !                                 signal=-1 coluMnwise   storage
    !                         output: spectral representation of a
    !                                 global field.
    !                                 signal=+1 coluMnwise   storage
    !                                 signal=-1 diagonalwise storage

     class(legendre),  intent(in   ) :: self
     integer,          intent(in   ) :: itype
     integer,          intent(in   ) :: MnWv
     real (kind=p_r8), intent(inout) :: a(MnWv)

     integer :: l, lx, n, m, mn
     integer :: mMax
     integer, pointer :: idx(:,:) => null()
     real (kind=p_r8), allocatable :: qwork(:)
   

     if(MnWv == self%mnwv2)then
        idx  => self%la0
        mMax = self%Mend2

     elseif(MnWv == self%mnwv3)then
        idx  => self%la1
        mMax = self%Mend2 + 1

     else
        write(*,'(A,1x,I6)')'>> error: wrong size of spectral coefficents:',MnWv
        write(*,'(A,1x,I6,1x,A,1x,I6)')'Size should be:',self%MnWv2,'or',self%MnWv3
        return
     endif

     allocate(qwork(MnWv))
     qwork=0.0_p_r8

     if (itype == +1) then

        l=0
        do m=1,self%mend1
          do n=1,mMax-m
            l=l+1
            lx=idx(m,n)
            qwork( 2*l-1 ) = a( 2*lx-1 )
            qwork( 2*l   ) = a( 2*lx   )
           end do
        end do
        do mn=1,mnwv
          a(mn) = qwork(mn)
        end do

     else if(itype == -1)then

        l=0
        do m=1,self%mend1
          do n=1,mMax-m
            l=l+1
            lx=idx(m,n)
            qwork( 2*lx-1 ) = a( 2*l-1 )
            qwork( 2*lx   ) = a( 2*l   )
          end do
        end do

        do mn=1,mnwv
          a(mn) = qwork(mn)
        end do
     end if

     deallocate(qwork)
   
   end subroutine SpecTransp_
  
   subroutine Fourier2SpecCoef_ (self, fp, fm, lat, MnWv, fln)
    
     !     calculates spectral representations of
     !     global fields from fourier representations
     !     of symmetric and anti-symmetric portions
     !     of global fields.
     ! 
     !     argument(dimensions)        description
     ! 
     !     fp(Imax+2)           input: fourier representation of
     !                                 symmetric portion of a
     !                                 global field at one gaussian
     !                                 latitude.
     !     fm(Imax+2)           input: fourier representation of
     !                                 anti-symmetric portion of a
     !                                 global field at one gaussian
     !                                 latitude.
     !     fln(Mnwv2) or        input: spectral representation of
     !     fln(MnWv3)                  (the laplacian of) a global
     !                                 field. includes contributions
     !                                 from gaussian latitudes up to
     !                                 but not including current
     !                                 iteration of gaussian loop
     !                                 in calling routine.
     !                         output: spectral representation of
     !                                 (the laplacian of) a global
     !                                 field. includes contributions
     !                                 from gaussian latitudes up to
     !                                 and including current
     !                                 iteration of gaussian loop
     !                                 in calling routine.
     !     lat                  input: current index of gaussian
     !                                 loop in calling routine.
   
     class(legendre),   intent (in   ) :: self
     real (kind=p_r8),  intent (in   ) :: fp(:)
     real (kind=p_r8),  intent (in   ) :: fm(:)
     integer,           intent (in   ) :: lat
     integer,           intent (in   ) :: MnWv
     real (kind=p_r8),  intent (inout) :: fln(MnWv)
   
     integer :: k, l, nn, mmax, mm, mn
   
     real (kind=p_r8), dimension (self%mnwv2) :: s
   
     l=0
     do nn=1,self%mend1
        mmax=2*(self%mend2-nn)
        if (mod(nn-1,2) == 0) then
           do mm=1,mmax
              l=l+1
              s(l)=fp(mm)
           end do
        else
           do mm=1,mmax
              l=l+1
              s(l)=fm(mm)
           end do
        end if
     end do
     do mn=1,self%mnwv0
        fln(2*mn-1)=fln(2*mn-1)+s(2*mn-1)*self%legS2F(mn,lat)*self%wgt(lat)
        fln(2*mn  )=fln(2*mn  )+s(2*mn  )*self%legS2F(mn,lat)*self%wgt(lat)
     end do

   end subroutine Fourier2SpecCoef_

   subroutine SpecCoef2Fourier_(self, fln, ap, am, lat)
     !     calculates the fourier representation of a field at a
     !     pair of latitudes symmetrically located about the
     !     equator. the calculation is made using the spectral
     !     representation of the field and the values of the
     !     associated legendre functions at that latitude.
     !     it is designed to triangular truncation only and
     !     for scalar fields.
     !
     !     argument(dimensions)            description
     !
     !     fln(Mnwv2)             input: spectral representation of a
     !         MnWv3                     global field.
     !     ap(Imax+2)            output: fourier representation of
     !                                   a global field at the
     !                                   latitude in the northern
     !                                   hemisphere at which the
     !                                   associated legendre functions
     !                                   have been defined.
     !     am(Imax+2)            output: fourier representation of
     !                                   a global field at the
     !                                   latitude in the southern
     !                                   hemisphere at which the
     !                                   associated legendre functions
     !                                   have been defined.
     !     lat                    input: current index of gaussian
     !                                   loop in calling routine.
     character(len=17) :: myname_="**(Spec2Fourier_)**"
   
     class(legendre),   intent (in   ) :: self
     real (kind=p_r8),  intent (in   ) :: fln(:)
     integer,           intent (in   ) :: lat
     real (kind=p_r8),  intent (  out) :: ap(:)
     real (kind=p_r8),  intent (  out) :: am(:)
   
     integer          :: l, k, mm, mn, nn, mstr, mmax1, mend1d
     integer          :: nHarm
     integer          :: MnWv
     integer          :: nMax
     integer          :: mMax
     real (kind=p_r8), allocatable :: s(:)
     real (kind=p_r8), pointer     :: leg(:,:) => null()

     nHarm = size(fln)
     if(nHarm .eq. self%MnWv2)then
        nMax   = self%Mend1
        MnWv   = self%MnWv0
        mend1d = 2 * self%mend1
        mmax1  = 2 * self%mend
        leg => self%legS2F
        allocate(s(self%MnWv2))
     elseif (nHarm .eq. self%MnWv3)then
        nMax   = self%Mend2
        MnWv   = self%MnWv1
        mend1d = 2 * self%mend1
        mmax1  = 2 * self%mend1
        leg => self%legExtS2F
        allocate(s(self%MnWv3))
     else
        write(unit=p_nferr, fmt='(2(A,1x),2(I6,1x,A,1x),I6)')&
              trim(myname_),&
              'unknow size of spectral field ',mnwv,&
              'should be',self%mnwv2,'or',self%mnwv3
         stop 3022
     endif
   
     l      = mend1d + mmax1
     do mn=1,MnWv
        s(2*mn-1) = leg(mn,lat)*fln(2*mn-1)
        s(2*mn  ) = leg(mn,lat)*fln(2*mn  )
     end do

     do nn=3,nMax
        mMax = 2 * (nMax+1-nn)
        if (mod(nn-1,2) == 0) then
           mstr=0
        else
           mstr=mend1d
        end if
        do mm=1,mMax
           l=l+1
           s(mm+mstr)=s(mm+mstr)+s(l)
        end do
     end do
     do mm=1,mend1d
        ap(mm)=s(mm)
        am(mm)=s(mm)
     end do
     do mm=1,mmax1
        ap(mm)=ap(mm)+s(mm+mend1d)
        am(mm)=am(mm)-s(mm+mend1d)
     end do

     deallocate(s)

   end subroutine


   subroutine SymAsy_ (self, a, b)
   
     !     converts the fourier representations of a field at two
     !     parallels at the same latitude in the northern and
     !     southern hemispheres into the fourier representations
     !     of the symmetric and anti-symmetric portions of that
     !     field at the same distance from the equator as the
     !     input latitude circles.
     ! 
     !     argument(dimensions)        description
     ! 
     !     a(Imx,Ldim)          input: fourier representation of one
     !                                 latitude circle of a field
     !                                 from the northern hemisphere
     !                                 at "n" levels in the vertical.
     !                         output: fourier representation of the
     !                                 symmetric portion of a field
     !                                 at the same latitude as the
     !                                 input, at "n" levels in 
     !                                 the vertical.
     !     b(Imx,Ldim)          input: fourier representation of one
     !                                 latitude circle of a field
     !                                 from the southern hemisphere
     !                                 at "n" levels in the vertical.
     !                         output: fourier representation of the
     !                                 anti-symmetric portion of a
     !                                 field at the same latitude as
     !                                 the input, at "n" levels in
     !                                 the vertical.
     !     t(Imx,Ldim)                 temporary storage
     ! 
    
   
     class(legendre),  intent (in   ) :: self
     real (kind=p_r8), intent (inout) :: a(self%xMx)
     real (kind=p_r8), intent (inout) :: b(self%xMx)
   
     integer :: i, k
   
     real (kind=p_r8) :: t(self%xMx)
   
     do i=1,self%xMax ! rever essa dimensao, acredito que deveria se xMx
        t(i)=a(i)
        a(i)=t(i)+b(i)
        b(i)=t(i)-b(i)
     end do
   
   end subroutine SymAsy_

! Eduardo Khamis
  subroutine getxMaxyMax (mend, xMax, yMax, linG)

    implicit none

    integer, intent(in)            :: mend
    integer, intent(out)           :: xMax, yMax
    logical, intent(in), optional  :: linG       ! Flag for linear (T) or quadratic (F) triangular truncation

    logical           :: linearGrid
    integer           :: nx, nm, n2m, n3m, n5m, n2, n3, n5, j, n, check, yfft
    integer, save     :: lfft = 40000

    integer, dimension(:), allocatable, save :: xfft

    if (present(linG)) then
       linearGrid = linG
    else
       linearGrid = .false.
    end if

    n2m = ceiling(log(real(lfft, p_r8)) / log(2.0_p_r8))
    n3m = ceiling(log(real(lfft, p_r8)) / log(3.0_p_r8))
    n5m = ceiling(log(real(lfft, p_r8)) / log(5.0_p_r8))
    nx  = n2m * (n3m + 1) * (n5m + 1)

    allocate(xfft (nx))
    xfft = 0

    n = 0
    do n2 = 1, n2m
      yfft = (2**n2)
      if (yfft > lfft) exit
      do n3 = 0, n3m
        yfft = (2**n2) * (3**n3)
        if (yfft > lfft) exit
        do n5 = 0, n5m
          yfft = (2**n2) * (3**n3) * (5**n5)
          if (yfft > lfft) exit
          n = n + 1
          xfft(n) = yfft
        end do
      end do
    end do
    nm = n

    n = 0
    do
      check = 0
      n = n + 1
      do j = 1, nm - 1
        if (xfft(j) > xfft(j + 1)) then
          yfft = xfft(j)
          xfft(j) = xfft(j + 1)
          xfft(j + 1) = yfft
          check = 1
        end if
      end do
      if (check == 0) exit
    end do

    if (linearGrid) then
       yfft=2
    else
       yfft=3
    end if
    xMax = yfft * mend + 1
    do n = 1, nm
      if (xfft(n) >= xMax) then
        xMax = xfft(n)
        exit
      end if
    end do
    yMax = xMax / 2
    if (mod(yMax, 2) /= 0) yMax = yMax + 1

    deallocate(xfft)

  end subroutine getxMaxyMax


  function GetSize(self,sname)result(sval)
     class(legendre)  :: self
     character(len=*) :: sname
     integer          :: sval
     !
     ! Colocar aqui um lower case em sname!
     ! sname = lower(sname)
     !
     select case(trim(sname))
        case('mend')
           sval = self%Mend
        case('mend1')
           sval = self%Mend1
        case('mend2')
           sval = self%Mend2
        case('mnwv0')
           sval = self%MnWv0
        case('mnwv1')
           sval = self%MnWv1
        case('mnwv2')
           sval = self%MnWv2
        case('mnwv3')
           sval = self%MnWv3
        case('xmax')
           sval = self%xMax
        case('ymax')
           sval = self%yMax
        case('ymaxhf')
           sval = self%yMaxHf
        case default
           write(*,'(A,1x,A)')'unknown size var name:',trim(sname)
           sval = -1
     end select
  end function

  subroutine DestroySpectralRep( self )
     class(legendre), intent(inout) :: self

     deallocate( self%la0 )
     deallocate( self%la1 )
     deallocate( self%eps )
     deallocate( self%e0  )
     deallocate( self%e1  )
     deallocate( self%snnp)

  end subroutine

  subroutine DestroyGaussRep( self )
     class(legendre), intent(inout) :: self

     deallocate ( self%colRad )
     deallocate ( self%rCs2   )
     deallocate ( self%wgt    )
     deallocate ( self%gLats  )

  end subroutine

  subroutine DestroyLegPols( self )
     class(legendre), intent(inout) :: self

     deallocate ( self%legS2F    )
     deallocate ( self%legDerS2F )
     deallocate ( self%legExtS2F )
     deallocate ( self%legDerNS  )
     deallocate ( self%legDerEW  )

  end subroutine


  subroutine print_size(label)
  
  !Cray-JNT!  
  ! Esta subrotina reporta o Virtual Memory High Water Mark (VmHWM)  
  ! presente no arquivo /proc/self/status do sistema Linux.
  
  !
  character(len=256),intent(in) :: label

  integer           :: iunit  
  character(len=80) :: linha
  
    iunit = 7
  
    open(iunit,file="/proc/self/status")
  
    do
  
      read(iunit,'(a)',END=999) linha  
      if (INDEX(linha,'VmHWM').eq.1) then  
         print *,trim(label), linha  
      endif
  
    enddo
  
  999 close(iunit)
  
  end subroutine print_size

end module LegendreTransform
