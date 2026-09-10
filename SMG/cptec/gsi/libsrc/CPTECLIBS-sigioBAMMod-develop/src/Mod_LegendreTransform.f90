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
module Mod_LegendreTransform
  use ModConstants, only: p_r8 => r8
  use ModConstants, only: emRad
  use ModConstants, only: emRad1
  use ModConstants, only: emRad12
  implicit none
  private

  integer, parameter :: p_nferr = 0


  interface splitTrans
     module procedure splitTrans2D, splitTrans3D
  end interface

  interface Four2Spec
     module procedure Four2Spec1D, Four2Spec2D
  end interface

!  interface spec2Four
!     module procedure spec2Four1D, spec2Four2D
!  end interface
!
!  interface four2Spec
!     module procedure four2Spec1D, four2Spec2D
!  end interface

  type, public :: legendre
     private

     logical :: created

     integer :: Mend
     integer :: Mend1
     integer :: Mend2
     integer :: MnWv0
     integer :: MnWv1
     integer :: MnWv2
     integer :: MnWv3
     integer :: xMax
     integer :: yMax
     integer :: yMaxHf
     
     integer, pointer, dimension(:)           :: lenDiag         => null()
     integer, pointer, dimension(:)           :: lenDiagExt      => null()
     integer, pointer, dimension(:)           :: lastPrevDiag    => null()
     integer, pointer, dimension(:)           :: lastPrevDiagExt => null()
     integer, pointer, dimension(:,:)         :: la0       => null()
     integer, pointer, dimension(:,:)         :: la1       => null()
     real(kind=p_r8), pointer, dimension(:)   :: e0        => null()
     real(kind=p_r8), pointer, dimension(:)   :: e1        => null()
     real(kind=p_r8), pointer, dimension(:)   :: eps       => null()
     real(kind=p_r8), pointer, dimension(:)   :: colRad    => null()
     real(kind=p_r8), pointer, dimension(:)   :: rCs2      => null()
     real(kind=p_r8), pointer, dimension(:)   :: wgt       => null()
     real(kind=p_r8), pointer, dimension(:)   :: gLats     => null()
     real(kind=p_r8), pointer, dimension(:,:) :: legS2F    => null()
     real(kind=p_r8), pointer, dimension(:,:) :: legExtS2F => null()
     real(kind=p_r8), pointer, dimension(:,:) :: legDerS2F => null()
     real(kind=p_r8), pointer, dimension(:,:) :: legF2S    => null()
     real(kind=p_r8), pointer, dimension(:,:) :: legDerNS  => null()
     real(kind=p_r8), pointer, dimension(:,:) :: legDerEW  => null()
     contains
        procedure, public  :: init => initialization_
        procedure, public  :: destroy => destroylegendre
        procedure, public  :: transp => specTransp_

!        procedure, public  :: Four2Spec
        procedure, private :: Four2Spec1D, Four2Spec2D
        generic,   public  :: Four2Spec => Four2Spec1D, Four2Spec2D

        procedure, private :: Spec2Four1D, Spec2Four2D
        generic,   public  :: Spec2Four => Spec2Four1D, Spec2Four2D

        procedure, public  :: DivgVortToUV

!        generic,   public  :: splitTrans => splitTrans1D, splitTrans2D
!        procedure, private :: splitTrans1D, splitTrans2D


        procedure, public  :: sumSpec
        procedure, public  :: GetSize

        procedure, private :: createSpectralRep
        procedure, private :: createGaussRep
        procedure, private :: createLegTrans

        procedure, private :: DestroySpectralRep
        procedure, private :: DestroyGaussRep
        procedure, private :: DestroyLegTrans

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

     self%MnWv2 = (Mend + 1) * (Mend + 2)
     self%MnWv0 = self%MnWv2 / 2
     self%MnWv3 = self%MnWv2 + 2 * (self%Mend1)
     self%MnWv1 = self%MnWv3 / 2

     call getxMaxyMax(Mend, self%xMax, self%yMax)
     self%yMaxHf = self%yMax/2

     self%created = .false.
     call self%createSpectralRep( )
     call self%createGaussRep( )
     call self%createLegTrans( )

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
     call self%destroyLegTrans( )     

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

    allocate(self%la0(self%Mend1,self%Mend1))
    allocate(self%la1(self%Mend1,self%Mend2))
    allocate(self%eps(self%mnwv1))
    allocate(self%e0(self%mnwv1))
    allocate(self%e1(self%mnwv0))

!    l = 0
!    do n=1, self%mend1
!      do m=1, self%mend2 - n
!        l = l + 1
!        self%la0(m,n) = l
!      end do
!    end do
!    l = 0
!    do m=1, self%mend1
!      l = l + 1
!      self%la1(m,1) = l
!    end do
!    do n=2, self%mend2
!      do m=1, self%mend1 + 2 - n
!        l = l + 1
!        self%la1(m,n) = l
!      end do
!    end do
!
!    do l=1, self%mend1
!      self%eps(l) = 0.0_p_r8
!    end do
!    l = self%mend1
!    do n=2, self%mend2
!      do m=1, self%mend1 + 2 - n
!        l = l + 1
!        am = m - 1
!        an = m + n - 2
!        self%eps(l) = sqrt((an * an - am * am) / (4.0_p_r8 * an * an - 1.0_p_r8))
!       end do
!    end do

!    do l=1, self%mend1
!      self%eps(l) = 0.0_p_r8
!      self%e0(l)  = 0.0_p_r8
!      self%e1(l)  = emRad/REAL(l,p_r8)
!    end do
!    self%e1(1) = 0.0_p_r8
!
!    l = self%mend1
!    do n=2, self%mend2
!      do m=1, self%mend1 + 2 - n
!        l = l + 1
!        am = m - 1
!        an = m + n - 2
!        self%eps(l) = sqrt((an * an - am * am) / (4.0_p_r8 * an * an - 1.0_p_r8))
!        self%e0(l)  = emRad*self%Eps(l)/REAL(n+m-2,p_r8)
!       end do
!    end do

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
  subroutine gaussianLatitudes(self)
!    note: pgi failed to compile epsil, scal and dgColIn as parameter
!    real(kind=p_r8), parameter :: epsil   = epsilon(1.0_p_r8) * 100.0_p_r8    !<
!    real(kind=p_r8), parameter :: scal    = 2.0_p_r8 / (real(yMax, p_r8)&
!                                          * real(yMax, p_r8))                 !< scale
!    real(kind=p_r8), parameter :: dgColIn = atan(1.0_p_r8) / real(yMax, p_r8) !<
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
  subroutine legendrePolynomial(n, colatitude, pln)
    integer,         intent(in)  :: n          !<
    real(kind=p_r8), intent(in)  :: colatitude !<
    real(kind=p_r8), intent(out) :: pln        !<

    integer         :: i  !< Loop iterator
    real(kind=p_r8) :: x  !<
    real(kind=p_r8) :: y1 !<
    real(kind=p_r8) :: y2 !<
    real(kind=p_r8) :: y3 !<
    real(kind=p_r8) :: g  !<

    x  = cos(colatitude)
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
  !> @brief createLegTrans
  !!
  !! @details
  !!
  !!
  !! @author Jose P. Bonatti \n
  !> Daniel M. Lamosa (last revision) \n
  !> Barbara A. G. P. Yamada (last revision)
  !> Eduardo Khamis 15-04-2019
  !!
  !! @date mar/2018
  ! ---------------------------------------------------------------------------
  subroutine createLegTrans(self)
    class(legendre) :: self

    character(len=20) :: h = "**(createLegTrans)**" !<

    integer :: diag !<

    if(self%created) then
      write(unit=p_nferr, fmt='(2a)') h, ' already created'
      stop
    else
       self%created = .true.
    end if

    ! Associated Legendre Functions
    allocate(self%legS2F   (self%mnwv2, self%yMaxHf))
    allocate(self%legDerS2F(self%mnwv2, self%yMaxHf))
    allocate(self%legExtS2F(self%mnwv3, self%yMaxHf))
    allocate(self%legF2S   (self%mnwv2, self%yMaxHf))
    allocate(self%legDerNS (self%mnwv2, self%yMaxHf))
    allocate(self%legDerEW (self%mnwv2, self%yMaxHf))

    call legPols(self)

    ! diagonal length
    allocate(self%lenDiag(self%mend1))
    allocate(self%lenDiagExt(self%mend2))

    do diag=1, self%mend1
      self%lenDiag(diag) = 2 * (self%mend1 + 1 - diag)
    end do
    self%lenDiagExt(1) = 2 * self%mend1
    do diag=2, self%mend2
      self%lenDiagExt(diag) = 2 * (self%mend1 + 2 - diag)
    end do

    ! last element previous diagonal
    allocate(self%lastPrevDiag(self%mend1))
    allocate(self%lastPrevDiagExt(self%mend2))

    do diag=1, self%mend1
      self%lastPrevDiag(diag) = (diag-1) * (2 * self%mend1 + 2 - diag)
    end do
    self%lastPrevDiagExt(1) = 0
    do diag=2, self%mend2
      self%lastPrevDiagExt(diag) = (diag-1) * (2 * self%mend1 + 4 - diag) - 2
    end do
    
  end subroutine createLegTrans

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
    real(kind=p_r8) :: dpln(self%mnwv0)
    real(kind=p_r8) :: der(self%mnwv0)
    real(kind=p_r8) :: plnwcs(self%mnwv0)
!    character(len=80) :: linha

    do j=1, self%yMaxHf
      call pln2(self, pln, j)
      l = 0
      do nn=1, self%mend1
        do mm=1, self%mend2 - nn
          l  = l + 1
          lx = self%la1(mm,nn)
          self%legS2F(2*l-1,j) = pln(lx)
          self%legS2F(2*l,  j) = pln(lx)
        end do
      end do
      do mn=1, self%mnwv2
        self%legF2S(mn,j) = self%legS2F(mn,j) * self%wgt(j)
       end do
      do mn=1, self%mnwv1
        self%legExtS2F(2*mn-1,j) = pln(mn)
        self%legExtS2F(2*mn,  j) = pln(mn)
      end do
      call plnder(self, j, pln, dpln, der, plnwcs)
      do mn=1, self%mnwv0
        self%legDerS2F(2*mn-1,j) = dpln(mn)
        self%legDerS2F(2*mn,  j) = dpln(mn)

        self%legDerNS(2*mn-1,j)  = der(mn)
        self%legDerNS(2*mn,  j)  = der(mn)

        self%legDerEW(2*mn-1,j)  = plnwcs(mn)
        self%legDerEW(2*mn,  j)  = plnwcs(mn)
      end do
!  open(7,file="/proc/self/status")
!
!  do
!
!    read(7,'(a)',END=990) linha
!
!    if (INDEX(linha,'VmHWM').eq.1) then
!
!       print *,'legPols', j, linha
!
!    endif
!
!  enddo
!
!990 close(7)
    end do

  end subroutine legPols


! Eduardo Khamis
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

!    deallocate(x)
!    deallocate(y)

  end subroutine pln2


! Eduardo Khamis
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
!
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

!    deallocate(an)

  end subroutine plnder

!   subroutine transScalar(self, itype, a)
!      class(legendre) :: self
!      integer, intent(in) :: itype
!      real (kind=p_r8), intent (inout) :: a(self%mnWv2) ! MnWv2
! 
!      call SpecTransp_(self, itype, a)
! 
!   end subroutine
!
!   subroutine transVector(self,itype,a)
!      class(legendre) :: self
!      integer, intent(in) :: itype
!      real (kind=p_r8), intent (inout) :: a(self%MnWv3) ! MnWv3
! 
!      call SpecTransp_(self, itype, a)
! 
!   end subroutine


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
  

! Eduardo Khamis
!  subroutine sumSpec (nMax, mMax, mnwv, xMax, yMax, yMaxHf, zMax, &
!                      spec, leg, four, dLength, lastPrev)
  subroutine sumSpec (self, spec, four, flag)

    ! fourier representation from spectral representation

    class(legendre), intent(in   ) :: self
    real(kind=p_r8), intent(in   ) :: spec(:)  ! spectral field
    real(kind=p_r8), intent(  out) :: four(:,:)! full fourier field
    logical,         intent(in   ) :: flag
    character(len=17) :: h="**(sumSpec)**"
    
    integer :: mnwv
    integer :: nMax
    integer :: mMax
    integer :: xMax, yMax

    integer :: j, ele, diag

    real(kind=p_r8), allocatable :: oddDiag(:) !(2*mMax,yMaxHf,zMax)
    real(kind=p_r8), allocatable :: evenDiag(:) !(2*mMax,yMaxHf,zMax)

    real(kind=p_r8), pointer :: leg(:,:) => null()
    integer,         pointer :: dLength(:) => null()
    integer,         pointer :: lastPrev(:) => null()

    mnwv = size(spec)
    mMax = self%Mend1
    xMax = size(four,1)
    yMax = size(four,2)

    if(mnwv .eq. self%mnwv2)then

        nMax     =  self%Mend1
        if(flag)then
           leg      => self%legDerS2F
        else
           leg      => self%legS2F
        endif
        dlength  => self%lenDiag
        lastPrev => self%lastPrevDiag

    elseif(mnwv .eq. self%mnwv3)then

        nMax     =  self%Mend2
        leg      => self%legExtS2F
        dlength  => self%lenDiagExt
        lastPrev => self%lastPrevDiagExt

    else
       write(unit=p_nferr, fmt='(2(A,1x),2(I6,1x,A,1x),I6)')&
             trim(h),&
             'unknow size of spectral field ',mnwv,&
             'should be',self%mnwv2,'or',self%mnwv3

    endif

    allocate(oddDiag(2*self%Mend1))
    allocate(evenDiag(2*self%Mend1))

    do j=1,self%yMaxHf

       oddDiag  = 0.0_p_r8
       evenDiag = 0.0_p_r8

    ! sum odd diagonals (n+m even)

       do diag=1,nMax,2
          do ele=1,dLength(diag)
             oddDiag(ele)=oddDiag(ele)+ &
                  leg(ele+lastPrev(diag),j)*spec(ele+lastPrev(diag))
          end do
       end do


    ! sum even diagonals (n+m odd)

       do diag=2,nMax,2
          do ele=1,dLength(diag)
             evenDiag(ele)=evenDiag(ele)+ &
                  leg(ele+lastPrev(diag),j)*spec(ele+lastPrev(diag))
          end do
       end do


    ! use even-odd properties

       do ele=1,2*self%Mend1
          four(ele,j)        = oddDiag(ele)+evenDiag(ele)
          four(ele,yMax+1-j) = oddDiag(ele)-evenDiag(ele)
       end do

    end do

    four(2*self%Mend1+1:xMax,:) = 0.0_p_r8

    deallocate(oddDiag)
    deallocate(evenDiag)

  end subroutine sumSpec

! Eduardo Khamis
  subroutine sumFour (self, four, spec)

    ! spectral representation from fourier representation

    type(legendre),  intent(in   ) :: self
    real(kind=p_r8), intent(in   ) :: four(:,:) ! full fourier field
    real(kind=p_r8), intent(  out) :: spec(:)   ! spectral field

    integer :: mnwv
    integer :: nMax

    integer :: j, jj, ele, diag

    real(kind=p_r8), allocatable:: fourEven(:)
    real(kind=p_r8), allocatable:: fourOdd(:)

    ! initialize result

    mnwv = size(spec)
    if(mnwv .eq. self%mnwv2)then
       nMax = self%Mend1
    else if(mnwv .eq. self%mnwv3)then
       nMax = self%Mend2
    endif

    spec=0.0_p_r8

    allocate(fourEven(2*self%Mend1))
    allocate(fourOdd(2*self%Mend1))

    do j=1,self%yMaxHf

       fourEven = 0.0_p_r8
       fourOdd  = 0.0_p_r8

       jj=self%yMax-j+1

       do ele=1,2*self%Mend1
          fourEven(ele) = four(ele,j)+four(ele,jj)
          fourOdd (ele) = four(ele,j)-four(ele,jj)
       end do

!  sum odd diagonals (n+m even)

       do diag=1,nMax,2
          do ele=1,self%lenDiag(diag)
             spec(ele+self%lastPrevDiag(diag)) = spec(ele+self%lastPrevDiag(diag))+ &
                  fourEven(ele)*self%legF2S(ele+self%lastPrevDiag(diag),j)
          end do
       end do

!  sum even diagonals (n+m odd)

       do diag=2,nMax,2
          do ele=1,self%lenDiag(diag)
             spec(ele+self%lastPrevDiag(diag)) = spec(ele+self%lastPrevDiag(diag))+ &
                  fourOdd(ele)*self%legF2S(ele+self%lastPrevDiag(diag),j)
          end do
       end do

    end do

    deallocate(fourEven)
    deallocate(fourOdd)

  end subroutine sumFour


! Eduardo Khamis 16-04-2019
  subroutine spec2Four1D (self, spec, four, der)
    implicit none

    class(legendre),   intent(in   ) :: self
    real(kind=p_r8),   intent(in   ) :: spec(:)
    real(kind=p_r8),   intent(  out) :: four(:,:)
    logical, optional, intent(in   ) :: der

    integer :: s1, f1, f2

    logical :: extended, derivate

    character(len=17) :: h="**(spec2Four1D)**"

    if (.NOT. self%created) then
       write(unit=p_nferr, fmt='(2A)') h, &
             ' Module not created; invoke InitLegTrans prior to this call'
       stop
    end if

    s1=size(spec,1)
    f1=size(four,1); f2=size(four,2)

    if (s1 .eq. self%MnWv2)then
          extended = .false.
    else if(s1 .eq. self%MnWv3)then
       extended = .true.
    else
       write(unit=p_nferr, fmt='(2A,I10)') h, &
       ' wrong first dim of spec : ', s1
       stop
    end if

    if (s1 /= self%mnwv2 .and. s1 /= self%mnwv3) then
       write(unit=p_nferr, fmt='(2A,3I10)') h, &
             ' wrong first dim of spec: ', s1
    end if

    if (f1 < 2*self%mend1) then
       write(unit=p_nferr, fmt='(2A,2I10)') h, &
             ' first dimension of four too small: ', f1, 2*self%mend1
       stop
    end if

    if (f2 /= 2*self%yMaxHf) then
       write(unit=p_nferr, fmt='(2A,I10,A,I10)') h, &
             ' second dimension of four is ', f2, '; should be ', 2*self%yMaxHf
       stop
    end if

    if (present(der)) then
       derivate = der
    else
       derivate = .false.
    end if

    if (derivate .and. extended) then
       write(unit=p_nferr, fmt='(2A)') h, &
             ' derivative cannot be applied to extended gaussian field'
       stop
    end if

    call self%sumSpec( spec, four, derivate)

  end subroutine spec2Four1D
  
! Eduardo Khamis  17-04-2019
  subroutine spec2Four2D (self, spec, four)
  
    class(legendre), intent(in   ) :: self
    real(kind=p_r8), intent(in   ) :: spec(:,:)
    real(kind=p_r8), intent(  out) :: four(:,:,:)

    integer :: specZ
    integer :: fourZ
    integer :: k

    character(len=17) :: h="**(spec2Four2D)**"

    if (.NOT. self%created) then
       write(unit=p_nferr, fmt='(2A)') h, &
             ' Module not created; invoke Init() prior to this call'
       stop
    end if

    
    specZ = size(spec,2)
    fourZ = size(four,3)


    if (specZ /= fourZ) then
       write(unit=p_nferr, fmt= '(2A,2I10)') h, &
             ' vertical layers of spec and four dissagre :', specZ, fourZ
       stop
    end if

    do k=1,specZ
       call self%spec2Four(spec(:,k), four(:,:,k))
    enddo


  end subroutine spec2Four2D

! Eduardo Khamis 17-04-2019
  subroutine four2Spec1D (self, four, spec)
    implicit none
    class(legendre)                :: self
    real(kind=p_r8), intent(in   ) :: four(:,:)
    real(kind=p_r8), intent(  out) :: spec(:)

    integer :: s1, f1, f2

    character(len=17) :: h="**(four2Spec1D)**"

    if (.NOT. self%created) then
       write(unit=p_nferr, fmt='(2A)') h, &
             ' Module not created; invoke ``Init( Mend )'' prior to this call'
       stop
    end if

    s1=size(spec,1)
    f1=size(four,1); f2=size(four,2)

    if (s1 /= self%mnwv2) then
       write(unit=p_nferr, fmt='(2A,I10)') h, &
             ' wrong first dim of spec: ', s1
    end if

    if (f1 < 2*self%mend1) then
       write(unit=p_nferr, fmt='(2A,2I10)') h, &
             ' first dimension of four too small: ', f1, 2*self%mend1
       stop
    end if

    if (f2 /= 2*self%yMaxHf) then
       write(unit=p_nferr, fmt='(2A,I10,A,I10)') h, &
             ' second dimension of four is ', f2, '; should be ', 2*self%yMaxHf
       stop
    end if

    call sumFour (self, four, spec)

  end subroutine four2Spec1D


! Eduardo Khamis 17-04-2019
  subroutine four2Spec2D (self, four, spec)
    implicit none
    class(legendre)                :: self
    real(kind=p_r8), intent(in   ) :: four(:,:,:)
    real(kind=p_r8), intent(  out) :: spec(:,:)

    integer :: specZ
    integer :: fourZ
    integer :: k

    character(len=17) :: h="**(four2Spec2D)**"

    if (.not. self%created) then
       write(unit=p_nferr, fmt= '(2A)') h, &
             ' Module not created; invoke InitLegTrans prior to this call'
       stop
    end if

    specZ = size(spec,2)
    FourZ = size(four,3)

    if (specZ /= fourZ) then
       write(unit=p_nferr, fmt='(2A,2I10)') h, &
             ' vertical layers of spec and four dissagre: ', specZ, fourZ
       stop
    end if

    do k = 1, specZ
       call self%four2Spec(four(:,:,k),spec(:,k))
    enddo

  end subroutine four2Spec2D



! Eduardo Khamis 17-04-2019
  subroutine splitTrans2D (full, north, south)
    implicit none

    real(kind=p_r8), intent(in   ) :: full (:,:)
    real(kind=p_r8), intent(  out) :: north(:,:)
    real(kind=p_r8), intent(  out) :: south(:,:)

    character(len=18), parameter :: h = "**(splitTrans2D)**"

    integer :: if1, in1, is1
    integer :: if2, in2, is2
    integer :: i, j

    if1=size(full,1); in1=size(north,1); is1=size(south,1)
    if2=size(full,2); in2=size(north,2); is2=size(south,2)

    if (in1 /= is1) then
       write(unit=p_nferr, fmt='(2A,2I10)') h, &
             ' dim 1 of north and south dissagree: ', in1, is1
       stop
    end if
    if (in2 /= is2) then
       write(unit=p_nferr, fmt='(2A,2I10)') h, &
             ' dim 2 of north and south dissagree: ', in2, is2
       stop
    end if

    if (in1 < if1) then
       write(unit=p_nferr, fmt='(2A,2I10)') h, &
             ' first dimension of north too small: ', in1, if1
       stop
    end if
    if (if2 /= 2*in2) then
       write(unit=p_nferr, fmt='(2A,2I10)') h, &
             ' second dimension of full /= 2*second dimension of north: ', if2, 2*in2
       stop
    end if

    north=0.0_p_r8
    south=0.0_p_r8
    do j=1,in2
       do i=1,if1
          north(i,j)=full(i,j)
          south(i,j)=full(i,if2-j+1)
       end do
    end do

  end subroutine splitTrans2D

! Eduardo Khamis 17-04-2019
  subroutine splitTrans3D (full, north, south)
    implicit none

    real(kind=p_r8), intent(in   ) :: full (:,:,:)
    real(kind=p_r8), intent(  out) :: north(:,:,:)
    real(kind=p_r8), intent(  out) :: south(:,:,:)

    integer :: if1, in1, is1
    integer :: if2, in2, is2
    integer :: if3, in3, is3
    integer :: k

    character(len=18), parameter :: h="**(splitTrans3D)**"

    if1=size(full,1); in1=size(north,1); is1=size(south,1)
    if2=size(full,2); in2=size(north,2); is2=size(south,2)
    if3=size(full,3); in3=size(north,3); is3=size(south,3)

    if (in3 /= is3) then
       write(unit=p_nferr, fmt='(2A,2I10)') h, &
             ' dim 3 of north and south dissagree: ', in3, is3
       stop
    end if

    if (if3 /= in2) then
       write(unit=p_nferr, fmt='(2A,2I10)') h, &
             ' second dimension of north and third dimension of full dissagree: ', in2, if3
       stop
    end if

    north=0.0_p_r8
    south=0.0_p_r8
    do k=1,in2
       call splitTrans2D(full(:,:,k), north(:,:,k), south(:,:,k))
    end do

  end subroutine splitTrans3D

!  subroutine DivgVortToUV_ ( Mend, qDivg, qVort, qUvel, qVvel)
  subroutine DivgVortToUV ( self, qDivg, qVort, qUvel, qVvel)
  
      ! Calculates Spectral Representation of Cosine-Weighted
      ! Wind Components from Spectral Representation of
      ! Vorticity and Divergence.
    
      ! qDivg Input:  Divergence (Spectral)
      ! qVort Input:  Vorticity  (Spectral)
      ! qUvel Output: Zonal Pseudo-Wind (Spectral)
      ! qVvel Output: Meridional Pseudo-Wind (Spectral)
    
      class(legendre), intent(in) :: self

      real(kind=p_r8), intent(inout) :: qdivg(2,self%MnWv0) ! (2,mnwv0) = (MnWv2)
      real(kind=p_r8), intent(inout) :: qvort(2,self%MnWv0) ! (2,mnwv0) = (MnWv2)
      real(kind=p_r8), intent(  out) :: quvel(2,self%MnWv1) ! (2,mnwv1) = (MnWv3)
      real(kind=p_r8), intent(  out) :: qvvel(2,self%MnWv1) ! (2,mnwv1) = (MnWv3)
      
      integer :: mm, nn, l, l0, l0p, l0m, l1, l1p, Nmax

      qUvel=0.0_p_r8
      qVvel=0.0_p_r8

      qDivg(2,1:self%Mend1)=0.0_p_r8
      qVort(2,1:self%Mend1)=0.0_p_r8

      CALL self%transp ( -1, qDivg, size(qDivg))
      CALL self%transp ( -1, qVort, size(qVort))
   
    !cdir novector
      DO mm=1,self%Mend1

        Nmax=self%Mend2+1-mm

        qUvel(1,mm) =  self%e1(mm)*qDivg(2,mm)
        qUvel(2,mm) = -self%e1(mm)*qDivg(1,mm)

        qVvel(1,mm) =  self%e1(mm)*qVort(2,mm)
        qVvel(2,mm) = -self%e1(mm)*qVort(1,mm)

        IF (Nmax >= 3) THEN
          l=self%Mend1
          qUvel(1,mm) = qUvel(1,mm)+self%e0(mm+l)*qVort(1,mm+l)
          qUvel(2,mm) = qUvel(2,mm)+self%e0(mm+l)*qVort(2,mm+l)

          qVvel(1,mm) = qVvel(1,mm)-self%e0(mm+l)*qDivg(1,mm+l)
          qVvel(2,mm) = qVvel(2,mm)-self%e0(mm+l)*qDivg(2,mm+l)
        END IF

        IF (Nmax >= 4) THEN
          DO nn=2,Nmax-2
            l0  = self%la0(mm,nn)
            l0p = self%la0(mm,nn+1)
            l0m = self%la0(mm,nn-1)
            l1  = self%la1(mm,nn)
            l1p = self%la1(mm,nn+1)
 
            qUvel(1,l1) = -self%e0(l1)*qVort(1,l0m)+self%e0(l1p)*qVort(1,l0p)+ &
                             self%e1(l0)*qDivg(2,l0)
            qUvel(2,l1) = -self%e0(l1)*qVort(2,l0m)+self%e0(l1p)*qVort(2,l0p)- &
                             self%e1(l0)*qDivg(1,l0)

            qVvel(1,l1) =  self%e0(l1)*qDivg(1,l0m)-self%e0(l1p)*qDivg(1,l0p)+ &
                             self%e1(l0)*qVort(2,l0)
            qVvel(2,l1) =  self%e0(l1)*qDivg(2,l0m)-self%e0(l1p)*qDivg(2,l0p)- &
                              self%e1(l0)*qVort(1,l0)
          END DO
        END IF

        IF (Nmax >= 3) THEN

          nn  = Nmax-1
          l0  = self%la0(mm,nn)
          l0m = self%la0(mm,nn-1)
          l1  = self%la1(mm,nn)

          qUvel(1,l1) = -self%e0(l1)*qVort(1,l0m)+self%e1(l0)*qDivg(2,l0)
          qUvel(2,l1) = -self%e0(l1)*qVort(2,l0m)-self%e1(l0)*qDivg(1,l0)

          qVvel(1,l1) =  self%e0(l1)*qDivg(1,l0m)+self%e1(l0)*qVort(2,l0)
          qVvel(2,l1) =  self%e0(l1)*qDivg(2,l0m)-self%e1(l0)*qVort(1,l0)
          
        END IF
        IF (Nmax >= 2) THEN
          nn=Nmax
          l0m=self%la0(mm,nn-1)
          l1 =self%la1(mm,nn)
          qUvel(1,l1)=-self%e0(l1)*qVort(1,l0m)
          qUvel(2,l1)=-self%e0(l1)*qVort(2,l0m)

          qVvel(1,l1)= self%e0(l1)*qDivg(1,l0m)
          qVvel(2,l1)= self%e0(l1)*qDivg(2,l0m)
        END IF
      END DO
    
      CALL self%transp (+1, qUvel, size(qUvel))
      CALL self%transp (+1, qVvel, size(qVvel))
      CALL self%transp (+1, qDivg, size(qDivg))
      CALL self%transp (+1, qVort, size(qVort))
  
  END SUBROUTINE DivgVortToUV

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

!    real(kind=p_r8) :: dl, dx, dKm = 112.0_p_r8

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


!!   For debuging :
!
!    if (linearGrid) then
!       write(unit=*, fmt='(/,A)') &
!             ' Linear Triangular Truncation : '
!    else
!       write(unit=*, fmt='(/,A)') &
!             ' Quadratic Triangular Truncation : '
!    end if
!
!    dl = 360.0_p_r8 / real(xMax, p_r8)
!    dx = dl * dKm
!    write(unit=*, fmt='(/,3(A,I5,/))') &
!          ' mend : ', mend, ' xMax : ', xMax, ' yMax : ', yMax
!    write(unit=*, fmt='(A,F13.9,A)')   ' dl: ', dl, ' degrees'
!    write(unit=*, fmt='(A,F13.2,A,/)') ' dx: ', dx, ' km'

  end subroutine getxMaxyMax

!  subroutine GetSizes(Mend, Mend1, Mend2, MnWv0, MnWv1, MnWv2, MnWv3)
!     integer, intent(in   ) :: Mend
!     integer, intent(  out) :: Mend1
!     integer, intent(  out) :: Mend2
!     integer, intent(  out) :: MnWv0
!     integer, intent(  out) :: MnWv1
!     integer, intent(  out) :: MnWv2
!     integer, intent(  out) :: MnWv3
!
!     Mend1 = Mend+1
!     Mend2 = Mend+2
!     Mnwv2 = Mend1*Mend2
!     Mnwv0 = Mnwv2/2
!     Mnwv3 = Mnwv2+2*Mend1
!     Mnwv1 = Mnwv3/2
!
!     
!  end subroutine

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

  end subroutine

  subroutine DestroyGaussRep( self )
     class(legendre), intent(inout) :: self

     deallocate ( self%colRad )
     deallocate ( self%rCs2   )
     deallocate ( self%wgt    )
     deallocate ( self%gLats  )

  end subroutine

  subroutine DestroyLegTrans( self )
     class(legendre), intent(inout) :: self

     deallocate( self%legS2F          )
     deallocate( self%legDerS2F       )
     deallocate( self%legExtS2F       )
     deallocate( self%legF2S          )
     deallocate( self%legDerNS        )
     deallocate( self%legDerEW        )
     deallocate( self%lenDiag         )
     deallocate( self%lenDiagExt      )
     deallocate( self%lastPrevDiag    )
     deallocate( self%lastPrevDiagExt )    
  
  end subroutine


end module Mod_LegendreTransform
