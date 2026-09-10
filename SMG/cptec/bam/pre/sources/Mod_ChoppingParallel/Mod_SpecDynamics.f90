!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_SpecDynamics </br></br>
!#
!# **Brief**: Module exports routines: </br>
!# InitDZtoUV: Should be invoked once, before any other module routine; sets up
!# local constants and mappings; </br>
!# DZtoUV: Velocity fields from Divergence and Vorticity; use values computed by
!# InitDZtoUV </br>
!#
!# InitUvtodz: Should be invoked once, before any other module routine; sets up
!# local constants and mappings; </br>
!# Uvtodz: Divergence and Vorticity from Velocity fields; use values computed by
!# InitUvtodz </br>
!#
!# Module require values from modules Sizes, AssocLegFunc and Constants </br></br>
!#
!# **Files in:**
!#
!# &bull; ? </br></br>
!#
!# **Files out:**
!#
!# &bull; ? </br></br>
!# 
!# **Author**: Paulo Kubota </br>
!#
!# **Version**: 2.0.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>07-04-2011 - Paulo Kubota  - version: 1.15.0 </li>
!#  <li>26-04-2019 - Denis Eiras   - version: 2.0.0 - some adaptations for modularizing Chopping </li>
!# </ul>
!# @endchanges
!#
!# @bug
!# <ul type="disc">
!#  <li>None items at this time </li>
!# </ul>
!# @endbug
!#
!# @todo
!# <ul type="disc">
!#  <li>None items at this time </li>
!# </ul>
!# @endtodo
!#
!# @documentation
!#
!# For theoretical information, please visit the following link: </br> <http://urlib.net/8JMKD3MGP3W34R/3SME6J2> </br>
!# **&#9993;**<mailto:atende.cptec@inpe.br> </br></br>
!# @enddocumentation
!#
!# @warning
!# Copyright Under GLP-3.0
!# &copy; https://opensource.org/licenses/GPL-3.0
!# @endwarning
!#
!#---


module Mod_SpecDynamics

  use Mod_Utils, only : &
    Epslon 
    ! intent(in)

  use Mod_Parallelism_Group_Chopping, only : &
    maxnodes, & 
    ! intent(in)
    myid

  use Mod_Parallelism_Fourier, only : &
    mygroup_four 
    ! intent(in)

  use Mod_InputParameters, only : &
    EMRad1, & 
    ! intent(in)
    EMRad 
    ! intent(in)

  use Mod_Sizes, only : &
    mMax, & ! intent(in)
    nMax, & ! intent(in)
    nExtMax, & ! intent(in)
    mnMax, & ! intent(in)
    mnExtMax, & ! intent(in)
    mnExtMap, & ! intent(in)
    mnMap, & ! intent(in)
    nMap, & ! intent(in)
    mMap, & ! intent(in)
    lm2m, & ! intent(in)
    kMax, & ! intent(in)
    kMaxloc, & ! intent(in)
    myfirstlev, & ! intent(in)
    mylastlev, & ! intent(in)
    mymMax, & ! intent(in)
    mymnMax, & ! intent(in)
    mymnExtMax, & ! intent(in)
    mymnExtMap, & ! intent(in)
    mymnMap, & ! intent(in)
    mynMap, & ! intent(in)
    mymextMap, & ! intent(in)
    mynextMap, & ! intent(in)
    mymMap, & ! intent(in)
    haveM1, & ! intent(in)
    havesurf, & ! intent(in)
    nodehasm, & ! intent(in)
    ngroups_four, & ! intent(in)
    nlevperg_four, & ! intent(in)
    rpi, & ! intent(in)
    del, & ! intent(in)
    ci              ! intent(in)

  implicit none
  private
  include 'precision.h'

  public :: InitDztouv
  public :: dztouv
  public :: InitUvtodz
  public :: uvtodz
  public :: Clear_SpecDynamics
  real(kind = p_r8) :: er, eriv

  real(kind = p_r8), allocatable :: alfa_dz(:)        
  ! er*Epslon(m,n)/n  for m<n<=nExtMax; 0 otherwise
  real(kind = p_r8), allocatable :: alfa_dzNp1(:)     
  ! alfa_dz(m,n+1)
  real(kind = p_r8), allocatable :: beta_dz(:)        
  ! m*er/(n*(n+1)) for m/=0 and m<=n<=nMax;
  ! er/(n+1)     for m=0;


  ! Observe, in the relation to be computed, that u and v are defined for
  ! 1<=mn<=mnExtMax, while Div and Vor for 1<=mn<=mnMax. Consequently, a
  ! mapping from 1:mnExtMax to 1:mnMax
  ! has to be computed.

  ! This mapping will have faults, since there is no Div or Vor at (m,nExtMax).

  ! Furthermore, the relation requires mappings from (m,nExt) to (m,n),
  ! (m,n+1) and (m,n-1)

  ! Mapping function mnp1_dz(1:2*mnExtMax) gives index of (m,nExt) on (m,n+1).
  ! It is faulty on (*,nMax:nExtMax). Since it is only used by the last term,
  ! faulty values have to be multipied by 0 (on alfa_dzNp1)

  ! Mapping function mnm1_dz(1:2*mnExtMax) gives index of (m,nExt) on (m,n-1).
  ! It is faulty on (m,m) for all m. Since it is only used by the second term,
  ! faulty values have to be multipied by 0 (on alfa_dz)

  ! Mapping function mnir_dz(1:2*mnExtMax) gives index of (m,nExt) on (m,n-1) and
  ! multiplies by i (trading imaginary by real and correcting sign). It is
  ! faulty on (m,nExtMax). To correct the fault, beta_dz(m,nExtMax) is set to 0.


  integer, allocatable :: mnir_dz(:)
  integer, allocatable :: mnm1_dz(:)
  integer, allocatable :: mnp1_dz(:)


  real(kind = p_r8), allocatable :: alfa_uv(:)        
  ! eriv*m
  ! for m<=n<=nMax
  real(kind = p_r8), allocatable :: beta_uv(:)        
  ! eriv*n*Epslon(m,n+1)
  ! for m<=n<=nMax;
  real(kind = p_r8), allocatable :: gama_uv(:)        
  ! eriv*(n+1)*Epslon(m,n)
  ! for m<=n<=nMax;

  ! Observe, in the relation to be computed, that u and v are defined for
  ! 1<=mn<=mnExtMax, while Div and Vor for 1<=mn<=mnMax. Consequently, a
  ! mapping from 1:mnMax to 1:mnExtMax has to be computed.

  ! In fact, the relation requires 3 mappings:
  !    1) (m,n) to (m,nExt)    implemented by index array mnir(mn);
  !    2) (m,n) to (m,nExt+1)  implemented by index array mnp1(mn);
  !    3) (m,n) to (m,nExt-1)  implemented by index array mnm1(mn);
  ! for m<=n<=nMax

  ! Mapping functions (1) and (2) are easily computed; mapping function (3)
  ! will be faulty for n=m; since it is used only at the last term,
  ! faulty values have to be multipied by 0 (on gama_uv(m,m))

  ! Mapping function mnir(1:2*mnMax) gives index of (m,n) on (m,nExt) and
  ! multiplies by i (trading immaginary by real). Other mappings keep in
  ! place the real and immaginary components.

  integer, allocatable :: mnir_uv(:)
  integer, allocatable :: mnm1_uv(:)
  integer, allocatable :: mnp1_uv(:)

  integer :: kGlob
  real(kind = p_r8) :: alphaGlob
  real(kind = p_r8) :: betaGlob
  integer, allocatable :: ncrit(:)

contains

  subroutine InitDZtoUV()
    !# Initializes Velocity fields from Divergence and Vorticity
    !# ---
    !# @info
    !# **Brief:** Should be invoked once, before any other module routine;
    !# Mapping Functions and Local Constants.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin  
    integer :: m, mglob, n, mn, mn2, indir, indnp1, indnm1
    real(kind = p_r8) :: aux

    ! mapping mnir_dz
    er = EMRad
    eriv = EMRad1
    allocate(mnir_dz(2 * mymnExtMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,indir)
    do m = 1, mymMax
      mglob = lm2m(m)
      do n = mglob, nMax
        mn = mymnExtMap(m, n)
        indir = mymnMap(m, n)
        mnir_dz(2 * mn - 1) = 2 * indir
        mnir_dz(2 * mn) = 2 * indir - 1
      end do
      mn = mymnExtMap(m, nExtMax)
      mnir_dz(2 * mn - 1) = 1     
      ! faulty mapping # 1
      mnir_dz(2 * mn) = 1     
      ! faulty mapping # 1
    end do
    !$OMP END PARALLEL DO

    ! mapping mnm1_dz

    allocate(mnm1_dz  (2 * mymnExtMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,indnm1)
    do m = 1, mymMax
      mglob = lm2m(m)
      mn = mymnExtMap(m, mglob)
      indnm1 = mymnMap(m, mglob)
      mnm1_dz(2 * mn - 1) = 1    
      ! faulty mapping # 2
      mnm1_dz(2 * mn) = 1    
      ! faulty mapping # 2
      do n = mglob + 1, nExtMax
        mn = mymnExtMap(m, n)
        indnm1 = mymnMap(m, n - 1)
        mnm1_dz(2 * mn - 1) = 2 * indnm1 - 1
        mnm1_dz(2 * mn) = 2 * indnm1
      end do
    end do
    !$OMP END PARALLEL DO

    ! mapping mnp1_dz

    allocate(mnp1_dz  (2 * mymnExtMax))
    mnp1_dz = 0
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,indnp1)
    do m = 1, mymMax
      mglob = lm2m(m)
      do n = mglob, nMax - 1
        mn = mymnExtMap(m, n)
        indnp1 = mymnMap(m, n + 1)
        mnp1_dz(2 * mn - 1) = 2 * indnp1 - 1
        mnp1_dz(2 * mn) = 2 * indnp1
      end do
      do n = nMax, nExtMax
        mn = mymnExtMap(m, n)
        mnp1_dz(2 * mn - 1) = 1 
        ! faulty mapping # 3
        mnp1_dz(2 * mn) = 1 
        ! faulty mapping # 3
      end do
    end do
    !$OMP END PARALLEL DO

    ! constant beta_dz

    allocate(beta_dz(2 * mymnExtMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,aux)
    do m = 1, mymMax
      mglob = lm2m(m)
      aux = er / real(mglob, p_r8)
      mn = mymnExtMap(m, mglob)
      beta_dz(2 * mn - 1) = aux
      beta_dz(2 * mn) = -aux
      do n = mglob + 1, nMax
        aux = real(mglob - 1, p_r8) * er / real((n - 1) * n, p_r8)
        mn = mymnExtMap(m, n)
        beta_dz(2 * mn - 1) = aux
        beta_dz(2 * mn) = -aux
      end do
      mn = mymnExtMap(m, nExtMax)
      beta_dz(2 * mn - 1) = 0.0_p_r8           
      ! corrects faulty mapping # 1
      beta_dz(2 * mn) = 0.0_p_r8           
      ! corrects faulty mapping # 1
    end do
    !$OMP END PARALLEL DO

    ! constant alfa_dz

    allocate(alfa_dz    (2 * mymnExtMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,aux)
    do m = 1, mymMax
      mglob = lm2m(m)
      mn = mymnExtMap(m, mglob)
      alfa_dz(2 * mn - 1) = 0.0_p_r8           
      ! corrects faulty mapping # 2
      alfa_dz(2 * mn) = 0.0_p_r8           
      ! corrects faulty mapping # 2
      do n = mglob + 1, nExtMax
        mn = mymnExtMap(m, n)
        aux = er * Epslon(mn) / real(n - 1, p_r8)
        alfa_dz(2 * mn - 1) = aux
        alfa_dz(2 * mn) = aux
      end do
    end do
    !$OMP END PARALLEL DO

    ! constant alfa_dz mapped to n-1

    allocate(alfa_dzNp1 (2 * mymnExtMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,mn2,aux)
    do m = 1, mymMax
      mglob = lm2m(m)
      do n = mglob, nMax - 1
        mn = mymnExtMap(m, n)
        mn2 = mymnExtMap(m, n + 1)
        aux = er * Epslon(mn2) / real(n, p_r8)
        alfa_dzNp1(2 * mn - 1) = aux
        alfa_dzNp1(2 * mn) = aux
      end do
      do n = nMax, nExtMax
        mn = mymnExtMap(m, n)
        alfa_dzNp1(2 * mn - 1) = 0.0_p_r8     
        ! corrects faulty mapping # 3
        alfa_dzNp1(2 * mn) = 0.0_p_r8     
        ! corrects faulty mapping # 3
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine InitDZtoUV


  subroutine DZtoUV(qdivp, qrotp, qup, qvp, mnRIExtFirst, mnRIExtLast)
    !# Velocity fields from Divergence and Vorticity
    !# ---
    !# @info
    !# **Brief:** Use values computed by InitDZtoUV. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
       
    ! Implements the following relations: </br>
    ! m            m          m       m         m         m     m       m     m </br>
    ! U  =CMPLX(Beta * Imag(Div), -Beta *Real(Div))  - alfa * Vor  + alfa * Vor </br>
    !  n            n          n       n         n         n    n-1     n+1   n+1 </br>
    !  m            m          m       m         m         m     m       m     m </br>
    ! V  =CMPLX(Beta * Imag(Vor), -Beta *Real(Vor))  + alfa * Div  - alfa * Div </br>
    !  n            n          n       n         n         n    n-1     n+1   n+1 </br>
    ! 
    ! for 0<=m<=mMax, m<=n<=nExtMax, where </br>
    !    m     m </br>
    ! Div = Vor = 0 for n > nMax or n < m </br>
    !    n     n  </br>

    real(kind = p_r8), intent(in) :: qdivp(2 * mymnMax, kMaxloc)
    real(kind = p_r8), intent(in) :: qrotp(2 * mymnMax, kMaxloc)
    real(kind = p_r8), intent(OUT) :: qup(2 * mymnExtMax, kMaxloc)
    real(kind = p_r8), intent(OUT) :: qvp(2 * mymnExtMax, kMaxloc)
    integer, intent(in) :: mnRIExtFirst
    integer, intent(in) :: mnRIExtLast
    integer :: mn, k

    do k = 1, kMaxloc
      do mn = mnRIExtFirst, mnRIExtLast
        qup(mn, k) = - &
          alfa_dz   (mn) * qrotp(mnm1_dz(mn), k) + &
          alfa_dzNp1(mn) * qrotp(mnp1_dz(mn), k) + &
          beta_dz   (mn) * qdivp(mnir_dz(mn), k)
        qvp(mn, k) = + &
          alfa_dz   (mn) * qdivp(mnm1_dz(mn), k) - &
          alfa_dzNp1(mn) * qdivp(mnp1_dz(mn), k) + &
          beta_dz   (mn) * qrotp(mnir_dz(mn), k)
      end do
    end do
  end subroutine DZtoUV


  subroutine InitUvtodz()
    !# Initializes DIVERGENCE AND VORTICITY FROM VELOCITY FIELDS
    !# ---
    !# @info
    !# **Brief:** Should be invoked once, before any other module routine;
    !# Mapping Functions and Local Constants </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin   
    integer :: m, mglob, n, mn, mnExt

    ! mapping mnir_uv

    er = EMRad
    eriv = EMRad1
    allocate(mnir_uv(2 * mymnMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,mnExt)
    do m = 1, mymMax
      mglob = lm2m(m)
      do n = mglob, nMax
        mn = mymnMap(m, n)
        mnExt = mymnExtMap(m, n)
        mnir_uv(2 * mn - 1) = 2 * mnExt
        mnir_uv(2 * mn) = 2 * mnExt - 1
      end do
    end do
    !$OMP END PARALLEL DO

    ! mapping mnm1_uv

    allocate(mnm1_uv(2 * mymnMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,mnExt)
    do m = 1, mymMax
      mglob = lm2m(m)
      mn = mymnMap(m, mglob)
      mnExt = mymnExtMap(m, mglob)
      mnm1_uv(2 * mn - 1) = 2 * mnExt - 1   
      ! faulty mapping
      mnm1_uv(2 * mn) = 2 * mnExt     
      ! faulty mapping
      do n = mglob + 1, nMax
        mn = mymnMap(m, n)
        mnExt = mymnExtMap(m, n - 1)
        mnm1_uv(2 * mn - 1) = 2 * mnExt - 1
        mnm1_uv(2 * mn) = 2 * mnExt
      end do
    end do
    !$OMP END PARALLEL DO

    ! mapping mnp1_uv

    allocate(mnp1_uv(2 * mymnMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,mnExt)
    do m = 1, mymMax
      mglob = lm2m(m)
      do n = mglob, nMax
        mn = mymnMap(m, n)
        mnExt = mymnExtMap(m, n + 1)
        mnp1_uv(2 * mn - 1) = 2 * mnExt - 1
        mnp1_uv(2 * mn) = 2 * mnExt
      end do
    end do
    !$OMP END PARALLEL DO

    ! constant alfa_uv

    allocate(alfa_uv(2 * mymnMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn)
    do m = 1, mymMax
      mglob = lm2m(m)
      do n = mglob, nMax
        mn = mymnMap(m, n)
        alfa_uv(2 * mn - 1) = -REAL(mglob - 1, p_r8)
        alfa_uv(2 * mn) = real(mglob - 1, p_r8)
      end do
    end do
    !$OMP END PARALLEL DO

    ! constant beta_uv

    allocate(beta_uv(2 * mymnMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,mnExt)
    do m = 1, mymMax
      mglob = lm2m(m)
      do n = mglob, nMax
        mn = mymnMap(m, n)
        mnExt = mymnExtMap(m, n + 1)
        beta_uv(2 * mn - 1) = real(n - 1, p_r8) * Epslon(mnExt)
        beta_uv(2 * mn) = real(n - 1, p_r8) * Epslon(mnExt)
      end do
    end do
    !$OMP END PARALLEL DO

    ! constant gama_uv

    allocate(gama_uv(2 * mymnExtMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,mnExt)
    do m = 1, mymMax
      mglob = lm2m(m)
      mn = mymnMap(m, mglob)
      gama_uv(2 * mn - 1) = 0.0_p_r8     
      ! corrects faulty mapping
      gama_uv(2 * mn) = 0.0_p_r8     
      ! corrects faulty mapping
      do n = mglob + 1, nMax
        mn = mymnMap(m, n)
        mnExt = mymnExtMap(m, n)
        gama_uv(2 * mn - 1) = real(n, p_r8) * Epslon(mnExt)
        gama_uv(2 * mn) = real(n, p_r8) * Epslon(mnExt)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine InitUvtodz


  subroutine Uvtodz(qup, qvp, qdivt, qrott, mnRIFirst, mnRILast)
    !# Divergence and Vorticity from Velocity fields
    !# ---
    !# @info
    !# **Brief:** Use values computed by InitUvtodz. Obtains divergence and
    !# vorticity tendencies.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
        
    ! Implements the following relations: </br>
    !   m              m        m       m        m         m   m       m   m </br>
    ! Div   =CMPLX(-Alfa * Imag(U ), Alfa * Real(U )) + Beta * V  - Gama * V </br>
    !    n              n        n       n        n        n+1 n+1      n  n-1 </br>
    !
    !    m              m        m       m        m         m   m       m   m </br>
    ! Vor   =CMPLX(-Alfa * Imag(V ), Alfa * Real(V )) + Beta * U  - Gama * U </br>
    !    n              n        n       n        n        n+1 n+1      n  n-1 </br>
    !
    ! for 0<=m<=mMax, m<=n<=nMax, where </br>
    !  m   m  </br>
    ! U = V = 0 for n < m </br>
    !  n   n </br>

    real(kind = p_r8), intent(in) :: qup(2 * mymnExtMax, kMaxloc)
    real(kind = p_r8), intent(in) :: qvp(2 * mymnExtMax, kMaxloc)
    real(kind = p_r8), intent(OUT) :: qdivt(2 * mymnMax, kMaxloc)
    real(kind = p_r8), intent(OUT) :: qrott(2 * mymnMax, kMaxloc)
    integer, intent(in) :: mnRIFirst
    integer, intent(in) :: mnRILast
    integer :: mn, k

    do k = 1, kMaxloc
      do mn = mnRIFirst, mnRILast
        qdivt(mn, k) = &
          alfa_uv(mn) * qup(mnir_uv(mn), k) + &
            beta_uv(mn) * qvp(mnp1_uv(mn), k) - &
            gama_uv(mn) * qvp(mnm1_uv(mn), k)
        qdivt(mn, k) = eriv * qdivt(mn, k)
        qrott(mn, k) = &
          alfa_uv(mn) * qvp(mnir_uv(mn), k) - &
            beta_uv(mn) * qup(mnp1_uv(mn), k) + &
            gama_uv(mn) * qup(mnm1_uv(mn), k)
        qrott(mn, k) = eriv * qrott(mn, k)
      end do
    end do
  end subroutine Uvtodz

  subroutine Clear_SpecDynamics()
    !# Cleans Spec Dynamics
    !# ---
    !# @info
    !# **Brief:** Cleans Spec Dynamics.  </br>
    !# **Authors**: </br>
    !# &bull; Paulo Kubota </br>
    !# **Date**: abr/2011 <br>
    !# @endin
    deallocate (alfa_dz)        
    ! er*Epslon(m,n)/n  for m<n<=nExtMax; 0 otherwise
    deallocate (alfa_dzNp1)     
    ! alfa_dz(m,n+1)
    deallocate (beta_dz)        
    ! m*er/(n*(n+1)) for m/=0 and m<=n<=nMax;
    deallocate (mnir_dz)
    deallocate (mnm1_dz)
    deallocate (mnp1_dz)
    deallocate (alfa_uv)        
    ! eriv*m
    deallocate (beta_uv)        
    ! eriv*n*Epslon(m,n+1)
    deallocate (gama_uv)        
    ! eriv*(n+1)*Epslon(m,n)
    deallocate (mnir_uv)
    deallocate (mnm1_uv)
    deallocate (mnp1_uv)
    ! DEALLOCATE ( ncrit )
  end subroutine Clear_SpecDynamics

end module Mod_SpecDynamics
