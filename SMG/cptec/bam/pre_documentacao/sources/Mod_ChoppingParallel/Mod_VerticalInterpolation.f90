!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_VerticalInterpolation </br></br>
!#
!# **Brief**: Module containig subroutines for Vertical Interpolation. </br></br>
!#
!# **Files in:**
!#
!# &bull; ? </br></br>
!#
!# **Files out:**
!#
!# &bull; ? </br></br>
!# 
!# **Author**: From NCEP </br>
!#
!# **Version**: 2.1.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>01-08-2007 - Simone Tomita  - version: 1.1.1 </li>
!#  <li>26-04-2019 - Denis Eiras    - version: 2.0.0 - some adaptations for modularizing Chopping </li>
!#  <li>09-10-2019 - Eduardo Khamis - version: 2.1.0 - changing for operational Chopping </li>
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

module Mod_VerticalInterpolation

  use Mod_InputParameters, only : Gama, Grav, Rd, Lc, Rv, Cp

  implicit none

  private
  include 'precision.h'
  include 'messages.h'

  public :: VertSigmaInter

contains


  subroutine VertSigmaInter (ibdim, Im, Km1, Km2, Nt, &
    p1, u1, v1, t1, q1, &
    p2, u2, v2, t2, q2)
    !# Vertical Sigma Interpolation  
    !# ---
    !# @info
    !# **Brief:** Vertically Interpolate Upper-Air Fields:
    !# <ul type="disc">
    !#  <li>Wind, Temperature, Humidity and other Tracers are Interpolated. </li>
    !#  <li>The Interpolation is Cubic Lagrangian in Log Pressure with a Monotonic
    !# Constraint in the Center of the Domain. </li>
    !#  <li>In the Outer Intervals it is Linear in Log Pressure. </li>
    !#  <li>Outside the Domain, Fields are Generally Held Constant, Except for
    !# Temperature and Humidity Below the Input Domain, Where the Temperature
    !# Lapse Rate is Held Fixed at -6.5 K/km and the Relative Humidity is Held
    !# Constant. </li>
    !# </ul>
    !# <ul type="disc">
    !#  Input Argument List: <ul type="disc">
    !#  <li>Im  - First Dimension </li>
    !#  <li>Km1 - Number of Input Levels </li>
    !#  <li>Km2 - Number of Output Levels </li>
    !#  <li>Nt  - Number of Tracers </li>
    !#  <li>p1  - Input Pressures (Ordered from Bottom to Top of Atmosphere) </li>
    !#  <li>u1  - Input Zonal Wind </li>
    !#  <li>v1  - Input Meridional Wind </li>
    !#  <li>t1  - Input Temperature (K) </li>
    !#  <li>q1  - Input Tracers (Specific Humidity First) </li>
    !#  <li>p2  - Output Pressures </li> 
    !#  </ul></ul>
    !# <ul type="disc">
    !#  Output Argument List: <ul type="disc">
    !#  <li>u2 - Output Zonal Wind </li>
    !#  <li>v2 - Output Meridional Wind </li>
    !#  <li>t2 - Output Temperature (K) </li>
    !#  <li>q2 - Output Tracers (Specific Humidity First) </li>
    !#  </ul></ul>
    !# **Authors**: </br>
    !# &bull; From NCEP </br>
    !# **Date**: Early 2003 <br>
    !# @endin
    implicit none

    integer, intent (in) :: ibdim, Im, Km1, Km2, Nt

    real (kind = p_r8), dimension (Ibdim, Km1), intent (in) :: p1, u1, v1, t1
    real (kind = p_r8), dimension (Ibdim, Km1, Nt), intent (in) :: q1
    real (kind = p_r8), dimension (Ibdim, Km2), intent (in) :: p2
    real (kind = p_r8), dimension (Ibdim, Km2), intent (OUT) :: u2, v2, t2
    real (kind = p_r8), dimension (Ibdim, Km2, Nt), intent (OUT) :: q2

    integer :: k, i, n, k1

    real (kind = p_r8) :: dltdz, dlpvdrt, dz, RdByCp
    real (kind = p_r8), dimension (Ibdim, Km1) :: z1
    real (kind = p_r8), dimension (Ibdim, Km2) :: z2
    real (kind = p_r8), dimension (Ibdim, Km1, 3 + Nt) :: c1
    real (kind = p_r8), dimension (Ibdim, Km2, 3 + Nt) :: c2
    integer           , dimension (Ibdim,Km2) :: k1s

    RdByCp = Rd / Cp

    dltdz = Gama * Rd / Grav
    dlpvdrt = -Lc / Rv

    ! Compute Log Pressure Interpolating Coordinate and
    ! Copy Input Wind, Temperature, Humidity and other Tracers

    do k = 1, Km1
      do i = 1, Im
        z1(i, k) = -LOG(p1(i, k))
        c1(i, k, 1) = u1(i, k)
        c1(i, k, 2) = v1(i, k)
        c1(i, k, 3) = t1(i, k)
        c1(i, k, 4) = q1(i, k, 1)
      enddo
    enddo
    do n = 2, Nt
      do k = 1, Km1
        do i = 1, Im
          c1(i, k, 3 + n) = q1(i, k, n)
        enddo
      enddo
    enddo
    do k = 1, Km2
      do i = 1, Im
        z2(i, k) = -log(p2(i, k))
      enddo
    enddo

    ! Perform Lagrangian One-Dimensional Interpolation that is
    ! 4th-Order in Interior,
    ! 2nd-Order in Outside Intervals and
    ! 1st-Order for Extrapolation

    call terp3 (ibdim, Im, Km1, Km2, 3 + Nt, z1, c1, z2, c2,k1s)

    ! Copy Output Wind, Temperature, Specific Humidity and other Tracers
    ! Except Below the Input Domain, Let Temperature Increase with a
    ! Fixed Lapse Rate and Let the Relative Humidity Remain Constant

    do k = 1, Km2
      do i = 1, Im
        k1=k1s(i,k)
        t2(i,k)  =c2(i,k,3)
        q2(i,k,1)=c2(i,k,4)
        u2(i,k)  =c2(i,k,1)
        v2(i,k)  =c2(i,k,2)
        dz=z2(i,k)-z1(i,1)
        if (dz >= 0.0_p_r8) then
          t2(i,k)   = c2(i,k,3)
          q2(i,k,1) = c2(i,k,4)
          u2(i,k)  =c2(i,k,1)
          v2(i,k)  =c2(i,k,2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if(Km2 == Km1)then
           !IF(k > 1.and.k < Km2  .and. p2(i,k) <= p1(i,1) .and. p2(i,k) >= p1(i,Km1)  )THEN
           if(p2(i,k) <= p1(i,1) .and. p2(i,k) >= p1(i,km1)  )then
              if(k > 1.and.k < km2)then
                 t2(i,k)  =0.25_p_r8*c2(i,k-1,3)+0.5_p_r8*c2(i,k,3)+0.25_p_r8*c2(i,k+1,3)
                 q2(i,k,1)=0.25_p_r8*c2(i,k-1,4)+0.5_p_r8*c2(i,k,4)+0.25_p_r8*c2(i,k+1,4)
                 u2(i,k)  =0.25_p_r8*c2(i,k-1,1)+0.5_p_r8*c2(i,k,1)+0.25_p_r8*c2(i,k+1,1)!c2(i,k,1)
                 v2(i,k)  =0.25_p_r8*c2(i,k-1,2)+0.5_p_r8*c2(i,k,2)+0.25_p_r8*c2(i,k+1,2)!c2(i,k,2)
              else if(k == 1 )then
                 t2(i,k)  =t1(i,1)
                 q2(i,k,1)=q1(i,1,1)
                 u2(i,k)  =u1(i,1)
                 v2(i,k)  =v1(i,1)
              else if(k == km2 )then
                 t2(i,k)  =0.25_p_r8*t2(i,k-2  )+0.5_p_r8*t2(i,k-1  )+0.25_p_r8*c2(i,k,3)	    !t1(i,k) 
                 q2(i,k,1)=0.25_p_r8*q2(i,k-2,1)+0.5_p_r8*q2(i,k-1,1)+0.25_p_r8*c2(i,k,4)	    !q1(i,k,1)
                 u2(i,k)  =0.25_p_r8*u2(i,k-2  )+0.5_p_r8*u2(i,k-1  )+0.25_p_r8*c2(i,k,1)!c2(i,k,1) !u1(i,k)
                 v2(i,k)  =0.25_p_r8*v2(i,k-2  )+0.5_p_r8*v2(i,k-1  )+0.25_p_r8*c2(i,k,2)!c2(i,k,2) !v1(i,k)
              end if
          else !if(p2(i,k) > p1(i,1) .and. p2(i,k) <= p1(i,km1)  )then
              if(k == 1 )then
                 t2(i,k)  =t1(i,1)
                 q2(i,k,1)=q1(i,1,1)
                 u2(i,k)  =u1(i,1)
                 v2(i,k)  =v1(i,1)
              else if( p2(i,k) < p1(i,km1) )then
                 t2(i,k)  =0.25_p_r8*t2(i,k-2  )+0.5_p_r8*t2(i,k-1  )+0.25_p_r8*c2(i,k,3)	    !t1(i,k) 
                 q2(i,k,1)=0.25_p_r8*q2(i,k-2,1)+0.5_p_r8*q2(i,k-1,1)+0.25_p_r8*c2(i,k,4)	    !q1(i,k,1)
                 u2(i,k)  =0.25_p_r8*u2(i,k-2  )+0.5_p_r8*u2(i,k-1  )+0.25_p_r8*c2(i,k,1)!c2(i,k,1) !u1(i,k)
                 v2(i,k)  =0.25_p_r8*v2(i,k-2  )+0.5_p_r8*v2(i,k-1  )+0.25_p_r8*c2(i,k,2)!c2(i,k,2) !v1(i,k)
              else
                 t2(i,k)  =c2(i,k,3)
                 q2(i,k,1)=c2(i,k,4)
                 u2(i,k)  =c2(i,k,1)
                 v2(i,k)  =c2(i,k,2)
              end if
	  end if

        else
           if(k > 1 .and. k < km2)then


              t2(i,k)  =0.25_p_r8*c2(i,k-1,3)+0.5_p_r8*c2(i,k,3)+0.25_p_r8*c2(i,k+1,3)
              q2(i,k,1)=0.25_p_r8*c2(i,k-1,4)+0.5_p_r8*c2(i,k,4)+0.25_p_r8*c2(i,k+1,4)
              u2(i,k)  =0.25_p_r8*c2(i,k-1,1)+0.5_p_r8*c2(i,k,1)+0.25_p_r8*c2(i,k+1,1)!c2(i,k,1)
              v2(i,k)  =0.25_p_r8*c2(i,k-1,2)+0.5_p_r8*c2(i,k,2)+0.25_p_r8*c2(i,k+1,2)!c2(i,k,2)

           else if(k == 1)then
              t2(i,k)  =t1(i,1)
              q2(i,k,1)=q1(i,1,1)
              u2(i,k)  =u1(i,1)
              v2(i,k)  =v1(i,1)
           else if(k==km2 .and. km2 < km1)then
              t2(i,k)  =c2(i,k,3)
              q2(i,k,1)=c2(i,k,4)
              u2(i,k)  =c2(i,k,1)
              v2(i,k)  =c2(i,k,2)
           else
              t2(i,k)  =t1(i,1)*exp(dltdz*dz)
              q2(i,k,1)=q1(i,1,1)*exp(dlpvdrt*(1.0_p_r8/t2(i,k)-1.0_p_r8/t1(i,1))-dz)
              u2(i,k)  =u1(i,1)
              v2(i,k)  =v1(i,1)
           end if
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(km2 == km1)then
           !if(k > 1.and.k < km2  .and. p2(i,k) <= p1(i,1) .and. p2(i,k) >= p1(i,km1)  )then
           if(p2(i,k) <= p1(i,1) .and. p2(i,k) >= p1(i,km1)  )then
              if(k > 1.and.k < km2)then
                 t2(i,k)  =0.25_p_r8*c2(i,k-1,3)+0.5_p_r8*c2(i,k,3)+0.25_p_r8*c2(i,k+1,3)
                 q2(i,k,1)=0.25_p_r8*c2(i,k-1,4)+0.5_p_r8*c2(i,k,4)+0.25_p_r8*c2(i,k+1,4)
                 u2(i,k)  =0.25_p_r8*c2(i,k-1,1)+0.5_p_r8*c2(i,k,1)+0.25_p_r8*c2(i,k+1,1)!c2(i,k,1)
                 v2(i,k)  =0.25_p_r8*c2(i,k-1,2)+0.5_p_r8*c2(i,k,2)+0.25_p_r8*c2(i,k+1,2)!c2(i,k,2)
              else if(k == 1 )then
                 t2(i,k)  =t1(i,1)
                 q2(i,k,1)=q1(i,1,1)
                 u2(i,k)  =u1(i,1)
                 v2(i,k)  =v1(i,1)
              else if(k == km2 )then
                 t2(i,k)  =0.25_p_r8*t2(i,k-2  )+0.5_p_r8*t2(i,k-1  )+0.25_p_r8*c2(i,k,3)	    !t1(i,k) 
                 q2(i,k,1)=0.25_p_r8*q2(i,k-2,1)+0.5_p_r8*q2(i,k-1,1)+0.25_p_r8*c2(i,k,4)	    !q1(i,k,1)
                 u2(i,k)  =0.25_p_r8*u2(i,k-2  )+0.5_p_r8*u2(i,k-1  )+0.25_p_r8*c2(i,k,1)!c2(i,k,1) !u1(i,k)
                 v2(i,k)  =0.25_p_r8*v2(i,k-2  )+0.5_p_r8*v2(i,k-1  )+0.25_p_r8*c2(i,k,2)!c2(i,k,2) !v1(i,k)
              end if
          else !if(p2(i,k) > p1(i,1) .and. p2(i,k) <= p1(i,km1)  )then
              if(k == 1 )then
                 t2(i,k)  =t1(i,1)
                 q2(i,k,1)=q1(i,1,1)
                 u2(i,k)  =u1(i,1)
                 v2(i,k)  =v1(i,1)
              else if( p2(i,k) < p1(i,km1) )then
                 t2(i,k)  =0.25_p_r8*t2(i,k-2  )+0.5_p_r8*t2(i,k-1  )+0.25_p_r8*c2(i,k,3)	    !t1(i,k) 
                 q2(i,k,1)=0.25_p_r8*q2(i,k-2,1)+0.5_p_r8*q2(i,k-1,1)+0.25_p_r8*c2(i,k,4)	    !q1(i,k,1)
                 u2(i,k)  =0.25_p_r8*u2(i,k-2  )+0.5_p_r8*u2(i,k-1  )+0.25_p_r8*c2(i,k,1)!c2(i,k,1) !u1(i,k)
                 v2(i,k)  =0.25_p_r8*v2(i,k-2  )+0.5_p_r8*v2(i,k-1  )+0.25_p_r8*c2(i,k,2)!c2(i,k,2) !v1(i,k)
              else
                 t2(i,k)  =c2(i,k,3)
                 q2(i,k,1)=c2(i,k,4)
                 u2(i,k)  =c2(i,k,1)
                 v2(i,k)  =c2(i,k,2)
              end if
	  end if

        else
           if(k > 1 .and. k < km2)then


              t2(i,k)  =0.25_p_r8*c2(i,k-1,3)+0.5_p_r8*c2(i,k,3)+0.25_p_r8*c2(i,k+1,3)
              q2(i,k,1)=0.25_p_r8*c2(i,k-1,4)+0.5_p_r8*c2(i,k,4)+0.25_p_r8*c2(i,k+1,4)
              u2(i,k)  =0.25_p_r8*c2(i,k-1,1)+0.5_p_r8*c2(i,k,1)+0.25_p_r8*c2(i,k+1,1)!c2(i,k,1)
              v2(i,k)  =0.25_p_r8*c2(i,k-1,2)+0.5_p_r8*c2(i,k,2)+0.25_p_r8*c2(i,k+1,2)!c2(i,k,2)

           else if(k == 1)then
              t2(i,k)  =t1(i,1)
              q2(i,k,1)=q1(i,1,1)
              u2(i,k)  =u1(i,1)
              v2(i,k)  =v1(i,1)
           else if(k==km2 .and. km2 < km1)then
              t2(i,k)  =0.25_p_r8*t2(i,k-2  )+0.5_p_r8*t2(i,k-1  )+0.25_p_r8*c2(i,k,3)	    !t1(i,k) 
              q2(i,k,1)=0.25_p_r8*q2(i,k-2,1)+0.5_p_r8*q2(i,k-1,1)+0.25_p_r8*c2(i,k,4)	    !q1(i,k,1)
              u2(i,k)  =0.25_p_r8*u2(i,k-2  )+0.5_p_r8*u2(i,k-1  )+0.25_p_r8*c2(i,k,1)!c2(i,k,1) !u1(i,k)
              v2(i,k)  =0.25_p_r8*v2(i,k-2  )+0.5_p_r8*v2(i,k-1  )+0.25_p_r8*c2(i,k,2)!c2(i,k,2) !v1(i,k)
           else
              t2(i,k)  =t1(i,1)*exp(dltdz*dz)
              q2(i,k,1)=q1(i,1,1)*exp(dlpvdrt*(1.0_p_r8/t2(i,k)-1.0_p_r8/t1(i,1))-dz)
              t2(i,k)  =0.25_p_r8*t2(i,k-2  )+0.5_p_r8*t2(i,k-1  )+0.25_p_r8*c2(i,k,3)	    !t1(i,k) 
              q2(i,k,1)=0.25_p_r8*q2(i,k-2,1)+0.5_p_r8*q2(i,k-1,1)+0.25_p_r8*c2(i,k,4)	    !q1(i,k,1)
              u2(i,k)  =0.25_p_r8*u2(i,k-2  )+0.5_p_r8*u2(i,k-1  )+0.25_p_r8*c2(i,k,1)!c2(i,k,1) !u1(i,k)
              v2(i,k)  =0.25_p_r8*v2(i,k-2  )+0.5_p_r8*v2(i,k-1  )+0.25_p_r8*c2(i,k,2)!c2(i,k,2) !v1(i,k)
           end if
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !  dltdz=Gama*Rd/Grav
        !  Gama=-6.5E-3_p_r8
!          if (Km2 == Km1) then
!             if (k > 1.and.k < Km2) then
!               ! t2(i,k)  =t1(i,1  )*EXP(0.5_p_r8*RdByCp*ABS(log(p1(i,k))-log(p2(i,k))))  
!               ! q2(i,k,1)=q1(i,1,1)*EXP(0.5_p_r8*RdByCp*ABS(log(p1(i,k))-log(p2(i,k))))  
!                t2(i,k)   = 0.25_p_r8*c2(i,k-1,3) + 0.5_p_r8*c2(i,k,3) + 0.25_p_r8*c2(i,k+1,3)
!                q2(i,k,1) = 0.25_p_r8*c2(i,k-1,4) + 0.5_p_r8*c2(i,k,4) + 0.25_p_r8*c2(i,k+1,4)
!             else if (k == 1) then 
!                t2(i,k)   = t1(i,1)!*EXP(dltdz*dz)
!                q2(i,k,1) = q1(i,1,1)!*EXP(dlpvdrt*(1.0_p_r8/t2(i,k)-1.0_p_r8/t1(i,1))-dz)
!             else if (k == Km2) then
!                t2(i,k)   = t1(i,k)  !*EXP(dltdz*dz)
!                q2(i,k,1) = q1(i,k,1)!*EXP(dlpvdrt*(1.0_p_r8/t2(i,k)-1.0_p_r8/t1(i,1))-dz)
!             endif
!          else
!             if (k > 1.and.k < Km2) then
!               ! t2(i,k)  =t1(i,1  )*EXP(0.5_p_r8*RdByCp*ABS(log(p1(i,k))-log(p2(i,k))))  
!               ! q2(i,k,1)=q1(i,1,1)*EXP(0.5_p_r8*RdByCp*ABS(log(p1(i,k))-log(p2(i,k))))  
!                t2(i,k)   = 0.25_p_r8*c2(i,k-1,3) + 0.5_p_r8*c2(i,k,3) + 0.25_p_r8*c2(i,k+1,3)
!                q2(i,k,1) = 0.25_p_r8*c2(i,k-1,4) + 0.5_p_r8*c2(i,k,4) + 0.25_p_r8*c2(i,k+1,4)
!             else if (k == 1) then
!                t2(i,k)   = t1(i,1)!*EXP(dltdz*dz)
!                q2(i,k,1) = q1(i,1,1)!*EXP(dlpvdrt*(1.0_p_r8/t2(i,k)-1.0_p_r8/t1(i,1))-dz)
!             else
!                t2(i,k)   =   t1(i,1)*exp(dltdz*dz)
!                q2(i,k,1) = q1(i,1,1)*exp(dlpvdrt*(1.0_p_r8/t2(i,k)-1.0_p_r8/t1(i,1)) - dz)
!             endif
!          endif


        endif
      enddo
    enddo
  !  DO n=2,nT
  !    DO k=1,Km2
  !      DO i=1,Im
  !        q2(i,k,n)=c2(i,k,3+n)
  !      ENDDO
  !    ENDDO
  !  ENDDO
  
    ! Copy Output Tracers
   
    do n = 2, nT
      do k = 1, Km2
        do i = 1, Im
          dz = z2(i,k) - z1(i,1)
          if (dz >= 0.0_p_r8) then
             q2(i,k,n) = c2(i,k,3+n)
          else
             if (Km2 == Km1) then
                if (k > 1 .and. k < Km2) then
                   q2(i,k,n) = 0.25_p_r8*c2(i,k-1,3+n) + 0.5_p_r8*c2(i,k,3+n) + 0.25_p_r8*c2(i,k+1,3+n)
                else if (k == 1) then
                   q2(i,k,n) = c1(i,1,3+n)
                else if (k == Km2) then
                   q2(i,k,n) = c1(i,k,3+n)
                endif
             else
                if (k > 1 .and. k < Km2) then
                   q2(i,k,n) = 0.25_p_r8*c2(i,k-1,3+n) + 0.5_p_r8*c2(i,k,3+n) + 0.25_p_r8*c2(i,k+1,3+n)
                else if (k == 1) then
                   q2(i,k,n) = c1(i,1,3+n)
                else
                   q2(i,k,n) = (c2(i,k,3+n) + c1(i,Km1,3+n)) / 2.0_p_r8
                endif
             endif
          endif
        enddo
      enddo
    enddo

  end subroutine VertSigmaInter


  subroutine terp3 (ibdim, Im, Km1, Km2, Nm, z1, q1, z2, q2, k1s)
    !# Cubically Interpolates Field(s)  
    !# ---
    !# @info
    !# **Brief:** Cubically interpolates field(s) in one dimension along the
    !# column(s). The interpolation is cubic Lagrangian with a monotonic
    !# constraint in the center of the domain. In the outer intervals it is
    !# linear. Outside the domain, fields are held constant. </br>
    !# <ul type="disc">
    !#   Input Argument List:
    !#   <table>
    !#   <tr> <td style="width:75px"><ul type="disc"><li>Im</li></ul></td>  <td>- Number of Columns</td></tr>
    !#   <tr> <td><ul type="disc"><li>Km1</li></ul></td> <td>- Number of Input Points in Each Column</td></tr>
    !#   <tr> <td><ul type="disc"><li>Km2</li></ul></td> <td>- Number of Output Points in Each Column</td></tr>
    !#   <tr> <td><ul type="disc"><li>Nm</li></ul></td>  <td>- Number of Fields per Column</td></tr>
    !#   <tr> <td><ul type="disc"><li>z1</li></ul></td>  <td>-  Input Coordinate Values in which to Interpolate (z1 Must Be
    !#           Strictly Monotonic in Either Direction)</td></tr>
    !#   <tr> <td><ul type="disc"><li>q1</li></ul></td>  <td>- Input Fields to Interpolate</td></tr>
    !#   <tr> <td><ul type="disc"><li>z2</li></ul></td>  <td>- Output Coordinate Values to which to Interpolate (z2 Need Not
    !#           Be Monotonic)</td></tr>
    !#   </table>
    !# </ul>
    !# **Authors**: </br>
    !# &bull; From NCEP </br>
    !# **Date**: Early 2003 <br>
    !# @endin  
    implicit none

    integer, intent (in) :: ibdim, Im, Km1, Km2, Nm

    real (kind = p_r8), dimension (ibdim, Km1), intent (in) :: z1
    real (kind = p_r8), dimension (ibdim, Km2), intent (in) :: z2
    real (kind = p_r8), dimension (ibdim, Km1, Nm), intent (in) :: q1
    real (kind = p_r8), dimension (ibdim, Km2, Nm), intent (out) :: q2
    integer,            dimension (ibdim, km2),intent (out) :: k1s

    integer :: k1, k2, i, n

    real (kind = p_r8) :: z1a, z1b, z1c, z1d, z2s, q1a, q1b, q1c, q1d, q2s

    !integer, dimension (Ibdim, Km2) :: k1s

    real (kind = p_r8), dimension (Ibdim) :: ffa, ffb, ffc, ffd

    ! Find the Surrounding Input Interval for Each Output Point
    k1s=0
    q2 =0.0_p_r8; z1a=0.0_p_r8; z1b=0.0_p_r8; z1c=0.0_p_r8; z1d=0.0_p_r8
    z2s=0.0_p_r8; q1a=0.0_p_r8; q1b=0.0_p_r8; q1c=0.0_p_r8; q1d=0.0_p_r8
    q2s=0.0_p_r8; ffa=0.0_p_r8; ffb=0.0_p_r8; ffc=0.0_p_r8; ffd=0.0_p_r8

    call rsearch (ibdim, Im, Km1, Km2, z1, z2, k1s)

    ! Generally Interpolate Cubically with Monotonic Constraint
    ! From Two Nearest Input Points on Either Side of the Output Point,
    ! But Within the Two Edge Intervals Interpolate Linearly.
    ! Keep the Output Fields Constant Outside the Input Domain.

    do k2 = 1, Km2
      do i = 1, Im
        k1 = k1s(i, k2)
        if (k1 == 1 .or. k1 == Km1 - 1) then
          z2s = z2(i, k2)
          z1a = z1(i, k1)
          z1b = z1(i, k1 + 1)
          ffa(i) = (z2s - z1b) / (z1a - z1b)
          ffb(i) = (z2s - z1a) / (z1b - z1a)
        else if (k2 == km2 .and. km2 < Km1-1 .and. k1 < Km1-1) then
          z2s=z2(i,k2)
          z1a=z1(i,k1-1)
          z1b=z1(i,k1)
          z1c=z1(i,k1+1)
          z1d=z1(i,k1+2)
          ffa(i)=(z2s-z1b)/(z1a-z1b)*(z2s-z1c)/(z1a-z1c)*(z2s-z1d)/(z1a-z1d)
          ffb(i)=(z2s-z1a)/(z1b-z1a)*(z2s-z1c)/(z1b-z1c)*(z2s-z1d)/(z1b-z1d)
          ffc(i)=(z2s-z1a)/(z1c-z1a)*(z2s-z1b)/(z1c-z1b)*(z2s-z1d)/(z1c-z1d)
          ffd(i)=(z2s-z1a)/(z1d-z1a)*(z2s-z1b)/(z1d-z1b)*(z2s-z1c)/(z1d-z1c)
        elseif (k1 > 1 .and. k1 < Km1 - 1) then
          z2s = z2(i, k2)
          z1a = z1(i, k1 - 1)
          z1b = z1(i, k1)
          z1c = z1(i, k1 + 1)
          z1d = z1(i, k1 + 2)
          ffa(i) = (z2s - z1b) / (z1a - z1b) * (z2s - z1c) / (z1a - z1c) * (z2s - z1d) / (z1a - z1d)
          ffb(i) = (z2s - z1a) / (z1b - z1a) * (z2s - z1c) / (z1b - z1c) * (z2s - z1d) / (z1b - z1d)
          ffc(i) = (z2s - z1a) / (z1c - z1a) * (z2s - z1b) / (z1c - z1b) * (z2s - z1d) / (z1c - z1d)
          ffd(i) = (z2s - z1a) / (z1d - z1a) * (z2s - z1b) / (z1d - z1b) * (z2s - z1c) / (z1d - z1c)
        endif
      enddo

      ! Interpolate

      do n = 1, Nm
        do i = 1, Im
          k1 = k1s(i, k2)
          if (k1 == 0) then
            q2s = q1(i, 1, n)
          elseif (k1 == Km1) then
            q2s = q1(i, Km1, n)
          elseif (k1 == 1 .or. k1 == Km1 - 1) then
            q1a = q1(i, k1, n)
            q1b = q1(i, k1 + 1, n)
            q2s = ffa(i) * q1a + ffb(i) * q1b
          else if (k2 == km2 .and. km2 < Km1-1 .and. k1 < Km1-1) then 
            q1a=q1(i,k1-1,n)
            q1b=q1(i,k1,n)
            q1c=q1(i,k1+1,n)
            q1d=q1(i,k1+2,n)
            q2s=0.50_p_r8*q1(i,k1-2,n) + 0.25_p_r8*q1(i,k1-1,n) + 0.25_p_r8*q1(i,k1,n)
          else
            q1a = q1(i, k1 - 1, n)
            q1b = q1(i, k1, n)
            q1c = q1(i, k1 + 1, n)
            q1d = q1(i, k1 + 2, n)
            q2s = min(max(ffa(i) * q1a + ffb(i) * q1b + ffc(i) * q1c + ffd(i) * q1d, &
              min(q1b, q1c)), max(q1b, q1c))
          endif
          q2(i, k2, n) = q2s
        enddo
      enddo
    enddo

  end subroutine terp3


  subroutine rsearch (ibdim, Im, Km1, Km2, z1, z2, l2)
    !# Search for a Surrounding Real Interval 
    !# ---
    !# @info
    !# **Brief:** Searches monotonic sequences of real numbers for intervals that
    !# surround a given search set of real numbers. The sequences must be
    !# monotonically ascending. The input sequences and sets and the output
    !# locations. May be arbitrarily dimensioned. </br>
    !# <ul type="disc">
    !#  Input Argument List:
    !#  <table>
    !#  <tr> <td style="width:75px"><ul type="disc"><li>Im</li></ul></td>  <td>- number of sequences to search</td></tr>
    !#  <tr> <td><ul type="disc"><li>Km1</li></ul></td> <td>- number of points in each sequence</td></tr>
    !#  <tr> <td><ul type="disc"><li>Km2</li></ul></td> <td>- number of points to search for in each respective sequence</td></tr>
    !#  <tr> <td><ul type="disc"><li>z1</li></ul></td>  <td>- sequence values to search (z1 must be monotonically ascending)</td></tr>
    !#  <tr> <td><ul type="disc"><li>z1</li></ul></td>  <td>- set of values to search for (z2 need not be monotonic</td></tr>
    !#  </table>
    !#  </ul>
    !# <ul type="disc">
    !#  Output Argument List: <ul type="disc">
    !#  <li>l2 - Interval locations having values from 0 to km1 (z2 will be between z1(l2) and z1(l2+1)) </li>
    !#  </ul></ul>
    !# <ul type="disc">
    !# Remarks: <ul type="disc">
    !#  <li>Returned values of 0 or km1 indicate that the given search value is
    !#  outside the range of the sequence. </li>
    !#  <li>If a search value is identical to one of the sequence values then the
    !#      Location returned points to the identical value. If the sequence is not
    !#      strictly monotonic and a search value is identical to more than one of the
    !#      sequence values, then the location returned may point to any of the 
    !#      Identical values. </li>
    !#  <li>To be exact, for each i from 1 to im and for each k from 1 to km2,
    !#      z=z2(i,k) is the search value and l=l2(i,k) is the location returned.
    !#      If l=0, then z is less than the start point z1(i,1). If l=km1, then z is
    !#      greater than or equal to the end point z1(i,km1). Otherwise z is between
    !#      the values z1(i,l) and z1(i,l+1) and may equal the former. </li>        
    !#  </ul></ul>
    !# **Authors**: </br>
    !# &bull; From NCEP </br>
    !# **Date**: Early 2003 <br>
    !# @endin  
    implicit none

    integer, intent(in) :: ibdim, Im, Km1, Km2

    real (kind = p_r8), dimension (Ibdim, Km1), intent(in) :: z1
    real (kind = p_r8), dimension (Ibdim, Km2), intent(in) :: z2

    integer, dimension (Ibdim, Km2), intent(out) :: l2
    integer :: i, k2
    integer, dimension (Km2) :: indx, rc
    l2=0;indx=0;rc=0

    ! Find the Surrounding Input Interval for Each Output Point

    do i = 1, Im

      if (z1(i, 1) <= z1(i, Km1)) then

        ! Input Coordinate is Monotonically Ascending

        call bsrch (Km1, Km2, z1(i, :), z2(i, :), indx, rc)

        do k2 = 1, Km2
          l2(i, k2) = indx(k2) - rc(k2)
        enddo

      else

        ! Input Coordinate is Monotonically Descending

        write (unit = p_nfprt, FMT = '(/,A)') ' Warnning: '
        write (unit = p_nfprt, FMT = '(A)')   ' Input Coordinate is Monotonically Descending'
        write (unit = p_nfprt, FMT = '(A)')   ' The Implemented Binary Search Does not Allowed That'
        write (unit = p_nfprt, FMT = '(A,/)') ' Stopping Computation at subroutine rsearch '
        stop

      endif

    enddo

  end subroutine rsearch


  subroutine bsrch (n, m, x, y, indx, rc)
    !# Localiza os 3 pontos mais próximos a serem interpolados
    !# ---
    !# @info
    !# **Brief:** Localiza os 3 pontos mais próximos a serem interpolados.  </br>
    !# **Authors**: </br>
    !# &bull; From NCEP </br>
    !# **Date**: Early 2003 <br>
    !# @endin
    implicit none

    integer, intent (in) :: m, n

    real (kind = p_r8), dimension (n), intent (in) :: x
    real (kind = p_r8), dimension (m), intent (in) :: y

    integer, dimension (m), intent (OUT) :: indx, rc
    integer :: i, j

    out : do j = 1, m

      if (y(j) < x(1)) then
        indx(j) = 1
        rc(j) = 0   !PK bug
        cycle out
      end if
      if (y(j) > x(n)) then
        indx(j) = n + 1
        rc(j) = 1
        cycle out
      end if
      do i = 1, n
        if (y(j) == x(i)) then
          indx(j) = i
          rc(j) = 0
          cycle out
        end if
      end do
      do i = 1, n - 1
        if (y(j) > x(i) .and. y(j) < x(i + 1)) then
          indx(j) = i + 1
          rc(j) = 1
          cycle out
        end if
      end do

    end do out

  end subroutine bsrch


end module Mod_VerticalInterpolation
