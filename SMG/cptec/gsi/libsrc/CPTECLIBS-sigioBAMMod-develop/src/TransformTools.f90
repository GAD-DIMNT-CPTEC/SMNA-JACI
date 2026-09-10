module TransformTools
   use LegendreTransform
   use Fourier
   use ModConstants, only: r8
   implicit none
   type, extends(legendre)  :: spectrans
      real (kind=r8), pointer     :: coskx(:,:)
      real (kind=r8), pointer     :: sinkx(:,:)
      real (kind=r8), allocatable :: trigs(:)
      integer,        allocatable :: ifax(:)
      contains
      procedure, public :: getVort => getVort_
      procedure, public :: getConv => getConv_
      procedure, public :: DivgVortToUV => DivgVortToUV_
      procedure, public :: spec2grid => spec2grid_
      procedure, public :: grid2spec => grid2spec_
      procedure, public :: initTransform => initTransform_
      procedure, public :: destroyTransform => destroyTransform_
   end type

   contains
      subroutine initTransform_( self, Mend )
         class(specTrans), intent(inout) :: self
         integer,          intent(in   ) :: Mend
         integer :: iret
         integer :: i
         integer :: m
         real (kind=r8) :: ri2pi, ang

         iret = self%initLegendre( Mend )
         call initFFT(self%xMax, self%ifax, self%trigs)

         ri2pi=8.0_r8*ATAN(1.0_r8)/REAL(self%xMax,r8)

         allocate(self%coskx(self%xMax,self%Mend))
         allocate(self%sinkx(self%xMax,self%Mend))

         do i = 1, self%xMax
            do m = 1, self%Mend
               ang             = real((i-1)*m,r8)*ri2pi
               self%coskx(i,m) = cos(ang)
               self%sinkx(i,m) = sin(ang)
            enddo
         enddo

      end subroutine

      subroutine destroyTransform_( self )
         class(specTrans), intent(inout) :: self
         integer :: iret
         real (kind=r8) :: ri2pi, ang

         iret = self%destroyLegendre(  )
         
         deallocate(self%coskx)
         deallocate(self%sinkx)
         deallocate(self%ifax)
         deallocate(self%trigs)

      end subroutine

      subroutine spec2grid_ (self, qCoef, gGrid)
    
    
         class(specTrans), intent (in)    :: self       
         real (kind=r8), intent (inout) :: qcoef(:) ! MnWv2 or MnWv3
         real (kind=r8), intent (  out) :: ggrid(self%xMax,self%yMax)

         real (kind=r8), allocatable :: ap(:)
         real (kind=r8), allocatable :: am(:)

         integer :: i, j, jj, m
       
         allocate (ap(self%Mends))
         allocate (am(self%Mends))
       
         call self%transp(-1, qCoef, size(qCoef))
         
         do j = 1, self%yMaxHf
            call self%Spec2Four( qcoef, ap, am, j)
            jj = self%yMax - j + 1
            do i = 1, self%xMax
               gGrid(i,j ) = ap(1)
               gGrid(i,jj) = am(1)

               do m = 1, self%Mend
                  gGrid(i,j ) = gGrid(i,j ) + 2.0_r8 *    &
                               (ap(2*m+1)*self%coskx(i,m) - &
                                ap(2*m+2)*self%sinkx(i,m))

                  gGrid(i,jj) = gGrid(i,JJ) + 2.0_r8 *    &
                               (am(2*m+1)*self%coskx(i,m) - &
                                am(2*m+2)*self%sinkx(i,m))
               enddo

            enddo
         enddo

         call self%transp(+1, qCoef, size(qCoef))
    
         deallocate (ap, am)
    
      end subroutine spec2grid_

      subroutine grid2spec_(self, gGrid, qCoef)
         class(specTrans), intent (in   ) :: self
         real (kind=r8),   intent (in   ) :: gGrid(self%xMax,self%yMax)
         real (kind=r8),   intent (inout) :: qCoef(:) ! MnWv2 or MnWv3

         real (kind=r8), dimension (self%xmx) :: fn, fs, gn, gs
         integer :: inc
         integer :: i, j
       
         qCoef = 0.0_r8
         inc   = 1
       
         do j=1,self%yMaxHF
            fn = 0.0_r8
            fs = 0.0_r8
            do i=1,self%xMax
               fn(i) = gGrid(i,j)
               fs(i) = gGrid(i,self%yMax+1-j)
            enddo

            CALL fft991 (fn, gn, inc, self%xmx, self%xMax, 1, self%ifax, self%trigs, -1)
            CALL fft991 (fs, gs, inc, self%xmx, self%xMax, 1, self%ifax, self%trigs, -1)

            CALL self%SymAsy (fn, fs)
            CALL self%Four2Spec (fn, fs, j, size(qCoef), qCoef)
       
         enddo

         call self%transp(+1, qCoef, size(qCoef))

      end subroutine

      
      subroutine DivgVortToUV_ ( self, qDivg, qVort, qUvel, qVvel)
      
          ! Calculates Spectral Representation of Cosine-Weighted
          ! Wind Components from Spectral Representation of
          ! Vorticity and Divergence.
        
          ! qDivg Input:  Divergence (Spectral)
          ! qVort Input:  Vorticity  (Spectral)
          ! qUvel Output: Zonal Pseudo-Wind (Spectral)
          ! qVvel Output: Meridional Pseudo-Wind (Spectral)
        
          class(specTrans), intent(in) :: self
    
          real(kind=r8), intent(inout) :: qdivg(2,self%MnWv0) ! (2,mnwv0) = (MnWv2)
          real(kind=r8), intent(inout) :: qvort(2,self%MnWv0) ! (2,mnwv0) = (MnWv2)
          real(kind=r8), intent(  out) :: quvel(2,self%MnWv1) ! (2,mnwv1) = (MnWv3)
          real(kind=r8), intent(  out) :: qvvel(2,self%MnWv1) ! (2,mnwv1) = (MnWv3)
          
          integer :: mm, nn, l, l0, l0p, l0m, l1, l1p, Nmax
    
          qUvel=0.0_r8
          qVvel=0.0_r8
    
          qDivg(2,1:self%Mend1)=0.0_r8
          qVort(2,1:self%Mend1)=0.0_r8
    
          call self%transp ( -1, qDivg, size(qDivg))
          call self%transp ( -1, qVort, size(qVort))
       
        !cdir novector
          do mm=1,self%Mend1
    
            Nmax=self%Mend2+1-mm
    
            qUvel(1,mm) =  self%e1(mm)*qDivg(2,mm)
            qUvel(2,mm) = -self%e1(mm)*qDivg(1,mm)
    
            qVvel(1,mm) =  self%e1(mm)*qVort(2,mm)
            qVvel(2,mm) = -self%e1(mm)*qVort(1,mm)
    
            if (Nmax >= 3) then
              l=self%Mend1
              qUvel(1,mm) = qUvel(1,mm)+self%e0(mm+l)*qVort(1,mm+l)
              qUvel(2,mm) = qUvel(2,mm)+self%e0(mm+l)*qVort(2,mm+l)
    
              qVvel(1,mm) = qVvel(1,mm)-self%e0(mm+l)*qDivg(1,mm+l)
              qVvel(2,mm) = qVvel(2,mm)-self%e0(mm+l)*qDivg(2,mm+l)
            end if
    
            if (Nmax >= 4) then
              do nn=2,Nmax-2
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
              end do
            end if
    
            if (Nmax >= 3) then
    
              nn  = Nmax-1
              l0  = self%la0(mm,nn)
              l0m = self%la0(mm,nn-1)
              l1  = self%la1(mm,nn)
    
              qUvel(1,l1) = -self%e0(l1)*qVort(1,l0m)+self%e1(l0)*qDivg(2,l0)
              qUvel(2,l1) = -self%e0(l1)*qVort(2,l0m)-self%e1(l0)*qDivg(1,l0)
    
              qVvel(1,l1) =  self%e0(l1)*qDivg(1,l0m)+self%e1(l0)*qVort(2,l0)
              qVvel(2,l1) =  self%e0(l1)*qDivg(2,l0m)-self%e1(l0)*qVort(1,l0)
              
            end if
            if (Nmax >= 2) then
              nn=Nmax
              l0m=self%la0(mm,nn-1)
              l1 =self%la1(mm,nn)
              qUvel(1,l1)=-self%e0(l1)*qVort(1,l0m)
              qUvel(2,l1)=-self%e0(l1)*qVort(2,l0m)
    
              qVvel(1,l1)= self%e0(l1)*qDivg(1,l0m)
              qVvel(2,l1)= self%e0(l1)*qDivg(2,l0m)
            end if
          end do
        
          CALL self%transp (+1, qUvel, size(qUvel))
          CALL self%transp (+1, qVvel, size(qVvel))
          CALL self%transp (+1, qDivg, size(qDivg))
          CALL self%transp (+1, qVort, size(qVort))
      
      END SUBROUTINE DivgVortToUV_
      
      SUBROUTINE UVtoDivgVort_ (self, u, v, qdiv, qrot)
      
         class(specTrans), intent(in) :: self
      
         real (kind=r8), intent (in   ) :: u (self%xMax,self%yMax)
         real (kind=r8), intent (in   ) :: v (self%xMax,self%yMax)
         real (kind=r8), intent (  out) :: qrot (self%Mnwv2)
         real (kind=r8), intent (  out) :: qdiv (self%Mnwv2)
      
         integer :: i, j, k, inc
      
         real (kind=r8) :: coslat(self%yMaxHf)
         Real (kind=r8) :: un (self%xmx)
         Real (kind=r8) :: us (self%xmx)
         Real (kind=r8) :: vn (self%xmx)
         Real (kind=r8) :: vs (self%xmx)
         Real (kind=r8) :: gw (self%xmx)
      
         qrot=0.0_r8
         qdiv=0.0_r8
         inc=1
      
         coslat=SIN(self%colrad)
         DO j=1,self%yMaxHF
            DO i=1,self%xMax
               un(i) = u(i,j)*coslat(j)
               us(i) = u(i,self%yMax+1-j)*coslat(j)
               vn(i) = v(i,j)*coslat(j)
               vs(i) = v(i,self%yMax+1-j)*coslat(j)
            ENDDO
            DO i=self%xMax+1,self%xmx
               un(i) = 0.0_r8
               us(i) = 0.0_r8
               vn(i) = 0.0_r8
               vs(i) = 0.0_r8
            ENDDO
            CALL fft991 (un, gw, inc, self%xmx, self%xMax, 1, self%ifax, self%trigs, -1)
            CALL fft991 (us, gw, inc, self%xmx, self%xMax, 1, self%ifax, self%trigs, -1)
            CALL fft991 (vn, gw, inc, self%xmx, self%xMax, 1, self%ifax, self%trigs, -1)
            CALL fft991 (vs, gw, inc, self%xmx, self%xMax, 1, self%ifax, self%trigs, -1)
            CALL self%SymAsy (un, us)
            CALL self%SymAsy (vn, vs)
            CALL self%GetConv  (us, un, vs, vn, qdiv, j)
            CALL self%GetVort  (us, un, vs, vn, qrot, j)
         ENDDO
        
      !  Get_Conv computes convergence, 
      !           signal of qdiv must be changed to get divergence
      
         qdiv=-qdiv
      
      END SUBROUTINE UVtoDivgVort_

      subroutine GetVort_ (self,am, ap, bm, bp, fln, lat)
      
        !     calculates the spectral representations of the horizontal
        !     vorticity   of pseudo-vector fields from the fourier
        !     representations of the symmetric and anti-symmetric
        !     portions of the two individual fields.
        ! 
        !     argument(dimensions)         description
        ! 
        !     am(Imax+2,Kmax)       input: fourier representation of
        !                                  anti-symmetric portion of
        !                                  zonal pseudo-wind field at
        !                                  one gaussian latitude.
        !     ap(Imax+2,Kmax)       input: fourier representation of
        !                                  symmetric portion of zonal
        !                                  pseudo-wind field at one
        !                                  gaussian latitude.
        !     bm(Imax+2,Kmax)       input: fourier representation of
        !                                  anti-symmetric portion of
        !                                  meridional pseudo-wind field
        !                                  at one gaussian latitude.
        !     bp(Imax+2,Kmax)       input: fourier representation of
        !                                  symmetric portion of
        !                                  meridional pseudo-wind field
        !                                  at one gaussian latitude.
        !     fln(Mnwv2,Kmax)       input: spectral representation of
        !                                  the vorticity  of the global
        !                                  wind field. includes
        !                                  contributions from gaussian
        !                                  latitudes up to but not
        !                                  including current iteration
        !                                  of gaussian loop in calling
        !                                  routine.
        !                          output: spectral representation of
        !                                  the vorticity  of the global
        !                                  wind field. includes
        !                                  contributions from gaussian
        !                                  latitudes up to and
        !                                  including current iteration
        !                                  of gaussian loop in calling
        !                                  routine.
        !     lat                   input: current index of gaussian
        !                                  loop in calling routine.
       
        class(specTrans), intent (in   ) :: self
        integer,          intent (in   ) :: lat      
        real (kind=r8),   intent (in   ) :: am(self%xMax+2)
        real (kind=r8),   intent (in   ) :: ap(self%xMax+2)
        real (kind=r8),   intent (in   ) :: bm(self%xMax+2)
        real (kind=r8),   intent (in   ) :: bp(self%xMax+2)
        real (kind=r8),   intent (inout) :: fln(self%mnwv2)
      
        integer :: k, l, nn, mmax, mm, mn
      
        real (kind=r8), dimension (self%mnwv2) :: s 
       
        l=0
        do nn=1,self%mend1
           mmax=2*(self%mend2-nn)
           if (mod(nn-1,2) == 0) then
              do mm=1,mmax
                 l=l+1
                 s(l)=am(mm)
              end do
           else
              do mm=1,mmax
                 l=l+1
                 s(l)=ap(mm)
              end do
           end if
        end do

        do mn=1,self%mnwv0
           fln(2*mn-1)=fln(2*mn-1)+s(2*mn-1)*self%legDerNS(mn,lat)
           fln(2*mn  )=fln(2*mn  )+s(2*mn  )*self%legDerNS(mn,lat)
        end do

        l=0
        do nn=1,self%mend1
           mmax=2*(self%mend2-nn)
           if (mod(nn-1,2) == 0) then
              do mm=1,mmax,2
                 l=l+1
                 s(2*l-1) = -bp(mm+1)
                 s(2*l  ) = +bp(mm  )
              end do
           else
              do mm=1,mmax,2
                 l=l+1
                 s(2*l-1) = -bm(mm+1)
                 s(2*l  ) = +bm(mm  )
              end do
           end if
        end do

        do mn=1,self%mnwv0
           fln(2*mn-1)=fln(2*mn-1)+s(2*mn-1)*self%legDerEW(mn,lat)
           fln(2*mn  )=fln(2*mn  )+s(2*mn  )*self%legDerEW(mn,lat)
        end do
      
      end subroutine GetVort_

      subroutine GetConv_ (self, am, ap, bm, bp, fln, lat)
      
        !     calculates the spectral representations of the horizontal
        !     convergence of pseudo-vector fields from the fourier
        !     representations of the symmetric and anti-symmetric
        !     portions of the two individual fields.
        ! 
        !     argument(dimensions)         description
        ! 
        !     am(Imax+2,Kmax)       input: fourier representation of
        !                                  anti-symmetric portion of
        !                                  zonal pseudo-wind field at
        !                                  one gaussian latitude.
        !     ap(Imax+2,Kmax)       input: fourier representation of
        !                                  symmetric portion of zonal
        !                                  pseudo-wind field at one
        !                                  gaussian latitude.
        !     bm(Imax+2,Kmax)       input: fourier representation of
        !                                  anti-symmetric portion of
        !                                  meridional pseudo-wind field
        !                                  at one gaussian latitude.
        !     bp(Imax+2,Kmax)       input: fourier representation of
        !                                  symmetric portion of
        !                                  meridional pseudo-wind field
        !                                  at one gaussian latitude.
        !     fln(Mnwv2,Kmax)       input: spectral representation of
        !                                  the divergence of the global
        !                                  wind field. includes
        !                                  contributions from gaussian
        !                                  latitudes up to but not
        !                                  including current iteration
        !                                  of gaussian loop in calling
        !                                  routine.
        !                          output: spectral representation of
        !                                  the divergence of the global
        !                                  wind field. includes
        !                                  contributions from gaussian
        !                                  latitudes up to and
        !                                  including current iteration
        !                                  of gaussian loop in calling
        !                                  routine.
        !     lat                   input: current index of gaussian
        !                                  loop in calling routine.
        class(specTrans), intent (in   ) :: self
        integer,          intent (in   ) :: lat      
        real (kind=r8),   intent (in   ) :: am(self%xMax+2)
        real (kind=r8),   intent (in   ) :: ap(self%xMax+2)
        real (kind=r8),   intent (in   ) :: bm(self%xMax+2)
        real (kind=r8),   intent (in   ) :: bp(self%xMax+2)
        real (kind=r8),   intent (inout) :: fln(self%mnwv2)
      
        integer :: k, l, nn, mmax, mm, mn
      
        real (kind=r8), dimension (self%mnwv2) :: s 
      
           l=0
           do nn=1,self%mend1
              mmax=2*(self%mend2-nn)
              if (mod(nn-1,2) == 0) then
                 do mm=1,mmax
                    l=l+1
                    s(l) = bm(mm)
                 end do
              else
                 do mm=1,mmax
                    l = l + 1
                    s(l) = bp(mm)
                 end do
              end if
           end do

           do mn=1,self%mnwv0
              fln(2*mn-1)=fln(2*mn-1)+s(2*mn-1)*self%legDerS2F(mn,lat)
              fln(2*mn  )=fln(2*mn  )+s(2*mn  )*self%legDerS2F(mn,lat)
           end do

           l=0
           do nn=1,self%mend1
              mmax=2*(self%mend2-nn)
              if (mod(nn-1,2) == 0) then
                 do mm=1,mmax,2
                    l=l+1
                    s(2*l-1) = +ap(mm+1)
                    s(2*l  ) = -ap(mm  )
                 end do
              else
                 do mm=1,mmax,2
                    l=l+1
                    s(2*l-1) = +am(mm+1)
                    s(2*l  ) = -am(mm  )
                 end do
              end if
           end do
           
           do mn=1,self%mnwv0
              fln(2*mn-1)=fln(2*mn-1)+s(2*mn-1)*self%legDerEW(mn,lat)
              fln(2*mn  )=fln(2*mn  )+s(2*mn  )*self%legDerEW(mn,lat)
           end do
     
      end subroutine GetConv_


end module
