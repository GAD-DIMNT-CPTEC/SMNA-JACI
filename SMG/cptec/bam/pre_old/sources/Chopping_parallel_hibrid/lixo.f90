 do k=1,nlevp
       do j=1,np
          do i=1,np
             pint(i,j,k) = hyai(k)*ps0 + hybi(k)*ps(i,j)
          end do
       end do
    end do
    !
    ! Set midpoint pressures and layer thicknesses
    !
    do k=1,nlev
       do j=1,np
          do i=1,np
             pmid(i,j,k) = hyam(k)*ps0 + hybm(k)*ps(i,j)
             pdel(i,j,k) = pint(i,j,k+1) - pint(i,j,k)
          end do
       end do
    end do
