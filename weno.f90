module weno

  implicit none

  private
  public :: weno5

contains

  subroutine weno5(gradphi, phi, advvel, &
       axis, bc0, bcn, &
       isize, jsize, ksize, &
       dx, dy, dz)
  
    implicit none
  
    integer, intent(in) :: axis
    integer, intent(in) :: bc0, bcn
    integer, intent(in) :: isize, jsize, ksize
    real(kind=8), intent(in) :: dx, dy, dz
    real(kind=8), dimension(isize, jsize, ksize), intent(in) :: phi
    real(kind=8), dimension(isize, jsize, ksize), intent(in) :: advvel
    
    real(kind=8), dimension(isize, jsize, ksize), intent(inout) :: gradphi
    
    integer :: i, j, k
    integer :: istep, jstep, kstep
    integer :: istart, jstart, kstart, iend, jend, kend
    integer :: im1, im2, im3, ip1, ip2
    integer :: jm1, jm2, jm3, jp1, jp2
    integer :: km1, km2, km3, kp1, kp2
    real(kind=8), parameter :: e = 1.0d-16
    real(kind=8), parameter :: zero = 0.d0, &
         one = 1.d0, &
         two = 2.d0, &
         three = 3.d0, &
         four = 4.d0, &
         five = 5.d0, &
         six = 6.d0, &
         seven = 7.d0, &
         ten = 10.d0, &
         eleven = 11.d0, &
         twelve = 12.d0, &
         thirteen = 13.d0
    
    real(kind=8) :: q1, q2, q3, q4, q5
    real(kind=8) :: a1, a2, a3
    real(kind=8) :: w1, w2, w3
    real(kind=8) :: is1, is2, is3
    real(kind=8) :: dsign
    real(kind=8) :: deltax
  
    !! Defaults
    istart = 1
    iend = isize
    jstart = 1
    jend = jsize
    kstart = 1
    kend = ksize
    
    istep = 0
    jstep = 0
    kstep = 0
    
    if (axis==1) then
       deltax = dx
    
       istart = 4
       iend = isize - 3
       istep = 1
    elseif (axis==2) then
       deltax = dy
    
       jstart = 4
       jend = jsize - 3
       jstep = 1
    elseif (axis==3) then
       deltax = dz
    
       kstart = 4
       kend = ksize - 3
       kstep = 1
    else
       print *, "ERROR: Invalid axis passed to WENO5"
       stop
    endif
  
    do k = kstart, kend
       do j = jstart, jend
          !! Note, if axis==2 and y is stretched, need to set deltax here
          do i = istart, iend
             if (advvel(i, j, k) > zero) then
                dsign = one
             
                istep = istep
                jstep = jstep
                kstep = kstep
             elseif (advvel(i, j, k) < zero) then
                dsign = -one
             
                istep = -istep
                jstep = -jstep
                kstep = -kstep
             else
                gradphi(i, j, k) = zero
                cycle
             endif
             
             im1 = i - 1 * istep
             im2 = i - 2 * istep
             im3 = i - 3 * istep
             ip1 = i + 1 * istep
             ip2 = i + 2 * istep
             
             jm1 = j - 1 * jstep
             jm2 = j - 2 * jstep
             jm3 = j - 3 * jstep
             jp1 = j + 1 * jstep
             jp2 = j + 2 * jstep
             
             km1 = k - 1 * kstep
             km2 = k - 2 * kstep
             km3 = k - 3 * kstep
             kp1 = k + 1 * kstep
             kp2 = k + 2 * kstep
             q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
             q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
             q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
             q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
             q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
             is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                  + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
             is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                  + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
             is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                  + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
             a1 = one / (e + is1)**2 / ten
             a2 = six / (e + is2)**2 / ten
             a3 = three / (e + is3)**2 / ten
             
             w1 = a1 / (a1 + a2 + a3)
             w2 = a2 / (a1 + a2 + a3)
             w3 = a3 / (a1 + a2 + a3)
             gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                  + w2 * (-q2 + five * q3 + two * q4) &
                  + w3 * (two * q3 + five * q4 - q5)
             gradphi(i, j, k) = gradphi(i, j, k) / six
          enddo
  
          if (axis==1) then
             jm1 = j
             jm2 = j
             jm3 = j
             jp1 = j
             jp2 = j
          
             km1 = k
             km2 = k
             km3 = k
             kp1 = k
             kp2 = k
          
             if ((bc0==0).and.(bcn==0)) then
                i = 1
                if (advvel(i, j, k) == zero) then
                   gradphi(i, j, k) = zero
                else
                   if (advvel(i, j, k) > zero) then
                      dsign = one
          
                      im1 = isize
                      im2 = isize - 1
                      im3 = isize - 2
                      ip1 = i + 1
                      ip2 = i + 2
                   else
                      dsign = -one
          
                      im1 = i + 1
                      im2 = i + 2
                      im3 = i + 3
                      ip1 = isize
                      ip2 = isize - 1
                   endif
                   q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
                   q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
                   q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
                   q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
                   q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
                   is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                        + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
                   is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                        + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
                   is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                        + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
                   a1 = one / (e + is1)**2 / ten
                   a2 = six / (e + is2)**2 / ten
                   a3 = three / (e + is3)**2 / ten
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                        + w2 * (-q2 + five * q3 + two * q4) &
                        + w3 * (two * q3 + five * q4 - q5)
                   gradphi(i, j, k) = gradphi(i, j, k) / six
                endif
          
                i = 2
                if (advvel(i, j, k) == zero) then
                   gradphi(i, j, k) = zero
                else
                   if (advvel(i, j, k) > zero) then
                      dsign = one
          
                      im1 = i - 1
                      im2 = isize
                      im3 = isize - 1
                      ip1 = i + 1
                      ip2 = i + 2
                   else
                      dsign = -one
          
                      im1 = i + 1
                      im2 = i + 2
                      im3 = i + 3
                      ip1 = i - 1
                      ip2 = isize
                   endif
                   q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
                   q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
                   q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
                   q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
                   q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
                   is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                        + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
                   is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                        + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
                   is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                        + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
                   a1 = one / (e + is1)**2 / ten
                   a2 = six / (e + is2)**2 / ten
                   a3 = three / (e + is3)**2 / ten
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                        + w2 * (-q2 + five * q3 + two * q4) &
                        + w3 * (two * q3 + five * q4 - q5)
                   gradphi(i, j, k) = gradphi(i, j, k) / six
                endif
          
                i = 3
                if (advvel(i, j, k) == zero) then
                   gradphi(i, j, k) = zero
                else
                   if (advvel(i, j, k) > zero) then
                      dsign = one
          
                      im1 = i - 1
                      im2 = i - 2
                      im3 = isize
                      ip1 = i + 1
                      ip2 = i + 2
                   else
                      dsign = -one
          
                      im1 = i + 1
                      im2 = i + 2
                      im3 = i + 3
                      ip1 = i - 1
                      ip2 = i - 2
                   endif
                   q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
                   q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
                   q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
                   q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
                   q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
                   is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                        + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
                   is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                        + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
                   is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                        + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
                   a1 = one / (e + is1)**2 / ten
                   a2 = six / (e + is2)**2 / ten
                   a3 = three / (e + is3)**2 / ten
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                        + w2 * (-q2 + five * q3 + two * q4) &
                        + w3 * (two * q3 + five * q4 - q5)
                   gradphi(i, j, k) = gradphi(i, j, k) / six
                endif
          
                i = isize
                if (advvel(i, j, k)==zero) then
                   gradphi(i, j, k) = zero
                else
                   if (advvel(i, j, k) > zero) then
                      dsign = one
          
                      im1 = i - 1
                      im2 = i - 2
                      im3 = i - 3
                      ip1 = 1
                      ip2 = 2
                   else
                      dsign = -one
          
                      im1 = 1
                      im2 = 2
                      im3 = 3
                      ip1 = i - 1
                      ip2 = i - 2
                   endif
                   q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
                   q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
                   q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
                   q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
                   q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
                   is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                        + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
                   is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                        + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
                   is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                        + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
                   a1 = one / (e + is1)**2 / ten
                   a2 = six / (e + is2)**2 / ten
                   a3 = three / (e + is3)**2 / ten
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                        + w2 * (-q2 + five * q3 + two * q4) &
                        + w3 * (two * q3 + five * q4 - q5)
                   gradphi(i, j, k) = gradphi(i, j, k) / six
                endif
          
                i = isize - 1
                if (advvel(i, j, k) == zero) then
                   gradphi(i, j, k) = zero
                else
                   if (advvel(i, j, k) > zero) then
                      dsign = one
          
                      im1 = i - 1
                      im2 = i - 2
                      im3 = i - 3
                      ip1 = i + 1
                      ip2 = 1
                   else
                      dsign = -one
          
                      im1 = i + 1
                      im2 = 1
                      im3 = 2
                      ip1 = i - 1
                      ip2 = i - 2
                   endif
                   q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
                   q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
                   q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
                   q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
                   q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
                   is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                        + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
                   is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                        + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
                   is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                        + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
                   a1 = one / (e + is1)**2 / ten
                   a2 = six / (e + is2)**2 / ten
                   a3 = three / (e + is3)**2 / ten
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                        + w2 * (-q2 + five * q3 + two * q4) &
                        + w3 * (two * q3 + five * q4 - q5)
                   gradphi(i, j, k) = gradphi(i, j, k) / six
                endif
          
                i = isize - 2
                if (advvel(i, j, k) == zero) then
                   gradphi(i, j, k) = zero
                else
                   if (advvel(i, j, k) > zero) then
                      dsign = one
          
                      im1 = i - 1
                      im2 = i - 2
                      im3 = i - 3
                      ip1 = i + 1
                      ip2 = i + 2
                   else
                      dsign = -one
          
                      im1 = i + 1
                      im2 = i + 2
                      im3 = 1
                      ip1 = i - 1
                      ip2 = i - 2
                   endif
                   q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
                   q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
                   q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
                   q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
                   q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
                   is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                        + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
                   is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                        + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
                   is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                        + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
                   a1 = one / (e + is1)**2 / ten
                   a2 = six / (e + is2)**2 / ten
                   a3 = three / (e + is3)**2 / ten
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                        + w2 * (-q2 + five * q3 + two * q4) &
                        + w3 * (two * q3 + five * q4 - q5)
                   gradphi(i, j, k) = gradphi(i, j, k) / six
                endif
             else
                !! Use second order
                i = 1
                if (bc0==1) then ! Zero grad
                   gradphi(i, j, k) = zero
                else ! Fixed value
                   gradphi(i, j, k) = (phi(i + 1, j, k) - phi(i, j, k)) / dx
                endif
                do i = 2, 3
                   gradphi(i, j, k) = (phi(i + 1, j, k) - phi(i - 1, j, k)) / (two * dx)
                enddo
          
                do i = isize - 2, isize - 1
                   gradphi(i, j, k) = (phi(i + 1, j, k) - phi(i - 1, j, k)) / (two * dx)
                enddo
                i = isize
                if (bcn==1) then ! Zero grad
                   gradphi(i, j, k) = zero
                else
                   gradphi(i, j, k) = (phi(i, j, k) - phi(i - 1, j, k)) / dx
                endif
             endif
          endif
       enddo
  
       if (axis==2) then
          km1 = k
          km2 = k
          km3 = k
          kp1 = k
          kp2 = k
       
          if ((bc0==0).and.(bcn==0)) then
             j = 1
       
             do i = 1, isize
                im1 = i
                im2 = i
                im3 = i
                ip1 = i
                ip2 = i
       
                if (advvel(i, j, k)==zero) then
                   gradphi(i, j, k) = zero
                else
                   if (advvel(i, j, k) > zero) then
                      dsign = one
       
                      jm1 = jsize
                      jm2 = jsize - 1
                      jm3 = jsize - 2
                      jp1 = j + 1
                      jp2 = j + 2
                   else
                      dsign = -one
       
                      jm1 = j + 1
                      jm2 = j + 2
                      jm3 = j + 3
                      jp1 = jsize
                      jp2 = jsize - 1
                   endif
                   q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
                   q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
                   q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
                   q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
                   q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
                   is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                        + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
                   is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                        + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
                   is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                        + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
                   a1 = one / (e + is1)**2 / ten
                   a2 = six / (e + is2)**2 / ten
                   a3 = three / (e + is3)**2 / ten
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                        + w2 * (-q2 + five * q3 + two * q4) &
                        + w3 * (two * q3 + five * q4 - q5)
                   gradphi(i, j, k) = gradphi(i, j, k) / six
                endif
       
                j = 2
                if (advvel(i, j, k)==zero) then
                   gradphi(i, j, k) = zero
                else
                   if (advvel(i, j, k) > zero) then
                      dsign = one
       
                      jm1 = j - 1
                      jm2 = jsize
                      jm3 = jsize - 1
                      jp1 = j + 1
                      jp2 = j + 2
                   else
                      dsign = -one
       
                      jm1 = j + 1
                      jm2 = j + 2
                      jm3 = j + 3
                      jp1 = j - 1
                      jp2 = jsize
                   endif
                   q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
                   q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
                   q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
                   q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
                   q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
                   is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                        + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
                   is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                        + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
                   is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                        + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
                   a1 = one / (e + is1)**2 / ten
                   a2 = six / (e + is2)**2 / ten
                   a3 = three / (e + is3)**2 / ten
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                        + w2 * (-q2 + five * q3 + two * q4) &
                        + w3 * (two * q3 + five * q4 - q5)
                   gradphi(i, j, k) = gradphi(i, j, k) / six
                endif
       
                j = 3
                if (advvel(i, j, k)==zero) then
                   gradphi(i, j, k) = zero
                else
                   if (advvel(i, j, k) > zero) then
                      dsign = one
       
                      jm1 = j - 1
                      jm2 = j - 2
                      jm3 = jsize
                      jp1 = j + 1
                      jp2 = j + 2
                   else
                      dsign = -one
       
                      jm1 = j + 1
                      jm2 = j + 2
                      jm3 = j + 3
                      jp1 = j - 1
                      jp2 = j - 2
                   endif
                   q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
                   q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
                   q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
                   q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
                   q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
                   is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                        + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
                   is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                        + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
                   is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                        + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
                   a1 = one / (e + is1)**2 / ten
                   a2 = six / (e + is2)**2 / ten
                   a3 = three / (e + is3)**2 / ten
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                        + w2 * (-q2 + five * q3 + two * q4) &
                        + w3 * (two * q3 + five * q4 - q5)
                   gradphi(i, j, k) = gradphi(i, j, k) / six
                endif
       
                j = jsize
                if (advvel(i, j, k) == zero) then
                   gradphi(i, j, k) = zero
                else
                   if (advvel(i, j, k) > zero) then
                      dsign = one
       
                      jm1 = j - 1
                      jm2 = j - 2
                      jm3 = j - 3
                      jp1 = j
                      jp2 = j
                   else
                      dsign = -one
       
                      jm1 = 1
                      jm2 = 2
                      jm3 = 3
                      jp1 = j - 1
                      jp2 = j - 2
                   endif
                   q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
                   q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
                   q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
                   q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
                   q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
                   is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                        + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
                   is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                        + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
                   is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                        + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
                   a1 = one / (e + is1)**2 / ten
                   a2 = six / (e + is2)**2 / ten
                   a3 = three / (e + is3)**2 / ten
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                        + w2 * (-q2 + five * q3 + two * q4) &
                        + w3 * (two * q3 + five * q4 - q5)
                   gradphi(i, j, k) = gradphi(i, j, k) / six
                endif
       
                j = jsize - 1
                if (advvel(i, j, k)==zero) then
                   gradphi(i, j, k) = zero
                else
                   if (advvel(i, j, k) > zero) then
                      dsign = one
       
                      jm1 = j - 1
                      jm2 = j - 2
                      jm3 = j - 3
                      jp1 = j + 1
                      jp2 = 1
                   else
                      dsign = -one
       
                      jm1 = j + 1
                      jm2 = 1
                      jm3 = 2
                      jp1 = j - 1
                      jp2 = j - 2
                   endif
                   q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
                   q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
                   q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
                   q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
                   q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
                   is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                        + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
                   is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                        + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
                   is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                        + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
                   a1 = one / (e + is1)**2 / ten
                   a2 = six / (e + is2)**2 / ten
                   a3 = three / (e + is3)**2 / ten
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                        + w2 * (-q2 + five * q3 + two * q4) &
                        + w3 * (two * q3 + five * q4 - q5)
                   gradphi(i, j, k) = gradphi(i, j, k) / six
                endif
       
                j = jsize - 2
                if (advvel(i, j, k)==zero) then
                   gradphi(i, j, k) = zero
                else
                   if (advvel(i, j, k) > zero) then
                      dsign = one
       
                      jm1 = j - 1
                      jm2 = j - 2
                      jm3 = j - 3
                      jp1 = j + 1
                      jp2 = j + 2
                   else
                      dsign = -one
       
                      jm1 = j + 1
                      jm2 = j + 2
                      jm3 = 1
                      jp1 = j - 1
                      jp2 = j - 2
                   endif
                   q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
                   q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
                   q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
                   q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
                   q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
                   is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                        + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
                   is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                        + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
                   is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                        + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
                   a1 = one / (e + is1)**2 / ten
                   a2 = six / (e + is2)**2 / ten
                   a3 = three / (e + is3)**2 / ten
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                        + w2 * (-q2 + five * q3 + two * q4) &
                        + w3 * (two * q3 + five * q4 - q5)
                   gradphi(i, j, k) = gradphi(i, j, k) / six
                endif
             enddo
          else
             do i = 1, isize
                !! Use second order
                j = 1
                if (bc0==1) then ! Zero grad
                   gradphi(i, j, k) = zero
                else ! Fixed value
                   gradphi(i, j, k) = (phi(i, j + 1, k) - phi(i, j, k)) / dy
                endif
                do j = 2, 3
                   gradphi(i, j, k) = (phi(i, j + 1, k) - phi(i, j - 1, k)) / (two * dy)
                enddo
       
                do j = jsize - 2, jsize - 1
                   gradphi(i, j, k) = (phi(i, j + 1, k) - phi(i, j - 1, k)) / (two * dy)
                enddo
                j = jsize
                if (bcn==1) then ! Zero grad
                   gradphi(i, j, k) = zero
                else
                   gradphi(i, j, k) = (phi(i, j, k) - phi(i, j - 1, k)) / dy
                endif
             enddo
          endif
       endif
    enddo
  
    if (axis==3) then
       if ((bc0==0).and.(bcn==0)) then
          do j = 1, jsize
             do i = 1, isize
                jm1 = j
                jm2 = j
                jm3 = j
                jp1 = j
                jp2 = j
    
                im1 = i
                im2 = i
                im3 = i
                ip1 = i
                ip2 = i
    
                k = 1
                if (advvel(i, j, k)==zero) then
                   gradphi(i, j, k) = zero
                else
                   if (advvel(i, j, k) > zero) then
                      dsign = one
    
                      km1 = ksize
                      km2 = ksize - 1
                      km3 = ksize - 2
                      kp1 = k + 1
                      kp2 = k + 2
                   else
                      dsign = -one
    
                      km1 = k + 1
                      km2 = k + 2
                      km3 = k + 3
                      kp1 = ksize
                      kp2 = ksize - 1
                   endif
                   q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
                   q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
                   q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
                   q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
                   q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
                   is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                        + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
                   is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                        + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
                   is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                        + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
                   a1 = one / (e + is1)**2 / ten
                   a2 = six / (e + is2)**2 / ten
                   a3 = three / (e + is3)**2 / ten
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                        + w2 * (-q2 + five * q3 + two * q4) &
                        + w3 * (two * q3 + five * q4 - q5)
                   gradphi(i, j, k) = gradphi(i, j, k) / six
                endif
    
                k = 2
                if (advvel(i, j, k)==zero) then
                   gradphi(i, j, k) = zero
                else
                   if (advvel(i, j, k) > zero) then
                      dsign = one
    
                      km1 = k - 1
                      km2 = ksize
                      km3 = ksize - 1
                      kp1 = k + 1
                      kp2 = k + 2
                   else
                      dsign = -one
    
                      km1 = k + 1
                      km2 = k + 2
                      km3 = k + 3
                      kp1 = k - 1
                      kp2 = ksize
                   endif
                   q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
                   q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
                   q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
                   q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
                   q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
                   is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                        + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
                   is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                        + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
                   is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                        + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
                   a1 = one / (e + is1)**2 / ten
                   a2 = six / (e + is2)**2 / ten
                   a3 = three / (e + is3)**2 / ten
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                        + w2 * (-q2 + five * q3 + two * q4) &
                        + w3 * (two * q3 + five * q4 - q5)
                   gradphi(i, j, k) = gradphi(i, j, k) / six
                endif
    
                k = 3
                if (advvel(i, j, k)==zero) then
                   gradphi(i, j, k) = zero
                else
                   if (advvel(i, j, k) > zero) then
                      dsign = one
    
                      km1 = k - 1
                      km2 = k - 2
                      km3 = ksize
                      kp1 = k + 1
                      kp2 = k + 2
                   else
                      dsign = -one
    
                      km1 = k + 1
                      km2 = k + 2
                      km3 = k + 3
                      kp1 = k - 1
                      kp2 = k - 2
                   endif
                   q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
                   q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
                   q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
                   q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
                   q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
                   is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                        + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
                   is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                        + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
                   is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                        + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
                   a1 = one / (e + is1)**2 / ten
                   a2 = six / (e + is2)**2 / ten
                   a3 = three / (e + is3)**2 / ten
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                        + w2 * (-q2 + five * q3 + two * q4) &
                        + w3 * (two * q3 + five * q4 - q5)
                   gradphi(i, j, k) = gradphi(i, j, k) / six
                endif
    
                k = ksize
                if (advvel(i, j, k) == zero) then
                   gradphi(i, j, k) = zero
                else
                   if (advvel(i, j, k) > zero) then
                      dsign = one
    
                      km1 = k - 1
                      km2 = k - 2
                      km3 = k - 3
                      kp1 = 1
                      kp2 = 2
                   else
                      dsign = -one
    
                      km1 = 1
                      km2 = 2
                      km3 = 3
                      kp1 = k - 1
                      kp2 = k - 2
                   endif
                   q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
                   q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
                   q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
                   q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
                   q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
                   is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                        + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
                   is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                        + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
                   is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                        + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
                   a1 = one / (e + is1)**2 / ten
                   a2 = six / (e + is2)**2 / ten
                   a3 = three / (e + is3)**2 / ten
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                        + w2 * (-q2 + five * q3 + two * q4) &
                        + w3 * (two * q3 + five * q4 - q5)
                   gradphi(i, j, k) = gradphi(i, j, k) / six
                endif
    
                k = ksize - 1
                if (advvel(i, j, k) == zero) then
                   gradphi(i, j, k) = zero
                else
                   if (advvel(i, j, k) > zero) then
                      dsign = one
    
                      km1 = k - 1
                      km2 = k - 2
                      km3 = k - 3
                      kp1 = k + 1
                      kp2 = 1
                   else
                      dsign = -one
    
                      km1 = k + 1
                      km2 = 1
                      km3 = 2
                      kp1 = k - 1
                      kp2 = k - 2
                   endif
                   q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
                   q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
                   q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
                   q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
                   q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
                   is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                        + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
                   is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                        + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
                   is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                        + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
                   a1 = one / (e + is1)**2 / ten
                   a2 = six / (e + is2)**2 / ten
                   a3 = three / (e + is3)**2 / ten
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                        + w2 * (-q2 + five * q3 + two * q4) &
                        + w3 * (two * q3 + five * q4 - q5)
                   gradphi(i, j, k) = gradphi(i, j, k) / six
                endif
    
                k = ksize - 2
                if (advvel(i, j, k) == zero) then
                   gradphi(i, j, k) = zero
                else
                   if (advvel(i, j, k) > zero) then
                      dsign = one
    
                      km1 = k - 1
                      km2 = k - 2
                      km3 = k - 3
                      kp1 = k + 1
                      kp2 = k + 2
                   else
                      dsign = -one
    
                      km1 = k + 1
                      km2 = k + 2
                      km3 = 1
                      kp1 = k - 1
                      kp2 = k - 2
                   endif
                   q1 = dsign * (phi(im2, jm2, km2) - phi(im3, jm3, km3)) / deltax
                   q2 = dsign * (phi(im1, jm1, km1) - phi(im2, jm2, km2)) / deltax
                   q3 = dsign * (phi(i, j, k) - phi(im1, jm1, km1)) / deltax
                   q4 = dsign * (phi(ip1, jp1, kp1) - phi(i, j, k)) / deltax
                   q5 = dsign * (phi(ip2, jp2, kp2) - phi(ip1, jp1, kp1)) / deltax
                   is1 = (thirteen / twelve) * (phi(im2,jm2,km2) - two * phi(im1,jm1,km1) + phi(i,j,k))**2 &
                        + (phi(im2,jm2,km2) - four * phi(im1,jm1,km1) + three * phi(i,j,k))**2 / four
                   is2 = (thirteen / twelve) * (phi(im1,jm1,km1) - two * phi(i,j,k) + phi(ip1,jp1,kp1))**2 &
                        + (phi(im1,jm1,km1) - phi(ip1,jp1,kp1))**2 / four
                   is3 = (thirteen / twelve) * (phi(i,j,k) - two * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 &
                        + (three * phi(i,j,k) - four * phi(ip1,jp1,kp1) + phi(ip2,jp2,kp2))**2 / four
                   a1 = one / (e + is1)**2 / ten
                   a2 = six / (e + is2)**2 / ten
                   a3 = three / (e + is3)**2 / ten
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   gradphi(i, j, k) = w1 * (two * q1 - seven * q2 + eleven * q3) &
                        + w2 * (-q2 + five * q3 + two * q4) &
                        + w3 * (two * q3 + five * q4 - q5)
                   gradphi(i, j, k) = gradphi(i, j, k) / six
                endif
             enddo
          enddo
       else
          do j = 1, jsize
             do i = 1, isize
                !! Use second order
                k = 1
                if (bc0==1) then ! Zero grad
                   gradphi(i, j, k) = zero
                else ! Fixed value
                   gradphi(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k)) / dz
                endif
                do k = 2, 3
                   gradphi(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k - 1)) / (two * dz)
                enddo
    
                do k = ksize - 2, ksize - 1
                   gradphi(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k - 1)) / (two * dz)
                enddo
                k = ksize
                if (bcn==1) then ! Zero grad
                   gradphi(i, j, k) = zero
                else
                   gradphi(i, j, k) = (phi(i, j, k) - phi(i, j, k - 1)) / dz
                endif
             enddo
          enddo
       endif
    endif
  
  endsubroutine weno5

endmodule weno
