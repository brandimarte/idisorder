!  *******************************************************************  !
!  I-Disorder Fortran Code 2007-2014                                    !
!                                                                       !
!  Written by Alexandre Reily Rocha (reilya@ift.unesp.br),              !
!             Pedro Brandimarte (brandimarte@gmail.com) and             !
!             Alberto Torres (alberto.trj@gmail.com).                   !
!                                                                       !
!  Copyright (c), All Rights Reserved                                   !
!                                                                       !
!  This program is free software. You can redistribute it and/or        !
!  modify it under the terms of the GNU General Public License          !
!  (version 3 or later) as published by the Free Software Foundation    !
!  <http://fsf.org/>.                                                   !
!                                                                       !
!  This program is distributed in the hope that it will be useful, but  !
!  WITHOUT ANY WARRANTY, without even the implied warranty of           !
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU     !
!  General Public License for more details (file 'LICENSE_GPL'          !
!  distributed along with this program or at                            !
!  <http://www.gnu.org/licenses/gpl.html>).                             !
!  *******************************************************************  !
!                         MODULE idsrdr_recipes                         !
!  *******************************************************************  !
!  Description: numerical subroutines based on W. H. Press, S. A.       !
!  Teukolsky, W. T. Vetterling and B. P. Flannery. "Numerical Recipes:  !
!  The Art of Scientific Computing", Cambridge Univ. Press, 3rd         !
!  edition (2007).                                                      !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2013                                   !
!  *******************************************************************  !

MODULE idsrdr_recipes

  implicit none
  
  PUBLIC :: RECPSsimpson


CONTAINS


!  *******************************************************************  !
!                             RECPSsimpson                              !
!  *******************************************************************  !
!  Description: Extended Simpson's rule to estimate definite            !
!  integrals. Given the lower and upper limits of integration 'x1' and  !
!  'x2', and given 'n', this routine returns a equidistant grid         !
!  'x(1, ..., n)' and the array 'w(1, ..., n)' with corresponding       !
!  weights.                                                             !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2013                                   !
!  ****************************** INPUT ******************************  !
!  real*8 x1                   : Lower limit of integration             !
!  real*8 x2                   : Upper limit of integration             !
!  integer n                   : Number of grid points                  !
!  ***************************** OUTPUT ******************************  !
!  real*8 x(n)                 : Equidistant grid points                !
!  real*8 w(n)                 : Weights of grid points                 !
!  *******************************************************************  !
  subroutine RECPSsimpson (x1, x2, n, x, w)

!   Input variables.
    integer, intent(in) :: n
    real(8), intent(in) :: x1, x2
    real(8), intent(out) :: x(n), w(n)

!   Local variables.
    integer :: i
    real(8) :: d

!   Grid step.
    d = (x2 - x1) / (1.D0 * n - 1.D0)
    d = DABS(d)

    if (MOD(n,2) .ne. 0) then ! n is odd
!      First and last points.
       x(1) = x1
       x(n) = x2
       w(1) = d / 3.D0
       w(n) = d / 3.D0

!      Middle points alternates between 2/3 and 4/3.
       do i = 2,n-1,2
          x(i) = x1 + 1.D0 * (i - 1.D0) * d
          w(i) = 4.D0 * d / 3.D0
       enddo
       do i = 3,n-1,2
          x(i) = x1 + 1.D0 * (i - 1.0) * d
          w(i) = 2.D0 * d / 3.D0
       enddo

    else ! n is even
!      First and last two points.
       x(1) = x1
       x(n) = x2
       x(n - 1) = x2 - d
       w(1) = d / 3.D0
       w(n) = d / 2.D0
       w(n - 1) = 5.D0 * d / 6.D0

!      Middle points alternates between 2/3 and 4/3.
       do i = 2,n-1,2
          x(i) = x1 + 1.D0 * (i - 1.D0) * d
          w(i) = 4.D0 * d / 3.D0
       enddo
       do i = 3,n-2,2
          x(i) = x1 + 1.D0 * (i - 1.D0) * d
          w(i) = 2.D0 * d / 3.D0
       enddo

    endif

  end subroutine RECPSsimpson


!  *******************************************************************  !


END MODULE idsrdr_recipes
