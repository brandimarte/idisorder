!  *******************************************************************  !
!  I-Disorder Fortran Code                                              !
!                                                                       !
!  Written by Alexandre Reily Rocha and Pedro Brandimarte, 2007-2013    !
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
!                          MODULE fortran_std                           !
!  *******************************************************************  !
!  Description: software environment commonly used.                     !
!                                                                       !
!  Based on P. Lignelet, "Structures de Donnees en Fortran 90/95",      !
!  Masson  (1996).                                                      !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !

MODULE fortran_std


  implicit none

  integer, parameter :: fMAX_INT = HUGE(0) ! + infinity

  character(*), parameter :: INVITE = " ? "


CONTAINS


!  *******************************************************************  !
!                                 INCR                                  !
!  *******************************************************************  !
!  Description: increment of 1 the parameter 'K'.                       !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !
  integer function INCR (K)

!   Input variables.
    integer, intent (inout) :: k

    K = K + 1
    INCR = K


  end function INCR


!  *******************************************************************  !
!                                istrue                                 !
!  *******************************************************************  !
!  Description: returns 'istrue' if 'L_OPT' is present, otherwise       !
!  returns '.false.'.                                                   !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !
  logical function istrue (L_OPT)

!   Input variables.
    logical, intent(in), optional :: L_OPT

    if (present(L_OPT)) then
       istrue = L_OPT
    else
       istrue = .false.
    endif


  end function istrue


!  *******************************************************************  !
!                                spause                                 !
!  *******************************************************************  !
!  Description: make a pause at the screen. The execution resume by     !
!  user intervention (hitting the 'Return' key).                        !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !
  subroutine spause (MSGLEVEL)

!   Input variables.
    integer, intent(in), optional :: MSGLEVEL

    if (present(MSGLEVEL)) print*, "Press Return to resume..."
    write (*,"(a)", ADVANCE="NO") INVITE
    read*


  end subroutine spause


!  *******************************************************************  !
!                                differ                                 !
!  *******************************************************************  !
!  Description: test if 2 reals 'A' and 'B' can be considered as        !
!  differents, taking into account the precision of the type.           !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !
  logical function differ (A, B)

!   Input variables.
    double precision, intent(in) :: A, B

!   Local variables.
    integer, parameter :: IANY = 5
    double precision :: ABSA, THRESHOLD
    double precision, parameter :: EPSIL = IANY * EPSILON(A)
    double precision, parameter :: ZERO_MIN = IANY * TINY(A) / EPSILON(A)

    ABSA = ABS(A)

    if (ABSA <= ZERO_MIN) then
       THRESHOLD = EPSIL
    else
       THRESHOLD = EPSIL * ABSA
    endif

    differ = ABS(A-B) > THRESHOLD


  end function differ


!  *******************************************************************  !


END MODULE fortran_std
