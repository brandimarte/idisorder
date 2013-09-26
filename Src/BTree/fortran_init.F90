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
!                          MODULE fortran_init                          !
!  *******************************************************************  !
!  Description: management of identification and initialization of      !
!  objects.                                                             !
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

MODULE fortran_init

!
! Modules
!
  use fortran_std,     only:

  implicit none

  PRIVATE ! default is private

  character(*), parameter :: INITIAL = '*'

  integer :: ID_LAST = 0

  TYPE T_INIT; PRIVATE
     integer :: ID ! internal identification
     character :: DATA ! initialization indicator
  END TYPE T_INIT

  TYPE(T_INIT), parameter, PUBLIC :: INIT_0 = T_INIT(0,INITIAL)

  PUBLIC :: T_INIT, IS_INIT, INITIALIZE, TEST_INIT, GET_ID


CONTAINS


!  *******************************************************************  !
!                                IS_INIT                                !
!  *******************************************************************  !
!  Description: check if 'INIT' is initialized.                         !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !
  logical function IS_INIT (INIT)

!   Input variables.
    TYPE(T_INIT), intent(in) :: INIT

    IS_INIT = INIT%DATA==INITIAL .and. INIT%ID>0 .and. INIT%ID<=ID_LAST


  end function IS_INIT


!  *******************************************************************  !
!                              INITIALIZE                               !
!  *******************************************************************  !
!  Description: initialize 'INIT'.                                      !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !
  subroutine INITIALIZE (INIT)

!
! Modules
!
    use fortran_std,     only: INCR

!   Input variables.
    TYPE(T_INIT), intent(out) :: INIT

    INIT = T_INIT(INCR(ID_LAST),INITIAL)


  end subroutine INITIALIZE


!  *******************************************************************  !
!                               TEST_INIT                               !
!  *******************************************************************  !
!  Description: check if 'INIT' is initialized and, if not,             !
!  initializes it.                                                      !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !
  subroutine TEST_INIT (INIT)

!   Input variables.
    TYPE(T_INIT), intent(inout) :: INIT

    if (.not. IS_INIT(INIT)) call INITIALIZE (INIT)


  end subroutine TEST_INIT


!  *******************************************************************  !
!                                GET_ID                                 !
!  *******************************************************************  !
!  Description: returns the 'ID' field from 'INIT'.                     !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !
  integer function GET_ID (INIT)

!   Input variables.
    TYPE(T_INIT), intent(in) :: INIT

    GET_ID = INIT%ID


  end function GET_ID


!  *******************************************************************  !


END MODULE fortran_init
