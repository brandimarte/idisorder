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
!                         MODULE fortran_except                         !
!  *******************************************************************  !
!  Description: management of failures/anomaly messages.                !
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

MODULE fortran_except

!
! Modules
!
  use fortran_std,     only:
  use fortran_init

  implicit none

  integer, parameter :: FAIL_RESEARCH = -1
  character(*), parameter :: MSG_FAIL_S = "Research failure"
  integer, parameter :: CONSTRAINT_ERROR = 1
  character(*), parameter :: MSG_CONST_ERR = "Parameter out of bounds"


CONTAINS


!  *******************************************************************  !
!                           EXCEPTION_MESSAGE                           !
!  *******************************************************************  !
!  Description: returns the 'EXCEPTION' via 'IOSTAT' if it is present,  !
!  otherwise displays 'MSG'.                                            !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  TYPE(T_INIT) OBJ        : Object where the exception occoured        !
!  integer EXCEPTION       : Number corresponding to an exception       !
!  character(*) MSG        : Message of exception                       !
!  ***************************** OUTPUT ******************************  !
!  integer IOSTAT          : Returns 'EXCEPTION' error code if it is    !
!                            present, otherwise returns 'fMAX_INT'      !
!  *******************************************************************  !
  subroutine EXCEPTION_MESSAGE (OBJ, IOSTAT, EXCEPTION, MSG)

!
! Modules
!
    use fortran_std,     only: fMAX_INT, spause

!   Input variables.
    TYPE(T_INIT), intent (in) :: OBJ
    integer, intent (out), optional :: IOSTAT
    integer, intent (in), optional :: EXCEPTION
    character(*), intent (in) :: MSG

    If (present(IOSTAT)) Then
       if (present(EXCEPTION)) then
          IOSTAT = EXCEPTION
       else
          IOSTAT = fMAX_INT
       endif
    Else
       write (*,"()") ! ensure result displayed in column 1
       print*, "*** OBJECT: ", GET_ID(OBJ), "  ***  ", MSG
       call spause (1)
    EndIf


  end subroutine EXCEPTION_MESSAGE


!  *******************************************************************  !
!                              IN_INTERVAL                              !
!  *******************************************************************  !
!  Description: check if an integer 'I' belongs to the interval         !
!  '[MIN,MAX]'. If 'MIN' is missing, then check if 'I <= MAX'. If       !
!  'MAX' is missing then check if 'MIN <= I'. If 'I' doesn't belong to  !
!  the interval it returns via 'IOSTAT' the exception                   !
!  'CONSTRAINT_ERROR'.                                                  !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  integer I               : Variable to be checked                     !
!  integer MIN             : Lower limit                                !
!  integer MAX             : Upper limit                                !
!  ***************************** OUTPUT ******************************  !
!  integer IOSTAT          : Returns 0 if succeeds, else returns a      !
!                            Fortran error code                         !
!  *******************************************************************  !
  logical function IN_INTERVAL (I, MIN, MAX, IOSTAT)

!   Input variables.
    integer, intent (in) :: I
    integer, intent (in), optional :: MIN, MAX
    integer, intent (out), optional :: IOSTAT

    IN_INTERVAL = .true.

    if (present(MIN)) IN_INTERVAL = MIN <= I
    if (present(MAX)) IN_INTERVAL = IN_INTERVAL .and. I<=MAX
    if (IN_INTERVAL) then
       if (present(IOSTAT)) IOSTAT = 0
    else
       call EXCEPTION_MESSAGE (INIT_0, IOSTAT, CONSTRAINT_ERROR,        &
            MSG_CONST_ERR)
    endif


  end function IN_INTERVAL


!  *******************************************************************  !
!                              ASSIGN_MSG                               !
!  *******************************************************************  !
!  Description: display a message (nonblocking) when attempting an      !
!  assignment for one type of 'OBJ_ID' which is forbidden (by           !
!  redefining).                                                         !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  TYPE(T_INIT) OBJ        : Object where the exception occoured        !
!  *******************************************************************  !
  subroutine ASSIGN_MSG (OBJ)

!   Input variables.
    TYPE(T_INIT), intent (in) :: OBJ

!   Local variables.
    TYPE(T_INIT) :: WOBJ

    if (IS_INIT(OBJ)) then
       WOBJ = OBJ
    else
       WOBJ = INIT_0
    endif
    call EXCEPTION_MESSAGE (WOBJ,                                       &
         MSG="Prohibted assignment for this type")


  end subroutine ASSIGN_MSG


!  *******************************************************************  !


END MODULE fortran_except
