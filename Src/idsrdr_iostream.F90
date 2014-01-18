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
!                        MODULE idsrdr_iostream                         !
!  *******************************************************************  !
!  Description: controlled opening/closing of a "stream" access file.   !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *******************************************************************  !

MODULE idsrdr_iostream

  implicit none

  PRIVATE ! default is private
  PUBLIC :: openstream, closestream


CONTAINS


!  *******************************************************************  !
!                              openstream                               !
!  *******************************************************************  !
!  Description: open the file with system name 'file' and logical unit  !
!  number 'lunit' for stream access (unformated) where 'LRECL' is the   !
!  size of the items to be read/write in the file. In case of failure,  !
!  the optional parameter 'IOSTAT' returns the Fortran code for the     !
!  first error detected and an error message is always displayed. If    !
!  the file is connected already to another logical unit number         !
!  different from 'lunit' or if 'lunit' is already associated with      !
!  another file, the connection is closed to avoid conflict. If the     !
!  file doesn't exist, it is created.                                   !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  ****************************** INPUT ******************************  !
!  character(*) file       : File system name                           !
!  integer lunit           : Data file logical unit number              !
!  ***************************** OUTPUT ******************************  !
!  integer IOSTAT          : [optional] Returns 0 if succeeds, else     !
!                            returns a Fortran error code               !
!  *******************************************************************  !
  subroutine openstream (file, lunit, IOSTAT)

!   Input variables.
    character(*), intent(in) :: file
    integer, intent(in) :: lunit
    integer, intent (out), optional :: IOSTAT

!   Local variables.
    logical :: F_open, U_open, F_exist
    character(11) :: acc, act, frm
    integer :: LUN, Koderr, i

    LUN = ABS(lunit)
    if (present(IOSTAT)) IOSTAT = 0
    inquire (LUN, OPENED=U_open)
    inquire (FILE=file, OPENED=F_open, NUMBER=i, EXIST=F_exist)
    IF (U_open) THEN
       If (F_open) Then
          if (LUN == i) then
             inquire (FILE=file, ACCESS=acc, ACTION=act, FORM=frm)
             if (acc(:1)//act(:1)//frm(:1) == "SRU") return
          endif
          close (i)
       EndIf
       close (LUN)
    ELSE
       If (F_open) close (i)
    ENDIF

    open (LUN, FILE=file, STATUS=merge("OLD","NEW",F_exist),            &
         ACCESS="STREAM", ACTION="READWRITE", POSITION='REWIND',        &
         IOSTAT=Koderr)

    if (Koderr == 0) return
    if (present(IOSTAT)) IOSTAT = Koderr

    print*, "ERROR loading the file ", file
    print*, "Fortran error code = ", Koderr


  end subroutine openstream


!  *******************************************************************  !
!                              closestream                              !
!  *******************************************************************  !
!  Description: close the file with system name 'file' and logical      !
!  unit number 'lunit'. In case of failure, the optional parameter      !
!  'IOSTAT' returns the Fortran code for the first error detected and   !
!  an error message is always displayed. If the file is connected to a  !
!  to another logical unit number different from 'lunit' then a error   !
!  message is displayed and the files lun is closed.                    !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  ****************************** INPUT ******************************  !
!  character(*) file       : File system name                           !
!  integer lunit           : Data file logical unit number              !
!  ***************************** OUTPUT ******************************  !
!  integer IOSTAT          : [optional] Returns 0 if succeeds, else     !
!                            returns a Fortran error code               !
!  *******************************************************************  !
  subroutine closestream (file, lunit, IOSTAT)

!   Input variables.
    character(*), intent(in) :: file
    integer, intent(in) :: lunit
    integer, intent (out), optional :: IOSTAT

!   Local variables.
    logical :: F_open, U_open, F_exist
    integer :: LUN, Koderr, i

    LUN = ABS(lunit)
    if (present(IOSTAT)) IOSTAT = 0
    inquire (LUN, OPENED=U_open)
    inquire (FILE=file, OPENED=F_open, NUMBER=i, EXIST=F_exist)
    IF (U_open) THEN
       If (F_open) Then
          if (LUN /= i) then
             print*, "ERROR: at closestream: the logical unit ", LUN
             print*, "       does not macth the file ", file
          endif
          close (i, IOSTAT=Koderr)
          if (Koderr /= 0) then
             print*, "ERROR closing the file ", file
             print*, "Fortran error code = ", Koderr
             if (present(IOSTAT)) IOSTAT = Koderr
             return
          endif
       EndIf
       close (LUN, IOSTAT=Koderr)
       if (Koderr /= 0) then
          print*, "ERROR closing the logical unit ", LUN
          print*, "Fortran error code = ", Koderr
          if (present(IOSTAT)) IOSTAT = Koderr
       endif
    ELSE
       If (F_open) Then
          close (i, IOSTAT=Koderr)
          if (Koderr /= 0) then
             print*, "ERROR closing the file ", file
             print*, "Fortran error code = ", Koderr
             if (present(IOSTAT)) IOSTAT = Koderr
          endif
       EndIf
    ENDIF


  end subroutine closestream


!  *******************************************************************  !


END MODULE idsrdr_iostream
