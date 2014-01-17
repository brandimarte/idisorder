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
!                       MODULE idsrdr_opendirect                        !
!  *******************************************************************  !
!  Description: controlled opening of a "direct" file.                  !
!                                                                       !
!  Based on P. Lignelet, "Structures de Donnees en Fortran 90/95",      !
!  Masson  (1996).                                                      !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *******************************************************************  !

MODULE idsrdr_opendirect

  implicit none


CONTAINS


!  *******************************************************************  !
!                                OPEN_DA                                !
!  *******************************************************************  !
!  Description: open the file with system name 'file' and logical unit  !
!  number 'unit' for direct access (unformated) where 'LRECL' is the    !
!  size of the items to be read/write in the file. In case of failure,  !
!  the optional parameter 'IOSTAT' returns the Fortran code for the     !
!  first error detected and an error message is always displayed. If    !
!  the file is connected already to another logical unit number         !
!  different from 'unit' or if 'unit' is already associated with        !
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
!  integer unit            : Data file logical unit number              !
!  integer LRECL           : Size of items in the file                  !
!  ***************************** OUTPUT ******************************  !
!  integer IOSTAT          : Returns 0 if succeeds, else returns a      !
!                            Fortran error code                         !
!  *******************************************************************  !
  subroutine OPEN_DA (file, unit, LRECL, IOSTAT)

!   Input variables.
    character(*), intent(in) :: file
    integer, intent(in) :: unit, LRECL
    integer, intent (out), optional :: IOSTAT

!   Local variables.
    logical :: F_open, U_open, F_exist
    character(11) :: acc, act, frm
    integer :: LUN, Koderr, i

    LUN = ABS(unit)
    if (present(IOSTAT)) IOSTAT = 0
    inquire (LUN, OPENED=U_open)
    inquire (FILE=file, OPENED=F_open, NUMBER=i, EXIST=F_exist)
    IF (U_open) THEN
       If (F_open) Then
          if (LUN == i) then
             inquire (FILE=file, ACCESS=acc, ACTION=act, FORM=frm)
             if (acc(:1)//act(:1)//frm(:1) == "DRU") return
          endif
          close (i)
       EndIf
       close (LUN)
    ELSE
       If (F_open) close (i)
    ENDIF

    open (LUN, FILE=file, STATUS=merge("OLD","NEW",F_exist),            &
         ACCESS="DIRECT", ACTION="READWRITE", RECL=LRECL, IOSTAT=Koderr)

    if (Koderr == 0) return
    if (present(IOSTAT)) IOSTAT = Koderr

    print*, "Error loading the file ", file
    print*, "Fortran error code = ", Koderr


  end subroutine OPEN_DA


!  *******************************************************************  !


END MODULE idsrdr_opendirect
