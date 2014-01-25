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
!  Description: controlled opening/closing of files.                    !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *******************************************************************  !

MODULE idsrdr_io

  implicit none

  integer :: min_lun = 10 ! minimum logical unit number
  integer :: max_lun = 99 ! maximum logical unit number

  logical, allocatable, dimension (:) :: lun_is_free ! lun is free?

  PUBLIC :: IOinit, IOassign, IOclose, IOopenStreamNew, IOopenStream,   &
            IOcloseStream, freeIO
  PRIVATE ! default is private


CONTAINS


!  *******************************************************************  !
!                                IOinit                                 !
!  *******************************************************************  !
!  Description: allocate and initializes global variables.              !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *******************************************************************  !
  subroutine IOinit

!   Local variables.
    integer :: nunits

!   Total number of units.
    nunits = max_lun - min_lun

!   Allocate and initialize 'lun_is_free' array.
    allocate (lun_is_free(nunits))
    lun_is_free = .true.


  end subroutine IOinit


!  *******************************************************************  !
!                               IOassign                                !
!  *******************************************************************  !
!  Description: looks for a free unit and assigns it to lun.            !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  ****************************** INPUT ******************************  !
!  integer lun             : Data file logical unit number              !
!  *******************************************************************  !
  subroutine IOassign (lun)

!
!   Modules
!
#ifdef MPI
    use parallel,        only: IOnode, MPI_Comm_MyWorld
#else
    use parallel,        only: IOnode
#endif

!   Input variables.
    integer, intent (out) :: lun

!   Local variables.
    integer :: iostat
    logical :: used
#ifdef MPI
    integer :: MPIerror
#endif

    do lun = min_lun,max_lun
       if (lun_is_free(lun)) then
          inquire (unit=lun, opened=used, iostat=iostat)
          if (iostat /= 0) used = .true.
          lun_is_free(lun) = .false.
          if (.not. used) return
       endif
    enddo

    if (IOnode) then
       write (6,'(/,a,/)') "ERROR: No luns available in io_assign"
#ifdef MPI
       call MPI_Abort (MPI_Comm_MyWorld, 1, MPIerror)
#endif
       stop
    endif


  end subroutine IOassign


!  *******************************************************************  !
!                                IOclose                                !
!  *******************************************************************  !
!  Description: closes a file and sets 'lun' as free.                   !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  ****************************** INPUT ******************************  !
!  integer lun             : Data file logical unit number              !
!  *******************************************************************  !
  subroutine IOclose (lun)

!   Input variables.
    integer, intent (in) :: lun

    close(lun)
    if (lun >= min_lun .and. lun < max_lun) then
       lun_is_free(lun) = .true.
    endif


  end subroutine IOclose
      

!  *******************************************************************  !
!                            IOopenStreamNew                            !
!  *******************************************************************  !
!  Description: open a new file with system name 'file' and logical     !
!  unit number 'lunit' for stream access (unformated) where 'LRECL' is  !
!  the size of the items to be read/write in the file. In case of       !
!  failure, the optional parameter 'IOSTAT' returns the Fortran code    !
!  for the first error detected and an error message is always          !
!  displayed. If the file is connected already to another logical unit  !
!  number different from 'lunit' or if 'lunit' is already associated    !
!  with another file, the connection is closed to avoid conflict. If    !
!  the file doesn't exist, it is created.                               !
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
  subroutine IOopenStreamNew (file, lunit, IOSTAT)

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
             if (acc(:1)//act(:1)//frm(:1) == 'SRU') return
          endif
          close (i)
       EndIf
       close (LUN)
    ELSE
       If (F_open) close (i)
    ENDIF

    open (LUN, FILE=file, STATUS='REPLACE', ACCESS='STREAM',            &
         ACTION='READWRITE', IOSTAT=Koderr)

    if (Koderr == 0) return
    if (present(IOSTAT)) IOSTAT = Koderr

    print*, "ERROR loading the file ", file
    print*, "Fortran error code = ", Koderr


  end subroutine IOopenStreamNew


!  *******************************************************************  !
!                             IOopenStream                              !
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
  subroutine IOopenStream (file, lunit, IOSTAT)

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
             if (acc(:1)//act(:1)//frm(:1) == 'SRU') return
          endif
          close (i)
       EndIf
       close (LUN)
    ELSE
       If (F_open) close (i)
    ENDIF

    open (LUN, FILE=file, STATUS='OLD', ACCESS='STREAM',                &
         ACTION='READWRITE', POSITION='REWIND', IOSTAT=Koderr)

    if (Koderr == 0) return
    if (present(IOSTAT)) IOSTAT = Koderr

    print*, "ERROR loading the file ", file
    print*, "Fortran error code = ", Koderr


  end subroutine IOopenStream


!  *******************************************************************  !
!                             IOcloseStream                             !
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
  subroutine IOcloseStream (file, lunit, IOSTAT)

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


  end subroutine IOcloseStream


!  *******************************************************************  !
!                                freeIO                                 !
!  *******************************************************************  !
!  Description: free allocated vectors.                                 !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *******************************************************************  !
  subroutine freeIO


!   Free memory.
    deallocate (lun_is_free)


  end subroutine freeIO


!  *******************************************************************  !


END MODULE idsrdr_io
