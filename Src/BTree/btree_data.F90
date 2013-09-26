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
!                           MODULE btree_data                           !
!  *******************************************************************  !
!  Description: management of data on a key-indexed B-Tree.             !
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

MODULE btree_data

!
! Modules
!
  use btree_generic,   only: values
  use alloc_direct
  use open_direct,     only:
  use fortran_except
  use fortran_init

  implicit none

  integer, parameter :: B_ERR_ACCESS = 17 ! error at an access

  integer, save :: B_Fdat ! B-Tree data file logical unit number
  integer, pointer, save :: B_HEAD(:) ! Head of B-Tree
  TYPE(DCB), save :: DATA_DCB
  TYPE(T_INIT), save, PRIVATE :: B_DATA_INIT
  integer, PRIVATE :: Koderr ! working variable
  TYPE(values) :: valinq

CONTAINS

!  *******************************************************************  !
!                           ACCESS_EXCEPTION                            !
!  *******************************************************************  !
!  Description: management of I/O errors.                               !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  TYPE(T_INIT) F_INIT   :                                              !
!  integer IO_Err        : Completion status of the I/O operation       !
!  character* MSG        : Error message                                !
!  integer numberNR      : Cell number of a record                      !
!  ***************************** OUTPUT ******************************  !
!  integer IOSTAT        : Return 0 if there is no error                !
!  *******************************************************************  !
  subroutine ACCESS_EXCEPTION (F_INIT, IO_Err, MSG, numberNR, IOSTAT)

!   Input variables.
    TYPE(T_INIT), intent(in) :: F_INIT
    integer, intent(in) :: IO_Err, numberNR
    character(*), intent(in) :: MSG
    integer, intent (out), optional :: IOSTAT
    
    if (IO_Err /= 0) then ! error
       call EXCEPTION_MESSAGE (F_INIT, IOSTAT, B_ERR_ACCESS, MSG)
       if (numberNR > 0) print*, "Error detected on cell number ",      &
            numberNR
       print*, "Fortran error code = ", IO_Err
    else
       if (present(IOSTAT)) IOSTAT = 0
    endif


  end subroutine ACCESS_EXCEPTION


!  *******************************************************************  !
!                             OPEN_and_INIT                             !
!  *******************************************************************  !
!  Description: open and initialize data file.                          !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  character* file       : DSN (Data Source Name) of the data file      !
!  integer unit          : LUN (Logical Unit Number)                    !
!  ***************************** OUTPUT ******************************  !
!  integer IOSTAT        : Return 0 if succeeds, else returns a         !
!                          Fortran error code                           !
!  *******************************************************************  !
  subroutine OPEN_and_INIT (file, unit, IOSTAT)

!
! Modules
!
    use open_direct,     only: OPEN_DA

!   Input variables.
    character(*), intent(in) :: file
    integer, intent(in) :: unit
    integer, intent (out), optional :: IOSTAT

!   Local variables.
    integer :: IO_Err, LRECL

    call INITIALIZE (B_DATA_INIT)
    inquire (IOLENGTH=LRECL) valinq
    call OPEN_DA (file, unit, MAX(IOlength(0),LRECL), IO_Err)
    B_Fdat = unit
    if (IO_Err == 0) then
       call LOAD_HEADER (unit, 0, DATA_DCB, B_HEAD, IOSTAT)
    else
       if (present(IOSTAT)) IOSTAT = IO_Err
    endif

  end subroutine OPEN_and_INIT


!  *******************************************************************  !
!                              BDATAwrite                               !
!  *******************************************************************  !
!  Description: insert in data file the data 'val' and returns its      !
!  position (its address) in the file.                                  !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  TYPE(values) val      : Data to be inserted                          !
!  ***************************** OUTPUT ******************************  !
!  integer number        : Is the cell number (address) of a record to  !
!                          be accessed directly in data file            !
!  integer IOSTAT        : Return 0 if succeeds, else returns a         !
!                          Fortran error code                           !
!  *******************************************************************  !
  subroutine BDATAwrite (val, number, IOSTAT)

!   Input variables.
    TYPE(values), intent(in) :: val
    integer, intent(out) :: number
    integer, intent (out), optional :: IOSTAT

    call ALLOC_DCB (DATA_DCB, number)
    write (B_Fdat, REC=number, IOSTAT=Koderr) val
    call ACCESS_EXCEPTION (B_DATA_INIT, Koderr,                         &
         "Error on data insertion", number, IOSTAT)


  end subroutine BDATAwrite


!  *******************************************************************  !
!                              BDATAremove                              !
!  *******************************************************************  !
!  Description: remove the item 'number' from data file and add it to   !
!  the liste of free items.                                             !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  integer number        : Item number (address) in data file           !
!  *******************************************************************  !
  subroutine BDATAremove (number)

!   Input variables.
    integer, intent(in) :: number

    call FREE_DCB (DATA_DCB, number)


  end subroutine BDATAremove


!  *******************************************************************  !
!                               BDATAread                               !
!  *******************************************************************  !
!  Description: recovers the value of the item 'number' in data file.   !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  integer number        : Item number (address) in data file           !
!  ***************************** OUTPUT ******************************  !
!  integer IOSTAT        : Return 0 if succeeds, else returns a         !
!                          Fortran error code                           !
!  *******************************************************************  !
  TYPE(values) function BDATAread (number, IOSTAT)

!   Input variables.
    integer, intent(in) :: number
    integer, intent (out), optional :: IOSTAT

    read (B_Fdat, REC=number, IOSTAT=Koderr) BDATAread
    call ACCESS_EXCEPTION (B_DATA_INIT, Koderr,                         &
         "Error on data reading", number, IOSTAT)


  end function BDATAread


!  *******************************************************************  !
!                              BDATAmodify                              !
!  *******************************************************************  !
!  Description: assign a new value 'new_val' to the item 'number' in    !
!  data file.                                                           !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  integer number        : Item number (address) in data file           !
!  TYPE(values) new_val  : New data                                     !
!  ***************************** OUTPUT ******************************  !
!  integer IOSTAT        : Return 0 if succeeds, else returns a         !
!                          Fortran error code                           !
!  *******************************************************************  !
   subroutine BDATAmodify (number, new_val, IOSTAT)

!   Input variables.
    integer, intent(in) :: number
    TYPE(values), intent(in) :: new_val
    integer, intent (out), optional :: IOSTAT

    if (IN_LIMITS(DATA_DCB, number)) then
       write (B_Fdat, REC=number, ERR=99, IOSTAT=Koderr) new_val
       if (present(IOSTAT)) IOSTAT = 0
       return
    endif
    Koderr = HUGE(0) ! modification out of the file bounds
99  call ACCESS_EXCEPTION (B_DATA_INIT, Koderr,                         &
         "Error on data modification", number, IOSTAT)


  end subroutine BDATAmodify


!  *******************************************************************  !


END MODULE btree_data
