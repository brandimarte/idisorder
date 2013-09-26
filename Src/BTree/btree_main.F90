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
!                           MODULE btree_main                           !
!  *******************************************************************  !
!  Description: main module of the abstract key-indexed B-Tree. It      !
!  defines operations like: initialize, return the item count, add a    !
!  new item, find an item with a given key, delete an item with a given !
!  key, visit the items in order of their keys etc.                     !
!  The file management is performed by a specific module 'btree_data'.  !
!  The search of a key in a page operates with a particularization of   !
!  the binary search procedure 'MAX_RANK_KEY' (locally renamed          !
!  'CALCUL_POS') imported from the module 'binary'.                     !
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

MODULE btree_main

!
! Modules
!
  use btree_generic,   only: ORDER, keys, values, operator(<)
  use btree_data
  use btree_binary,    only: KEYVAL=>Item, CALC_POS=>Max_Rank_Key
  use fortran_except,  only: KEY_MISSING=>FAIL_RESEARCH
  use fortran_init
  use alloc_direct

  implicit none

  integer, parameter :: KEY_PRESENT = 14, & ! exceptions
                        NON_LOADED = 15

  integer, parameter :: Max_Degree = 2*ORDER

! Each node 'Bnode' has an array of items ('KEYVAL'), its
! degree 'degree' and an array of links to its children 'son'. 
  TYPE Bnode
     integer :: degree, son(0:Max_Degree)
     TYPE(KEYVAL) :: content(Max_Degree)
  END TYPE Bnode

! The 'root' of B-Tree.
  TYPE(Bnode), save :: ROOT
  integer, parameter :: LOC_ROOT = 2

  logical :: ST_INIT = .false.
  integer, save :: B_INDEX ! LUN of Index file
  TYPE(DCB), save :: INDEX_DCB
  TYPE(T_INIT), save :: B_INDEX_INIT
  integer, PRIVATE :: Koderr, loc

CONTAINS


!  *******************************************************************  !
!                               MSG_NOKEY                               !
!  *******************************************************************  !
!  Description: write excepticon message if a key is missing.           !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  integer L_IOS         : Logical IOS number                           !
!  ************************** INPUT/OUTPUT ***************************  !
!  integer IOSTAT        : Returns 0 if succeeds, else returns a        !
!                          Fortran error code                           !
!  *******************************************************************  !
  subroutine MSG_NOKEY (L_IOS, IOSTAT)

!   Input variables.
    integer, intent(in) :: L_IOS
    integer, intent (inout), optional :: IOSTAT

    if (L_IOS == KEY_MISSING) then
       call EXCEPTION_MESSAGE (B_INDEX_INIT, IOSTAT, KEY_MISSING,       &
            "Key missing")
    else
       if (present(IOSTAT)) IOSTAT = L_IOS
    endif


  end subroutine MSG_NOKEY


!  *******************************************************************  !
!                               localize                                !
!  *******************************************************************  !
!  Description: returns the address at 'B_Fdat' of a given 'key' or     !
!  returns 'KEY_MISSING' (through the global variable 'Koderr').        !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  TYPE(keys) key        : Key to be localized                          !
!  ***************************** OUTPUT ******************************  !
!  integer localize      : Returns the address of 'key' otherwise       !
!                          returns 'KEY_MISSING'                        !
!  *******************************************************************  !
  integer function localize (key)

!   Input variables.
    TYPE(keys), intent(in) :: key

!   Local variables.
    integer :: rank, succ
    logical :: found
    TYPE(Bnode) :: head

    LOCALIZE = 0
    head = ROOT

    DO

!      Calculates the position.
       call CALC_POS (Less_or_Equal_than=key,                           &
                      Into=head%content(:head%degree),                  &
                      Rank=rank, Found=found)

       if (found) then
          LOCALIZE = head%content(rank)%addr
          Koderr = 0
          return
       endif
       succ = head%son(rank)
       if (succ == NULL) EXIT
       read (B_INDEX, REC=succ, ERR=99, IOSTAT=Koderr) head

    ENDDO
    Koderr = KEY_MISSING
    return

99  call ACCESS_EXCEPTION (B_INDEX_INIT, Koderr,                        &
         "Access to Index file", succ)


  end function localize


!  *******************************************************************  !
!                               BTREEopen                               !
!  *******************************************************************  !
!  Description: open the data file (system name 'DSN_NAME' and logical  !
!  unity number 'LUN_DATA') and the index file ('DSN_INDEX' and         !
!  'LUN_INDEX'). These files can exist already or not (if they don't    !
!  exist, they are created according to its standard). In case of       !
!  failure, the optional parameter 'IOSTAT' returns the Fortran code    !
!  for the first error detected and an error message is always          !
!  displayed.                                                           !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  character(*) DSN_DATA   : Data file system name                      !
!  character(*) DSN_INDEX  : Index file system name                     !
!  integer LUN_DATA        : Data file logical unit number              !
!  integer LUN_INDEX       : Index file logical unit number             !
!  ***************************** OUTPUT ******************************  !
!  integer IOSTAT          : Returns 0 if succeeds, else returns a      !
!                            Fortran error code                         !
!  *******************************************************************  !
  subroutine BTREEopen (DSN_DATA, LUN_DATA, DSN_INDEX, LUN_INDEX, IOSTAT)

!
! Modules
!
    use open_direct

!   Input variables.
    character(*), intent(in) :: DSN_DATA, DSN_INDEX
    integer, intent (in), optional :: LUN_DATA, LUN_INDEX
    integer, intent (out), optional :: IOSTAT

!   Local variables.
    integer :: ios, addr_root, root_err, LRECL, lundat
    logical :: opened

    ST_INIT = .false.
    ios = 0
    lundat = ABS(LUN_DATA)
    call OPEN_and_INIT (DSN_DATA, LUN_DATA, Koderr)

    call INITIALIZE (B_INDEX_INIT)
    B_INDEX = ABS(LUN_INDEX)
    inquire (IOLENGTH=LRECL) ROOT
    call OPEN_DA (DSN_INDEX, B_INDEX, MAX(IOLENGTH(0),LRECL), ios)
    IF (ios == 0) THEN
       call LOAD_HEADER (B_INDEX, 0, INDEX_DCB, B_HEAD, ios)
       If (ios == 0) Then
          read (B_INDEX, REC=LOC_ROOT, IOSTAT=root_err) ROOT
          if (root_err /= 0) then
             ROOT%degree = 0
             ROOT%son = NULL
             call ALLOC_DCB (INDEX_DCB, addr_root)
             write (B_INDEX, REC=LOC_ROOT, ERR=99, IOSTAT=ios) ROOT
          endif
          inquire (lundat, OPENED=opened)
          opened = opened .and. lundat/=B_INDEX
          if (.not. opened) then
             print*, "Impossible combined opening of files:"
             print*, "Data:", LUN_DATA, DSN_DATA
             print*, "Index", LUN_INDEX, DSN_INDEX
          endif
          ST_INIT = Koderr==0 .and. opened
          if (.not. ST_INIT) then
             close (B_INDEX)
             ios = NON_LOADED
          endif
       EndIf
    ENDIF

99  if (Koderr /= 0) ios = Koderr
    if (present(IOSTAT)) IOSTAT = ios


  end subroutine BTREEopen


!  *******************************************************************  !
!                               KEYvalue                                !
!  *******************************************************************  !
!  Description: recovers the value 'val_key' of a given key 'the_key'   !
!  in the B-Tree. If the key is not present, it evokes 'KEY_MISSING'    !
!  (and 'val_key' remains undeterminate).                               !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  TYPE(keys) the_key      : Key to be recovered                        !
!  ***************************** OUTPUT ******************************  !
!  TYPE(values) val_key    : Value recovered from 'the_key'             !
!  integer IOSTAT          : Returns 0 if succeeds, else returns a      !
!                            Fortran error code                         !
!  *******************************************************************  !
  subroutine KEYvalue (the_key, val_key, IOSTAT)

!   Input variables.
    TYPE(keys), intent(in) :: the_key
    TYPE(values), intent(out) :: val_key
    integer, intent (out), optional :: IOSTAT

    loc = localize (the_key)
    call MSG_NOKEY (Koderr, IOSTAT)
    if (loc > NULL) val_key = BDATAread (loc, IOSTAT)


  end subroutine KEYvalue


!  *******************************************************************  !
!                              BTREEmodify                              !
!  *******************************************************************  !
!  Description: assign a new value 'new_val' to the key 'key'. If the   !
!  key is not present, it evokes 'KEY_MISSING'.                         !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  TYPE(keys) key          : Key to be modified                         !
!  TYPE(values) new_val    : New value to be assigned to the key        !
!  ***************************** OUTPUT ******************************  !
!  integer IOSTAT          : Returns 0 if succeeds, else returns a      !
!                            Fortran error code                         !
!  *******************************************************************  !
  subroutine BTREEmodify (key, new_val, IOSTAT)

!   Input variables.
    TYPE(keys), intent(in) :: key
    TYPE(values), intent(in) :: new_val
    integer, intent (out), optional :: IOSTAT

    loc = localize (key)
    call MSG_NOKEY (Koderr, IOSTAT)
    if (loc > NULL) call BDATAmodify (loc, new_val, IOSTAT)


  end subroutine BTREEmodify


!  *******************************************************************  !
!                                BTREEmap                               !
!  *******************************************************************  !
!  Description: apply a common treatment once and only one to each      !
!  item from the B-Tree.                                                !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  ***************************** OUTPUT ******************************  !
!  integer IOSTAT          : Returns 0 if succeeds, else returns a      !
!                            Fortran error code                         !
!  *******************************************************************  !
    subroutine BTREEmap (treat, IOSTAT)

!     Input variables.
      integer, intent (out), optional :: IOSTAT
      interface
         subroutine treat (key, val)

           use btree_generic

           TYPE(keys), intent(in) :: key
           TYPE(values), intent(in) :: val

         end subroutine treat
      end interface

!     Local variables.
      integer :: IO_Err, succ
      logical :: DATA_open, INDEX_open

    
      if (present(IOSTAT)) IOSTAT = 0
      inquire (B_INDEX, OPENED=INDEX_open, IOSTAT=IO_Err)
      if (IO_Err == 0) inquire (B_Fdat, OPENED=DATA_open, IOSTAT=IO_Err)
      if (ST_INIT .and. INDEX_open .and. DATA_open .and. IO_Err==0) then
         call rec_enum (ROOT)
         if (IO_Err /= 0) &
              call ACCESS_EXCEPTION (B_INDEX_INIT, IO_Err,              &
              "Linking error at Index file", succ, IOSTAT)
      else
         call EXCEPTION_MESSAGE (B_INDEX_INIT, IOSTAT, NON_LOADED,      &
              "B-Tree not loaded or closed")
      endif


    CONTAINS


!     *************************************************************     !
!                                rec_enum                               !
!     *************************************************************     !
!     Description: recursive subroutine to run over all keys.           !
!     *************************************************************     !
      recursive subroutine rec_enum (roo)

!       Input variables.
        TYPE(Bnode), intent(in) :: roo

!       Local variables.
        TYPE(Bnode) :: head
        integer :: i

        DO i = 1,roo%degree+1
           succ = roo%son(i-1)
           if (succ /= NULL) then
              read (B_INDEX, REC=succ, ERR=99, IOSTAT=IO_Err) head
              call rec_enum (head)
              if (IO_Err /= 0) exit
           endif
           if (i > roo%degree) exit
           call treat (roo%content(i)%key,                              &
                       BDATAread(roo%content(i)%addr))
        ENDDO


99    end subroutine rec_enum


!     *************************************************************     !


    end subroutine BTREEmap


!  *******************************************************************  !
!                               BTREEclose                              !
!  *******************************************************************  !
!  Description: close the B-Tree.                                       !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !
  subroutine BTREEclose

    call CLOSE_DCB (INDEX_DCB)
    call CLOSE_DCB (DATA_DCB)
    ST_INIT = .false.


  end subroutine BTREEclose


!  *******************************************************************  !


END MODULE btree_main
