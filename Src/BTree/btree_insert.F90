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
!                          MODULE btree_insert                          !
!  *******************************************************************  !
!  Description: insertion of data in a key-indexed B-Tree.              !
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

MODULE btree_insert

!
! Modules
!
  use btree_main

  implicit none

  PRIVATE ! default is private
  PUBLIC :: BTREEinsert, keys, values

! For comunication between 'BTREEinsert' and 'rec_insert'.
  TYPE(values) :: val
  TYPE(keys) :: key
  logical :: rise, already_there
  integer :: p_right
  TYPE(KEYVAL) :: item

CONTAINS

!  *******************************************************************  !
!                              BTREEinsert                              !
!  *******************************************************************  !
!  Description: insert a new item in the B-Tree.                        !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  TYPE(keys) key_new    : New key to be inserted                       !
!  TYPE(values) val_new  : New value to be inserted                     !
!  ***************************** OUTPUT ******************************  !
!  integer IOSTAT        : Return 0 if there is no error                !
!  *******************************************************************  !
  subroutine BTREEinsert (key_new, val_new, IOSTAT)

!   Input variables.
    TYPE(keys), intent(in) :: key_new
    TYPE(values), intent(in) :: val_new
    integer, intent (out), optional :: IOSTAT

!   Local variables.
    integer :: addr_oldRoot

!   Initialize variables.
    key = key_new
    val = val_new
    rise = .true.
    already_there = .false.
    item%key = key
    p_right = NULL
    if (present(IOSTAT)) IOSTAT = 0

!   Recursive insertion.
    call rec_insert (ROOT, LOC_ROOT, IOSTAT)

    if (already_there) then
       call EXCEPTION_MESSAGE (B_INDEX_INIT, IOSTAT, KEY_PRESENT,       &
            "KEY PRESENT")
       return
    endif

    if (rise) then ! create a new root
       call alloc (INDEX_DCB, addr_oldRoot)
       write (B_INDEX, REC=addr_oldRoot) ROOT
       ROOT%degree = 1
       ROOT%content(1) = item
       ROOT%son(:1) = (/addr_oldRoot, p_right/)
       write (B_INDEX, REC=LOC_ROOT) ROOT
    endif


  end subroutine BTREEinsert


!  *******************************************************************  !
!                              rec_insert                               !
!  *******************************************************************  !
!  Description: preform a recursive insertion in the B-Tree.            !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  integer addr          : Address of the new node                      !
!  ************************** INPUT/OUTPUT ***************************  !
!  TYPE(Bnode) roo       : Node to be inserted                          !
!  integer IOSTAT        : Return 0 if succeeds, else returns a         !
!                          Fortran error code                           !
!  *******************************************************************  !
    recursive subroutine rec_insert (roo, addr, IOSTAT)

!     Input variables.
      TYPE(Bnode), intent(inout) :: roo
      integer, intent(in) :: addr
      integer, intent (inout), optional :: IOSTAT

!     Local variables.
      integer :: rank, succ
      integer, save :: IO_Err
      TYPE(Bnode) :: new_roo

!     Calculates the position.
      call CALC_POS (Less_or_Equal_than=key,                            &
                     Into=roo%content(:roo%degree),                     &
                     Rank=rank, Found=already_there)

      if (already_there) return

      succ = roo%son(rank)
      if (succ == NULL) then
         call BDATAwrite (val, item%addr)
         call expand_root
      else
         read (B_INDEX, REC=succ, ERR=99, IOSTAT=IO_Err) new_roo
         call rec_insert (new_roo, succ, IOSTAT)
         if (already_there) return
         if (rise) & ! insert 'item' in the right of 'cont(rank)'
              call expand_root
      endif
      return

99    call ACCESS_EXCEPTION (B_INDEX_INIT, IO_Err,                      &
           "Reading of index file", succ, IOSTAT)


    CONTAINS


!     *************************************************************     !
!                              expand_root                              !
!     *************************************************************     !
!     Description: insert a node in the right of 'cont(rank)'.          !
!     *************************************************************     !
      subroutine expand_root

!       Local variables.
        integer :: size, page_addr, copyaddr
        TYPE(Bnode) :: page
        TYPE(KEYVAL) :: copy

        size = roo%degree
        
        IF (size < Max_Degree) THEN

           size = size +1
           roo%degree = size
           rise = .false.
           roo%content(rank+2:size) = roo%content(rank+1:size-1)
           roo%son(rank+2:size) = roo%son(rank+1:size-1)
           roo%content(rank+1) = item
           roo%son(rank+1) = p_right

        ELSE ! the page is filled

           call alloc (INDEX_DCB, page_addr) ! empty page created
           
!          Backup.
           copy = item
           copyaddr = p_right

           If (rank <= ORDER) Then ! insert at left page

              if (rank < ORDER) then
                 item = roo%content(ORDER)
                 p_right = roo%son(ORDER)
                 roo%content(rank+2:ORDER) = roo%content(rank+1:ORDER-1)
                 roo%son(rank+2:ORDER) = roo%son(rank+1:ORDER-1)
                 roo%content(rank+1) = copy
                 roo%son(rank+1) = copyaddr
              endif
              page%content(:ORDER) = roo%content(ORDER+1:)
              page%son(1:ORDER) = roo%son(ORDER+1:)

           Else ! insert at right page

              rank = rank - ORDER
              item = roo%content(ORDER+1)
              p_right = roo%son(ORDER+1)
              page%content(:rank-1) = roo%content(ORDER+2:rank+ORDER)
              page%son(1:rank-1) = roo%son(ORDER+2:rank+ORDER)
              page%content(rank) = copy
              page%son(rank) = copyaddr
              page%content(rank+1:ORDER) = roo%content(ORDER+rank+1:)
              page%son(rank+1:ORDER) = roo%son(ORDER+rank+1:)

           EndIf

           roo%degree = ORDER
           page%degree = ORDER
           write (B_INDEX, REC=page_addr) page ! empty page

        ENDIF ! IF (size < Max_Degree)


      end subroutine expand_root


!     *************************************************************     !


    end subroutine rec_insert


!  *******************************************************************  !


END MODULE btree_insert
