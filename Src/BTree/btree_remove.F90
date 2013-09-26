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
!                          MODULE btree_remove                          !
!  *******************************************************************  !
!  Description: removal of data in a key-indexed B-Tree.                !
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

MODULE btree_remove

!
! Modules
!
  use btree_main

  implicit none

  PRIVATE ! default is private
  PUBLIC :: BTREEremove, keys

! For comunication between 'BTREEremove' and 'rec_remove'.
  TYPE(keys) :: key
  logical :: rise, change_son, rewrite, nokey

! Working variables for 'rec_remove' (to do not stack).
  logical :: found
  integer :: size, IO_Err

CONTAINS

!  *******************************************************************  !
!                              BTREEremove                              !
!  *******************************************************************  !
!  Description: insert an item from the B-Tree if its there.            !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  TYPE(keys) key_old    : Key to be removed                            !
!  ***************************** OUTPUT ******************************  !
!  integer IOSTAT        : Returns 0 if there is no error               !
!  *******************************************************************  !
  subroutine BTREEremove (key_old, IOSTAT)

!   Input variables.
    TYPE(keys), intent(in) :: key_old
    integer, intent (out), optional :: IOSTAT

!   Local variables.
    integer :: addr_newRoot

!   Initialize variables.
    key = key_old
    if (present(IOSTAT)) IOSTAT = 0
    rewrite = .true. ! current page has been modified
    nokey = .false. ! 'key_old' exists

!   Recursive removal.
    call rec_remove (ROOT, IOSTAT)

    if (nokey) then
       call MSG_NOKEY (KEY_MISSING, IOSTAT)
       return
    endif

    if (rise .and. ROOT%degree < 1) then
       addr_newRoot = ROOT%son(0)
       if (addr_newRoot /+ NULL) & ! switch the root
            read (B_INDEX, REC=addr_newRoot, ERR=99, IOSTAT=IO_Err) ROOT
       call FREE_DCB (INDEX_DCB, addr_newRoot)
    endif

    if (rewrite) write (B_INDEX, REC=LOC_ROOT) ROOT
    return

99  call ACCESS_EXCEPTION (B_INDEX_INIT, IO_Err,                        &
         "New root unavailable", addr_newRoot, IOSTAT)


  end subroutine BTREEremove


!  *******************************************************************  !
!                              rec_remove                               !
!  *******************************************************************  !
!  Description: preform a recursive removal in the B-Tree.              !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ************************** INPUT/OUTPUT ***************************  !
!  TYPE(Bnode) roo       : Node to be removed                           !
!  integer IOSTAT        : Returns 0 if succeeds, else returns a        !
!                          Fortran error code                           !
!  *******************************************************************  !
    recursive subroutine rec_remove (roo, IOSTAT)

!     Input variables.
      TYPE(Bnode), intent(inout) :: roo
      integer, intent (inout), optional :: IOSTAT

!     Local variables.
      integer :: rank, succ
      TYPE(Bnode) :: s_roo

!     Calculates the position.
      call CALC_POS (Less_or_Equal_than=key,                            &
                     Into=roo%content(:roo%degree),                     &
                     Rank=rank, Found=found)

      IF (found) THEN ! remove 'content(rank)'
         call remove (roo%content(rank)%addr)
         succ = roo%son(rank-1) ! last son
         if (succ == NULL) then
            size = roo%degree - 1
            roo%degree = size
            rise = size < ORDER
            roo%content(rank:size) = roo%content(rank+1:size+1)
            return
         endif

         read (B_INDEX, REC=succ, ERR=99, IOSTAT=IO_Err) s_roo
         call largest_leftSubtree (s_roo)
         change_son = rewrite
         rewrite = .true. ! ROOT modified by 'largest_leftSubtree'
         if (rise) call reorder (ROOT, s_roo, succ, rank-1)

      ELSE
         succ = roo%son(rank)
         if (succ /+ NULL) then
            read (B_INDEX, REC=succ, ERR=99, IOSTAT=IO_Err) s_roo
            call rec_remove (s_roo, IOSTAT)
            if (nokey) return
            change_son = rewrite
            rewrite = rise
            if (rise) call reorder (ROOT, s_roo, succ, rank)
         else
            nokey = .true. ! 'key_old' is not there
            return
         endif

      ENDIF

      if (change_son) write (B_INDEX, REC=succ) s_roo
      return

99    call ERR_ACCESS (succ)


    CONTAINS


!     *************************************************************     !
!                          largest_leftSubtree                          !
!     *************************************************************     !
!     Description: go down to the right from the key 'key_old' to       !
!     be removed (in the last left subtree before these key),           !
!     searching for the largest key < 'key_old'. On return, if          !
!     'rewrite' is .true. then the node 'node' has been modified on     !
!     memory.                                                           !
!     *********************** INPUT/OUTPUT ************************     !
!     TYPE(Bnode) node      : Node to be modified                       !
!     *************************************************************     !
      recursive subroutine largest_leftSubtree (node)

!       Input variables.
        TYPE(Bnode), intent(inout), target :: node

!       Local variables.
        TYPE(Bnode) :: r_son
        integer, pointer :: size
        integer :: succ_right

        size => node%degree
        succ_right = node%son(size)
        
        if (succ_right /= NULL) then

           read (B_INDEX, REC=succ_right, ERR=99, IOSTAT=IO_Err) r_son
           call largest_leftSubtree (r_son)
           change_son = rewrite
           rewrite = rise
           if (rise)  call reorder (node, r_son, succ_right, size)
           if (change_son) write (B_INDEX, REC=succ_right) r_son

        else ! 'node' = most right leaf form left subtree

           roo%content(rank) = node%content(size)
           size = size - 1
           rise = size < ORDER

        endif
        return

99      call ERR_ACCESS (succ_right)


      end subroutine largest_leftSubtree


!     *************************************************************     !
!                          largest_leftSubtree                          !
!     *************************************************************     !
!     Description: if the page 'son' is too small (it has no more       !
!     than 'ORDER-1' keys), then use a neighbor page 'brother' in       !
!     order to balance then if possible, or to merge otherwise.         !
!     *************************** INPUT ***************************     !
!     integer s_addr        : "son" address                             !
!     integer rank          : Size of the page                          !
!     *********************** INPUT/OUTPUT ************************     !
!     TYPE(Bnode) father    : "father" page                             !
!     TYPE(Bnode) son       : "son" page                                !
!     *************************************************************     !
      subroutine reorder (father, son, s_addr, rank)

!       Input variables.
        TYPE(Bnode), intent(inout), target :: father
        TYPE(Bnode), intent(inout) :: son
        integer, intent(in) :: s_addr, rank

!       Local variables.
        TYPE(Bnode), target :: brother
        integer :: b_addr, pos, avail_keys
        integer, pointer :: f_size, b_size

        f_size => father%degree
        b_size => brother%degree
        pos = rank

        IF (rank < f_size) THEN ! take 'brother' page on the right
           pos = pos + 1
           b_addr = father%son(pos)
           read (B_INDEX, REC=b_addr, ERR=99, IOSTAT=IO_Err) brother
           son%content(ORDER) = father%content(pos)
           son%son(ORDER) = brother%son(0)

           If (b_size > ORDER) Then ! balance pages 'brother' and 'son'

              avail_keys = (b_size - ORDER + 1) / 2
              son%content(ORDER+1:ORDER+avail_keys-1) =                 &
                   brother%content(:avail_keys-1)
              son%son(ORDER+1:ORDER+avail_keys-1) =                     &
                   brother%son(1:avail_keys-1)
              father%content(pos) = brother%content(avail_keys)
              father%son(pos) = b_addr
              brother%son(0) = brother%son(avail_keys)
              b_size = b_size - avail_keys
              brother%content(:b_size) =                                &
                   brother%content(1+avail_keys:b_size+avail_keys)
              brother%son(1:b_size) =                                   &
                   brother%son(1+avail_keys:b_size+avail_keys)
              son%degree = ORDER - 1 + avail_keys
              rise = .false.

           Else ! merge pages 'brother' and 'son'

              son%content(ORDER+1:) = brother%content(:ORDER)
              son%son(ORDER+1:) = brother%son(1:ORDER)
              father%content(pos:f_size-1) = father%content(pos+1:f_size)
              father%son(pos:f_size-1) = father%son(pos+1:f_size)
              son%degree = Max_Degree
              f_size = f_size - 1
              rise = f_size < ORDER
              call FREE_DCB (INDEX_DCB, b_addr) ! remove 'brother' page
              return

           EndIf

        ELSE ! take 'brother' page on the left

           b_addr = father%son(pos-1)
           read (B_INDEX, REC=b_addr, ERR=99, IOSTAT=IO_Err) brother

           If (b_size > ORDER) Then ! balance pages 'brother' and 'son'

              avail_keys = (b_size - ORDER + 1) / 2
              son%content(avail_keys+1:ORDER+avail_keys-1) =            &
                   brother%content(:ORDER-1)
              son%son(avail_keys+1:ORDER+avail_keys-1) =                &
                   brother%son(1:ORDER-1)
              son%content(avail_keys) = father%content(pos)
              son%son(avail_keys) = father%son(0)
              b_size = b_size - avail_keys + 1
              son%content(:avail_keys-1) =                              &
                   brother%content(1+b_size:b_size+avail_keys-1)
              son%son(1:avail_keys-1) =                                 &
                   brother%son(1+b_size:b_size+avail_keys-1)
              son%son(0) = brother%son(b_size)
              father%content(pos) = brother%content(b_size)
              father%son(pos) = b_addr
              b_size = b_size - 1
              son%degree = ORDER - 1 + avail_keys
              rise = .false.

           Else ! merge pages 'brother' and 'son'

              b_size = b_size + 1
              brother%content(b_size) = father%content(pos)
              brother%son(b_size) = son%son(0)
              brother%content(1+b_size:b_size+ORDER-1) =                &
                   son%content(:ORDER-1)
              brother%son(1+b_size:b_size+ORDER-1) = son%son(1:ORDER-1)
              b_size = Max_Degree
              f_size = f_size - 1
              change_son = .false.
              call FREE_DCB (INDEX_DCB, s_addr) ! remove 'son' page
              rise = f_size < ORDER

           EndIf

        ENDIF ! IF (rank < f_size)

        write (B_INDEX, REC=b_addr) brother
        return

99      call ERR_ACCESS (b_addr)


      end subroutine reorder


!     *************************************************************     !
!                              ERR_ACCESS                               !
!     *************************************************************     !
!     Description: management of errors when accessing the 'Index'      !
!     file (which then became unusable).                                !
!     *************************** INPUT ***************************     !
!     integer numbNR        : unavailable addres                        !
!     *************************************************************     !
      subroutine ERR_ACCESS (numbNR)

!       Input variables.
        integer, intent(in) :: numbNR

        call ACCESS_EXCEPTION (B_INDEX_INIT, IO_Err,                    &
            "Access to index file", numbNR, IOSTAT)


      end subroutine ERR_ACCESS


!     *************************************************************     !


    end subroutine rec_remove


!  *******************************************************************  !


END MODULE btree_remove
