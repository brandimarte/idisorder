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
!                          MODULE btree_binary                          !
!  *******************************************************************  !
!  Description: binary search and B-Tree specific settings.             !
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

MODULE btree_binary

!
! Modules
!
  use btree_generic,   only: keys, values, operator(<)

  implicit none

! An 'Item' is given by a key 'key' and an address 'addr'.
  TYPE Item
     TYPE(keys) :: key
     integer :: addr ! in the data file
  END TYPE Item

CONTAINS

!  *******************************************************************  !
!                                the_key                                !
!  *******************************************************************  !
!  Description: given an 'item' it returns its 'key'.                   !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  TYPE(Item) key_addr     : Address of the item                        !
!  *******************************************************************  !
  TYPE(keys) function the_key (key_addr)

!   Input variables.
    TYPE(Item), intent(in) :: key_addr

    the_key = key_addr%key


  end function the_key


!  *******************************************************************  !
!                             Max_Rank_Key                              !
!  *******************************************************************  !
!  Description: retunrs the 'Rank' of the biggest key in a sorted       !
!  vector 'Into' which is less or equal than 'Less_or_Equal_than'. The  !
!  'Rank' is 0 if the the key is <= than all the keys in vector         !
!  'Into'. Also returns the logical variable 'Found' with '.true.' if   !
!  the key is equal to 'Into(Rank)', otherwise returns '.false.'.       !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  TYPE(keys) Less_or_Equal_than  : Key of comparison                   !
!  TYPE(Item) Into(:)             : B-Tree to be search                 !
!  ***************************** OUTPUT ******************************  !
!  integer Rank                   : Rank of the biggest key in 'Into'   !
!                                   which is <= 'Less_or_Equal_than'    !
!  logical Found                  : the key is = 'Into(Rank)'?
!  *******************************************************************  !
  subroutine Max_Rank_Key (Less_or_Equal_than, Into, Rank, Found)

!   Input variables.
    TYPE(keys), intent(in), target :: Less_or_Equal_than
    TYPE(Item), intent(in) :: Into(:)
    integer, intent(out) :: Rank
    logical, intent(out) :: Found

!   Local variables.
    TYPE(keys), pointer :: key
    integer :: middle, sup

    key => Less_or_Equal_than
    Rank = 0
    sup = size(Into) + 1

!   Invariant propertie (with '0 <= Rank,sup <= n+1'):
!     - If 'Rank < sup' then 'Into(Rank) <= key < Into(sup)'
!       (note thet 'key' in 'Into[Rank,sup[')
!     - If 'sup <= Rank' then 'Into(Rank) <= key < Into(Rank+1)'
    DO WHILE (1+Rank < sup)
       middle = (Rank + sup) / 2 ! '0 <= Rank < middle < sup <= n+1'
       if (key < the_key(Into(middle))) then
          sup = middle ! restored
       else ! 'Into(middle) <= key'
          Rank = middle
       endif
    ENDDO

!   '0 <= (Rank = sup-1) < sup <= n+1'
!   then 'key' in 'Into[Rank,Rank+1['
!   then 'key = Into(Rank)' or 'key' is not in 'Into'
    if (Rank >= 1) then
       Found = .not. (the_key (Into(Rank)) < key)
    else
       Found = .false.
    endif


  end subroutine Max_Rank_Key


!  *******************************************************************  !


END MODULE btree_binary
