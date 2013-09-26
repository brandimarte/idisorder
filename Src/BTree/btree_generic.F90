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
!                         MODULE btree_generic                          !
!  *******************************************************************  !
!  Description: abstract symble-table managed by a key-indexed B-Tree.  !
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

MODULE btree_generic


  implicit none

  integer, parameter :: ORDER = 2 ! B-Tree order

  integer, parameter :: n_long = 20

! The 'keys' are given by the matrix index '(i,j)'.
  TYPE keys
     integer :: ij
  END TYPE keys

! The value 'values' of a give key is double complex number.
  TYPE values
     double complex :: na
  END TYPE values

! Order relation between the keys.  
  INTERFACE operator (<)
     MODULE procedure LT
  END INTERFACE operator (<)

CONTAINS


!  *******************************************************************  !
!                              BTREEremove                              !
!  *******************************************************************  !
!  Description: Order relation between the keys.                        !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  TYPE(keys) k1       : First key to be compared                       !
!  TYPE(keys) k2       : Seconde key to be compared                     !
!  *******************************************************************  !
  logical function LT (k1, k2)

!   Input variables.
    TYPE(keys), intent(in) :: k1, k2

    LT = k1%ij < k2%ij


  end function LT


!  *******************************************************************  !


END MODULE btree_generic
