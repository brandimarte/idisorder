!  *******************************************************************  !
!  I-Disorder Fortran Code 2007-2014                                    !
!                                                                       !
!  Written by Alexandre Reily Rocha (reilya@ift.unesp.br),              !
!             Pedro Brandimarte (brandimarte@gmail.com) and             !
!             Alberto Torres (alberto.trj@gmail.com).                   !
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
!                         MODULE idsrdr_arrays                          !
!  *******************************************************************  !
!  Description: allocate and initialize arrays to be used on different  !
!  parts of the program.                                                !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *******************************************************************  !

MODULE idsrdr_arrays

!
! Modules
!
  use idsrdr_green,    only: 
  use idsrdr_spectral, only: 
  use idsrdr_power,    only: 
  use idsrdr_current,  only: 

  implicit none

  PUBLIC  :: initarrays
  PRIVATE ! default is private


CONTAINS


!  *******************************************************************  !
!                              initarrays                               !
!  *******************************************************************  !
!  Description: allocate and initialize arrays to be used on different  !
!  parts of the program.                                                !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *******************************************************************  !
  subroutine initarrays

!
! Modules
!
    use idsrdr_green,    only: greeninit
    use idsrdr_spectral, only: spectralinit
    use idsrdr_power,    only: powerinit
    use idsrdr_current,  only: currentinit


!   Initialize Green's functions structures.
    call greeninit

!   Initialize spectral function and DOS arrays.
    call spectralinit

!   Initialize calculated power and occupation arrays.
    call powerinit

!   Initialize calculated current array.
    call currentinit


  end subroutine initarrays


!  *******************************************************************  !


END MODULE idsrdr_arrays
