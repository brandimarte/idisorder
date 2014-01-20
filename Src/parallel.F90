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
!                            MODULE parallel                            !
!  *******************************************************************  !
!  Description: set/store parallel variables.                           !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !

MODULE parallel


  implicit none

  PUBLIC ! default is public

  integer, save :: Node = 0    ! Actual node (rank)

  integer, save :: Nodes = 1   ! Total number of nodes (comm_size)

  logical, save :: IOnode = 0  ! True if it is the I/O node


!## Alberto:
  integer, save :: MPI_Comm_MyWorld = 0   ! MPI communicator

#ifdef MASTER_SLAVE
  integer, save :: Master = 0             ! Rank of the Master process

  logical, save :: IamMaster = .false.    ! Only one to rule them all
#endif


!  *******************************************************************  !


END MODULE parallel
