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
!                           MODULE idsrdr_end                           !
!  *******************************************************************  !
!  Description: closes the program properly                             !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !

MODULE idsrdr_end

  implicit none

  PUBLIC ! default is public


CONTAINS


!  *******************************************************************  !
!                               finalize                                !
!  *******************************************************************  !
!  Description: closes the program properly.                            !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !
  subroutine finalize

!
! Modules
!
    use parallel,        only: IOnode
    use idsrdr_init,     only: slabel, label_length

    include "mpif.h"

!   Local variables.
    character(len=label_length+4), external :: paste
    external :: timer
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

    call MPI_Barrier (MPI_Comm_World, MPIerror)

!   Free memory.


    if (IOnode) then
       write (6,'(/,a)') "End of program I-Disorder"
       write (6,'(/,a,a)') "Total transmission written to file: ",      &
            paste (slabel,'.TRC') 
    endif

!   Stop time counter.
    call timer ('i-disorder', 1)

!   Finalizes MPI.
#ifdef MPI
    call MPI_Finalize (MPIerror)
#endif


  end subroutine finalize


!  *******************************************************************  !


END MODULE idsrdr_end
