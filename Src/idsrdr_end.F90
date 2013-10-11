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

!
! Modules
!
  use parallel,        only: 
  use idsrdr_options,  only: 
  use idsrdr_leads,    only: 
  use idsrdr_engrid,   only: 
  use idsrdr_units,    only: 

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
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode                 : True if it is the I/O node          !
!  integer label_length           : Length of system label              !
!  character(label_length) slabel : System Label (for output files)     !
!  *******************************************************************  !
  subroutine finalize

!
! Modules
!
    use parallel,        only: IOnode
    use idsrdr_options,  only: label_length, slabel, freeopt
    use idsrdr_leads,    only: freeleads
    use idsrdr_engrid,   only: freegrid
    use idsrdr_units,    only: freeunits

    include "mpif.h"

!   Local variables.
    character(len=label_length+4), external :: paste
    external :: timer
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

    call MPI_Barrier (MPI_Comm_world, MPIerror)

!   Free memory.
    if (IOnode) write (6,'(/,30("*"),a,30("*"))')                       &
            ' Ending I-Disorder '

    if (IOnode) write (6,'(/,a)', ADVANCE='no')                         &
         'finalize: Freeing memory...'
    call freegrid
    call freeleads
    call freeopt
    call freeunits

    if (IOnode) then
       write (6,'(a,/)') ' done!'
       write (6,'(a,a)') "Total transmission written to file: ",        &
            paste (slabel,'.TRC') 
       write (6,'(/,a,/)') "End of program I-Disorder"
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
