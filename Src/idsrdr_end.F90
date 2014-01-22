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
#ifdef MASTER_SLAVE
#include "master-slave.h"
#endif

MODULE idsrdr_end

!
! Modules
!
  use parallel,        only: 
  use idsrdr_init,     only: 
  use idsrdr_options,  only: 
  use idsrdr_leads,    only: 
  use idsrdr_engrid,   only: 
  use idsrdr_units,    only: 
  use idsrdr_ephcoupl, only: 
  use idsrdr_green,    only: 
  use idsrdr_spectral, only: 
  use idsrdr_hilbert,  only: 

  implicit none

  PUBLIC  :: finalize
  PRIVATE ! default is private


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
!  real*8 time_begin              : Initial processor time in seconds   !
!  integer label_length           : Length of system label              !
!  character(label_length) slabel : System Label (for output files)     !
!  *******************************************************************  !
  subroutine finalize

!
! Modules
!
    use parallel,        only: IOnode
#ifdef MPI
    use parallel,        only: MPI_Comm_MyWorld
#ifdef MASTER_SLAVE
    use parallel,        only: Master
#endif
#endif
    use idsrdr_init,     only: time_begin
    use idsrdr_options,  only: label_length, slabel, freeopt
    use idsrdr_leads,    only: freeleads
    use idsrdr_engrid,   only: freegrid
    use idsrdr_units,    only: freeunits
    use idsrdr_ephcoupl, only: EPHfree
    use idsrdr_green,    only: freegreen
    use idsrdr_spectral, only: freespectral
    use idsrdr_hilbert,  only: freehilb

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    real(8) :: time_end
    character(len=label_length+4), external :: paste
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

#ifdef MPI
    call MPI_Barrier (MPI_Comm_MyWorld, MPIerror)
#endif

!   Free memory.
    if (IOnode) write (6,'(/,30("*"),a,30("*"))')                       &
            ' Ending I-Disorder '

    if (IOnode) write (6,'(/,a)', ADVANCE='no')                         &
         'finalize: Freeing memory...'
    call freegrid
    call freeleads
    call freeopt
    call freeunits
    call EPHfree
    call freegreen
    call freespectral
    call freehilb

    if (IOnode) then

!      Final time.
       call cpu_time (time_end)

       write (6,'(a,/)') ' done!'
       write (6,'(a,a)') "Total transmission written to file: ",        &
            paste (slabel,'.TRC') 
       write (6,'(/,a,f12.4,a)') "Time of calculation was ",            &
            time_end - time_begin, " seconds"
       write (6,'(/,a,/)') "End of program I-Disorder"
    endif

#ifdef MASTER_SLAVE
!   Send the EXIT message to the master, causing him to finalize his existence (by his own hands!)
    if (IOnode) &
       call MPI_Send(TASK_EXIT, 1, MPI_INTEGER, Master, TASK_TAG, MPI_Comm_World, MPIerror) 

    call MPI_Barrier (MPI_Comm_World, MPIerror)  ! Barrier for Master & Slaves
#endif

!   Finalizes MPI.
#ifdef MPI
    call MPI_Finalize (MPIerror)
#endif


  end subroutine finalize


!  *******************************************************************  !


END MODULE idsrdr_end
