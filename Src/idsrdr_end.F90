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
  use idsrdr_init,     only: 
  use idsrdr_leads,    only: 
  use idsrdr_engrid,   only: 
  use idsrdr_units,    only: 
  use idsrdr_ephcoupl, only: 
  use idsrdr_green,    only: 
  use idsrdr_spectral, only: 
  use idsrdr_hilbert,  only: 
  use idsrdr_current,  only: 
  use idsrdr_power,    only: 
  use idsrdr_conduct,  only: 
  use idsrdr_io,       only: 

  implicit none

#ifdef MASTER_SLAVE
  include "master-slave.h"
#endif

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
!  *******************************************************************  !
  subroutine finalize

!
! Modules
!
#ifdef MPI
#ifdef MASTER_SLAVE
    use parallel,        only: IOnode, MPI_Comm_MyWorld, Master
#else
    use parallel,        only: IOnode, MPI_Comm_MyWorld
#endif
#else
    use parallel,        only: IOnode
#endif
    use idsrdr_init,     only: time_begin
    use idsrdr_leads,    only: freeleads
    use idsrdr_engrid,   only: freegrid
    use idsrdr_units,    only: freeunits
    use idsrdr_ephcoupl, only: freeEPH
    use idsrdr_green,    only: freegreen
    use idsrdr_spectral, only: freespectral
    use idsrdr_hilbert,  only: freehilb
    use idsrdr_current,  only: freecurr
    use idsrdr_power,    only: freepower
    use idsrdr_conduct,  only: freedIdV
    use idsrdr_io,       only: freeIO

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    real(8) :: time_end
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
    call freeunits
    call freeEPH
    call freegreen
    call freespectral
    call freehilb
    call freecurr
    call freepower
    call freedIdV
    call freeIO

    if (IOnode) then

!      Final time.
       call cpu_time (time_end)

       write (6,'(a,/)') ' done!'
       write (6,'(a,f12.4,a)') "Time of calculation was ",              &
            time_end - time_begin, " seconds"
       write (6,'(/,a,/)') "End of program I-Disorder"
    endif

#ifdef MASTER_SLAVE
!   Send the EXIT message to the master, causing him
!   to finalize his existence (by his own hands!).
    if (IOnode)                                                         &
       call MPI_Send (TASK_EXIT, 1, MPI_INTEGER, Master, TASK_TAG,      &
                      MPI_Comm_World, MPIerror) 

    call MPI_Barrier (MPI_Comm_World, MPIerror) ! Barrier for 'Master'
                                                ! and 'Slaves'
#endif

!   Finalizes MPI.
#ifdef MPI
    call MPI_Finalize (MPIerror)
#endif


  end subroutine finalize


!  *******************************************************************  !


END MODULE idsrdr_end

