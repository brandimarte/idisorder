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
!                          MODULE idsrdr_init                           !
!  *******************************************************************  !
!  Description: intialize the program properly and read input options.  !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !

MODULE idsrdr_init

!
! Modules
!
  use parallel,        only:
  use idsrdr_options,  only: 
  use idsrdr_leads,    only: 
  use fdf
  use master_slave,    only: 
  use idsrdr_io,       only: 

  implicit none

  PUBLIC  :: init, nsc, time_begin
  PRIVATE :: header, initread

! Program working variables.
  integer :: nsc(2)

! Initial processor time in seconds (processor-dependent approximation).
  real(8) :: time_begin


CONTAINS


!  *******************************************************************  !
!                                 init                                  !
!  *******************************************************************  !
!  Description: initializes timer and the MPI execution environment.    !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *********************** INPUT FROM MODULES ************************  !
!  integer Node                : Actual node (MPI_Comm_rank)            !
!  integer Nodes               : Total number of nodes (MPI_Comm_size)  !
!  logical IOnode              : True if it is the I/O node             !
!  ***************************** OUTPUT ******************************  !
!  real*8 time_begin              : Initial processor time in seconds   !
!  integer nsc(2)              : Number of unit cells along parallel    !
!                                directions                             !
!  *******************************************************************  !
  subroutine init

!
! Modules
!
#ifdef MPI
#ifdef MASTER_SLAVE
    use parallel,        only: Node, Nodes, IOnode, MPI_Comm_MyWorld,   &
                               Master, IamMaster
#else
    use parallel,        only: Node, Nodes, IOnode, MPI_Comm_MyWorld
#endif
#else
    use parallel,        only: Node, Nodes, IOnode
#endif
    use idsrdr_options,  only: readopt
    use idsrdr_leads,    only: readleads
    use idsrdr_io,       only: IOinit

#ifdef MPI
    include "mpif.h"

!   Local variables.
    integer :: MPIerror ! Return error code in MPI routines
#endif
#ifdef MASTER_SLAVE
    integer :: group_World, group_Slave ! MPI groups
#endif


!   Initialise MPI and set processor number.
#ifdef MPI
    call MPI_Init (MPIerror)
    call MPI_Comm_rank (MPI_Comm_World, Node, MPIerror)
    call MPI_Comm_size (MPI_Comm_World, Nodes, MPIerror)
#endif

    IOnode = (Node == 0)

!   Print version information.
    if (IOnode) then

!      Print header.
       call header (Nodes)

!      Initial time.
       call cpu_time (time_begin)

    endif

!   Init logical units control.
    call IOinit

!## BEGIN Alberto

!   Use 'MPI_Comm_World' only to communicate with
!   the 'Master' otherwise use 'MPI_Comm_MyWorld'.
#ifdef MASTER_SLAVE
    Master = Nodes - 1 ! defining Master as the last process (rank)
    if (Node == Master) IamMaster=.true. ! Slaves get the default .false.

!   Create "world" communicator for slaves: 'MPI_Comm_MyWorld'.
    call MPI_Comm_group (MPI_Comm_World, group_World, MPIerror)
    call MPI_Group_excl (group_World, 1, Master, group_Slave, MPIerror)
    call MPI_Comm_create (MPI_Comm_World, group_Slave,                  &
                          MPI_Comm_MyWorld, MPIerror)

    Nodes = Nodes - 1
    call MPI_Barrier (MPI_Comm_World, MPIerror)
    if (IamMaster) then
       call Init_Master
       call exit(0)
    endif
#elif defined MPI
!   No master, no slave; and we are all equals:
    MPI_Comm_MyWorld = MPI_Comm_World
#endif

!## END Alberto

!   Initialise read.
    call initread

!   Read simulation data.
    call readopt

!   Number of unit cells along parallel directions.
    nsc = 1

!   Read leads input data.
    call readleads (nsc)


  end subroutine init


!  *******************************************************************  !
!                               initread                                !
!  *******************************************************************  !
!  Description: initialize the reading of the data for I-Disorder.      !
!                                                                       !
!  Uses FDF (Flexible Data Fromat) package of J.M.Soler and A.Garcia.   !
!                                                                       !
!  Taken from 'reinit' subroutine from SIESTA package.                  !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  *******************************************************************  !
  subroutine initread

!
! Modules
!
    use parallel,        only: IOnode
    use idsrdr_io,       only: IOassign, IOclose

!   Local variables.
    character string*20
    character filein*20, fileout*20, line*150 
    integer :: count, length, lun, lun_tmp
    logical :: debug_input, file_exists

!   Print welcome and presentation.
    if (IOnode) then

       write (6,'(/,27(1h*),a,27(1h*))') ' Dump of input data file '

!      Dump data file to output file and generate scratch file
!      for FDF to read from (except if INPUT_DEBUG exists).
       inquire (file='INPUT_DEBUG', exist=debug_input)
       if (debug_input) then
          write (6,'(a)') 'WARNING: ' //                                &
               'I-Disorder is reading its input from file INPUT_DEBUG'
           
          call IOassign(lun)
          filein = 'INPUT_DEBUG'
          open (lun, file='INPUT_DEBUG', form='formatted', status='old')
          rewind (lun)
       else
          write (6,'(a,/)') 'initread: Reading from standard input'
          lun = 5
          call IOassign (lun_tmp)
          do
             call system_clock (count)
             write (string,*) count
             filein = 'INPUT_TMP.' // adjustl(string)
             inquire (file=filein, exist=file_exists)
             if (.not. file_exists) exit
          end do
          open (lun_tmp, file=filein, form='formatted', status='replace')
          rewind (lun_tmp)
       endif

10     continue
       read (lun, err=20, end=20, fmt='(a)') line
       call chrlen (line, 0, length)
       if (length /= 0) then
          write(6,'(a,a)') 'initread: ', line(1:length)
          if (.not. debug_input) write (lun_tmp,'(a)') line(1:length)
       endif
       goto 10
20     continue

       write (6,'(2a)') 'initread: ', repeat('*', 69)

!      Choose proper file for fdf processing.
       if (debug_input) then
          call IOclose (lun)
       else
          call IOclose (lun_tmp)
       endif

!      Set up fdf.
       fileout = 'fdf.log'
       call fdf_init (filein, fileout)

    endif


  end subroutine initread


!  *******************************************************************  !
!                                header                                 !
!  *******************************************************************  !
!  Description: prints welcome message, date and copyright.             !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  integer Nodes               : Total number of nodes (MPI_Comm_size)  !
!  *******************************************************************  !
  subroutine header (Nodes)

!   Input variables.
    integer, intent(in) :: Nodes

!   Local variables.
    integer, dimension(8) :: values

    call date_and_time (VALUES=values)

    write (6,'(/,a,73(1h*),/)') '   '
    write (6,'(a,/)')                                                   &
         '                         *  WELCOME TO I-DISORDER  *'
    write (6,'(a,i4,a,i2,a,i2,a,i2,a,i2,a,i2)')                         &
         '                            ',                                &
         values(1), '.', values(2), '.', values(3), ' , ',              &
         values(5), ':', values(6), ':', values(7)
    write (6,'(/,a,a)') '      ', 'Written by Alexandre Reily Rocha' // &
         ' and Pedro Brandimarte, 2007-2013'
    write (6,'(/,a,a)') '      ', 'CUDA implementation by Alberto '  // &
         'Torres, 2014'
    write (6,'(/,a,a)') '      ', 'Copyright (c), All Rights Reserved'
    write (6,'(/,a,a)') '      ', 'This program is free software. '  // &
         'You can redistribute it and/or'
    write (6,'(a,a)') '      ', 'modify it under the terms of the '  // &
         'GNU General Public License'
    write (6,'(a,a)') '      ', '(version 3 or later) as published ' // &
         'by the Free Software Foundation'
    write (6,'(a,a)') '      ', '<http://fsf.org/>. See the GNU '    // &
         'General Public License for details.'
    write (6,'(/,a,73(1h*))') '   '

#ifdef MPI
    if (Nodes > 1) then
       write(6,'(/,a,i4,a)')                                            &
            '* Running on ', Nodes, ' nodes in parallel'
    else
       write(6,'(/,a,i4,a)') '* Running in serial mode with MPI'
    endif
#else
    write(6,'(/,a,i4,a)') '* Running in serial mode'
#endif


  end subroutine header


!  *******************************************************************  !


END MODULE idsrdr_init
