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
!                         MODULE master_slave                           !
!  *******************************************************************  !
!  Description: treats the nodes distribuition (one node 'Master'       !
!  which distributes the work over the other nodes, the 'Slaves').      !
!                                                                       !
!  Written by Alberto Torres, Jan 2014.                                 !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: alberto.trj@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *******************************************************************  !

MODULE master_slave

#include "master-slave.h"

!
! Modules
!
  use parallel,        only:

  implicit none

  PUBLIC  :: Init_Master, Kill_Master, Master_SetupLoop, Slave_AskWork
  PRIVATE :: Master_Distribute


CONTAINS


!  *******************************************************************  !
!                               Init_Master                             !
!  *******************************************************************  !
!  Description: Awakes the master process, which will distribute        !
!  work among the slaves.                                               !
!                                                                       !
!  Written by Alberto Torres, Jan 2014.                                 !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: alberto.trj@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IamMaster           : True if the process is the master      !
!  integer Node                : Number of the process/node/rank        !
!  *******************************************************************  !
  subroutine Init_Master

!
! Modules
!
    use parallel,        only: IamMaster, Node
      
    include "mpif.h"

!   Local variables.
    integer :: MPIerror, MPIstatus(MPI_STATUS_SIZE)
    integer :: Task, Ntotal
!   logical :: Got_EnergyLoop_data = .false.

!   Only Master can do.
    if (.not. IamMaster) then
       write (*,'(a,i4)')                                               &
            "Init_Master: ERROR: Called from slave process ", Node
       call MPI_Abort (MPI_Comm_World, 666, MPIerror)
    end if

!   Wait until work needs to be distributed to slaves,
!   and proceed accordingly, or until the program ends.
    do while (.true.)

!      What kind of task is required from me?
!      To distribute the interations of a loop, perhaps?
       call MPI_Recv (Task,       & ! Task to be distributed
            1, MPI_INTEGER,       & ! size of one integer
            MPI_ANY_SOURCE,       & ! receive from any process
            TASK_TAG,             & ! "I want work" type of message
            MPI_Comm_World,       & ! Global communicator 
            MPIstatus, MPIerror)

!      Let me see what I am going to do.
       select case (Task)

!        The program is ending. Time to exit.
         case(TASK_EXIT)
            call Kill_Master
            return

!        I will distribute energy loop iterations!
         case (TASK_ENERGYLOOP)
!           How many iterations?
            call MPI_Recv (Ntotal,     & ! Number of points in the loop
                 1, MPI_INTEGER,       & ! size of one integer
                 MPI_ANY_SOURCE,       & ! receive from any process
                 INFO_TAG,             & ! "I want work" type of message
                 MPI_Comm_World,       & ! Global communicator 
                 MPIstatus, MPIerror)

!           Prepare your backs, slaves!
!           Assuming the loop begins at 1, in steps of 1.
            call Master_Distribute (1, Ntotal, 1)

       end select

    end do ! while


  end subroutine Init_Master


!  *******************************************************************  !
!                               Kill_Master                             !
!  *******************************************************************  !
!  Description: Finalize master process.                                !
!                                                                       !
!  Written by Alberto Torres, Jan 2014.                                 !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: alberto.trj@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer IamMaster           : True if the process is the master      !
!  integer Node                : Number of the process/node/rank        !
!  *******************************************************************  !
  subroutine Kill_Master

!
! Modules
!
    use parallel,        only: IamMaster, Node

    include "mpif.h"

!   Local variables.
    integer :: MPIerror

    if (.not. IamMaster) then
       write (*,'(a,i4)')                                               &
            "Kill_Master: ERROR: Called from slave process ", Node
       call MPI_Abort (MPI_Comm_World, 666, MPIerror)
    end if
    call MPI_Barrier (MPI_Comm_World, MPIerror) ! Barrier for Master
                                                ! and Slaves
    write (*,'(a)') "Master: exiting"
    call MPI_Finalize (MPIerror)


  end subroutine Kill_Master


!  *******************************************************************  !
!                            Master_Distribute                          !
!  *******************************************************************  !
!  Description: Distributes work (loop interations) among the slaves.   !
!                                                                       !
!  Written by Alberto Torres, Jan 2014.                                 !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: alberto.trj@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer Nodes               : Total number of processes/nodes/ranks  !
!  integer Node                : Number of the process/node/rank        !
!  ****************************** INPUT ******************************  !
!  integer iBegin, iEnd, iStep : Loop limits and step                   !
!  *******************************************************************  !
  subroutine Master_Distribute (iBegin, iEnd, iStep)

!
! Modules
!
    use parallel,        only: Nodes, Node

    include "mpif.h"

!   Input variables.
    integer, intent(in) :: iBegin, iEnd, iStep

!   Local variables.
    integer :: i, WhoWantsWork
    integer :: MPIerror, MPIstatus(MPI_STATUS_SIZE)

    do i = iBegin, iEnd, iStep
!      Wait for some process willing to work.
       call MPI_Recv (WhoWantsWork,   & ! Who asks?
            1, MPI_INTEGER,           & ! one integer
            MPI_ANY_SOURCE,           & ! receive from any sender
            IWANTWORK_TAG,            & ! "I want work" type of message
            MPI_Comm_World,           & ! global communicator 
            MPIstatus, MPIerror)

!       Give command!
        call MPI_Send (i, 1, MPI_INTEGER,  & ! Data, size and type
             WhoWantsWork,                 & ! destination
             WORK_TAG,                     & ! message tag
             MPI_Comm_World, MPIerror)
        write (*,'(a,i4,a,i2)') 'Master: Sending i= ', i,               &
             ' to node ', WhoWantsWork
    end do

!   I must tell my slaves the work is done.
!   Send a terminate message, EndWork=-1, to slaves.
    do i = 1, Nodes
       call MPI_Recv (WhoWantsWork,    & ! Who asks?
            1, MPI_INTEGER,            & ! one integer
            MPI_ANY_SOURCE,            & ! receive from any sender
            IWANTWORK_TAG,             & ! "I want work" type of message
            MPI_Comm_World,            & ! global communicator 
            MPIstatus, MPIerror)

       call MPI_Send (ENDWORK_MSG, 1, MPI_INTEGER, & ! End Work message
            WhoWantsWork,                          & ! destination
            WORK_TAG,                              & ! message tag
            MPI_Comm_World, MPIerror)
    end do


  end subroutine Master_Distribute


!  *******************************************************************  !
!                            Master_SetupLoop                           !
!  *******************************************************************  !
!  Description: Informs the Master the size of the loop.                !
!                                                                       !
!  Written by Alberto Torres, Jan 2014.                                 !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: alberto.trj@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer Master              : Number of Master process/node/rank     !
!  ****************************** INPUT ******************************  !
!  integer Ntotal              : Number of iterations in the loop       !
!                                (points in the energy grid)            !
!  *******************************************************************  !
  subroutine Master_SetupLoop (Ntotal)

!
! Modules
!
    use parallel,        only: IOnode, Master, MPI_Comm_MyWorld

    include "mpif.h"

!   Input variables.
    integer, intent(in) :: Ntotal   ! Number of iterations of the loop

!   Local variables.
    integer :: MPIerror

!   Only 'IOnode' sends info about the loop.
    if (IOnode) then
       call MPI_Send             &
            (TASK_ENERGYLOOP,    & ! Master, we want to work in a loop...
            1, MPI_INTEGER,      & ! data type
            Master,              & ! destination = Master
            TASK_TAG,            & ! user chosen message tag
            MPI_Comm_World, MPIerror)
       call MPI_Send (Ntotal,          & ! ...with Ntotal iterations
            1, MPI_INTEGER,            & ! data type
            Master,                    & ! destination = Master
            INFO_TAG,                  & ! user chosen message tag
            MPI_Comm_World, MPIerror)
    end if

!   Syncronize before requesting work.
!    call MPI_Barrier(MPI_Comm_MyWorld, MPIerror) ! Precisa? Testar!


  end subroutine Master_SetupLoop


!  *******************************************************************  !
!                            Master_Distribute                          !
!  *******************************************************************  !
!  Description: Distributes work (loop interations) among the slaves.   !
!                                                                       !
!  Written by Alberto Torres, Jan 2014.                                 !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: alberto.trj@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer Node                : Number of the process/node/rank        !
!  integer Master              : Number of Master process/node/rank     !
!  ***************************** OUTPUT ******************************  !
!  integer i                   : Itaration to work with                 !
!  *******************************************************************  !
  subroutine Slave_AskWork (i)

!
! Modules
!
    use parallel,        only: Node, Master

    include "mpif.h"

!   Input variables.
    integer, intent(out) :: i

!   Local variables.
    integer :: MPIerror, MPIstatus(MPI_STATUS_SIZE)

!   I must ask Master for work.
    call MPI_Send (Node,                     & ! It's me, Master!
                   1, MPI_INTEGER,           & ! data type
                   Master,                   & ! destination = Master
                   IWANTWORK_TAG,            & ! user chosen message tag
                   MPI_Comm_World, MPIerror)

!   I listen and obey.
    call MPI_Recv (i, 1, MPI_INTEGER,        & ! data, size and type
                   Master,                   & ! receive from Master
                   WORK_TAG,                 & ! "work" type of message
                   MPI_Comm_World, MPIstatus, MPIerror)


  end subroutine Slave_AskWork


!  *******************************************************************  !


END MODULE master_slave
