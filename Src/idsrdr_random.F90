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
!                         MODULE idsrdr_random                          !
!  *******************************************************************  !
!  Description:                                                         !
!                                                                       !
!  Written by Alexandre Reily Rocha, Feb 2007.                          !
!  Instituto de Fisica Teorica                                          !
!  Universidade Estadual de Sao Paulo                                   !
!  e-mail: reilya@ift.unesp.br                                          !
!  ***************************** HISTORY *****************************  !
!  Original version:    February 2007                                   !
!  *******************************************************************  !

MODULE idsrdr_random

!
!   Modules
!
  use parallel,        only: 

  implicit none
  
  PUBLIC  :: random_d, irandomizedefects, irandomize_index

CONTAINS


!  *******************************************************************  !
!                               random_d                                !
!  *******************************************************************  !
!  Description:                                                         !
!                                                                       !
!  Written by Alexandre Reily Rocha, Feb 2007.                          !
!  Instituto de Fisica Teorica                                          !
!  Universidade Estadual de Sao Paulo                                   !
!  e-mail: reilya@ift.unesp.br                                          !
!  ***************************** HISTORY *****************************  !
!  Original version:    February 2007                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode               : True if it is the I/O node            !
!  ****************************** INPUT ******************************  !
!  integer NDefects             : Number of defects blocks              !
!  real*8 avgdist               : Average defect distance               !
!  ***************************** OUTPUT ******************************  !
!  real*8 dist(NDefects)        : Defects distribution                  !
!  real*8 theta(NDefects)       :                                       !
!  *******************************************************************  !
  subroutine random_d (NDefects, avgdist, dist, theta)

!
!   Modules
!
#ifdef MPI
    use parallel,        only: IOnode, MPI_Comm_MyWorld
#else
    use parallel,        only: IOnode
#endif
    use ifport

#ifdef MPI
    include "mpif.h"
#endif

!   Input variables.
    integer, intent(in) :: NDefects
    real(8), intent(in) :: avgdist
    real(8), dimension (NDefects), intent(out) :: dist, theta

!   Local variables.
    integer :: I
    integer, dimension (3) :: now
    real(8) :: Pi, aux
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

    Pi = 4.0d0 * ATAN(1.0d0)

    if (IOnode) then

       call itime (now)
       aux = drand (3600*now(1)+60*now(2)+now(3))

       DO I = 1,NDefects
          theta(i) = 2.0d0 * Pi * drand (0)
       ENDDO

       call itime (now)
       aux = drand (3600*now(1)+60*now(2)+now(3))

       DO I = 1,NDefects       
          dist(i) = 2.d0 * avgdist * drand (0)
       ENDDO

    endif ! if (IOnode)

#ifdef MPI
    call MPI_Bcast (theta, NDefects, MPI_Double_Precision, 0,           &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (dist, NDefects, MPI_Double_Precision, 0,            &
                    MPI_Comm_MyWorld, MPIerror)
#endif


  end subroutine random_d


!  *******************************************************************  !
!                           irandomizedefects                           !
!  *******************************************************************  !
!  Description:                                                         !
!                                                                       !
!  Written by Alexandre Reily Rocha, Feb 2007.                          !
!  Instituto de Fisica Teorica                                          !
!  Universidade Estadual de Sao Paulo                                   !
!  e-mail: reilya@ift.unesp.br                                          !
!  ***************************** HISTORY *****************************  !
!  Original version:    February 2007                                   !
!  ****************************** INPUT ******************************  !
!  integer idefecttypes              :                                  !
!  real*8 unitweight(idefecttypes)   :                                  !
!  logical start                     :                                  !
!  *******************************************************************  !
  integer function irandomizedefects (idefecttypes, unitweight, start)

!
!   Modules
!
    use ifport

!   Input variables.
    integer, intent(in) :: idefecttypes
    real(8), dimension (idefecttypes), intent(in) :: unitweight
    logical, intent(in) :: start

!   Local variables.
    integer :: ii
    integer, dimension (3) :: now
    real(8) :: aux, interval

    if (start) then

       call itime (now)
       aux = drand (3600*now(1)+60*now(2)+now(3))

    else

       aux = drand (0)
       print*, "irandomizedefects: aux 1", aux

    endif

    if (idefecttypes == 0) then 
       irandomizedefects = 1
    else
       ii = 1
       interval = unitweight(ii)
       do while ((aux > interval) .and. (ii <= idefecttypes))
          interval = interval + unitweight(ii+1)
          ii = ii + 1
       enddo
       irandomizedefects = ii
    endif


  end function irandomizedefects


!  *******************************************************************  !
!                           irandomize_index                            !
!  *******************************************************************  !
!  Description:                                                         !
!                                                                       !
!  Written by Alexandre Reily Rocha, Feb 2007.                          !
!  Instituto de Fisica Teorica                                          !
!  Universidade Estadual de Sao Paulo                                   !
!  e-mail: reilya@ift.unesp.br                                          !
!  ***************************** HISTORY *****************************  !
!  Original version:    February 2007                                   !
!  ****************************** INPUT ******************************  !
!  integer nunits                    :                                  !
!  logical start                     :                                  !
!  *******************************************************************  !
  integer function irandomize_index (nunits, start)

!
!   Modules
!
    use ifport

!   Input variables.
    integer, intent(in) :: nunits
    logical, intent(in) :: start

!   Local variables.
    integer, dimension (3) :: now
    real(8) :: aux

    if (start) then

       call itime (now)
       aux = rand (3600*now(1)+60*now(2)+now(3))

    else

       aux = drand (0)
       print*, "irandomize_index: aux 2", aux

    endif

    irandomize_index = NINT(aux*(nunits-1)+2)


  end function irandomize_index


!  *******************************************************************  !


END MODULE idsrdr_random
