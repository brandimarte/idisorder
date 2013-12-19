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
    use parallel,        only: IOnode
#ifndef IBM
    use ifport
#endif

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
#ifdef IBM
    integer, allocatable, dimension (:) :: seed
    integer, external :: time
#endif
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

    Pi = 4.0d0 * ATAN(1.0d0)

    if (IOnode) then

#ifdef IBM
       now = 0
       now(1) = time()
       call random_seed (SIZE=now(2))
       allocate (seed(now(2)))
       seed = now(1) + now(2)
       call random_seed (PUT=seed)
       call random_number (aux)
       print*, "aux1", aux
       deallocate (seed)
#else
       call itime (now)
       aux = drand (3600*now(1)+60*now(2)+now(3))
#endif

       DO I = 1,NDefects
#ifdef IBM
          call random_number (aux)
          print*, "aux random", aux
          theta(i) = 2.0d0 * Pi * aux
#else
          theta(i) = 2.0d0 * Pi * drand (0)
#endif
       ENDDO

#ifdef IBM
       now(1) = time()
       call random_seed (GENERATOR=2)
       allocate (seed(now(2)))
       seed = now(1) + now(2)
       call random_seed (PUT=seed)
       call random_number (aux)
       print*, "aux", aux
       deallocate (seed)
#else
       call itime (now)
       aux = drand (3600*now(1)+60*now(2)+now(3))
#endif

       DO I = 1,NDefects       
#ifdef IBM
          call random_number (aux)
          dist(i) = 2.d0 * avgdist * aux
#else
          dist(i) = 2.d0 * avgdist * drand (0)
#endif
       ENDDO

    endif ! if (IOnode)

#ifdef MPI
    call MPI_Bcast (theta, NDefects, MPI_Double_Precision, 0,           &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (dist, NDefects, MPI_Double_Precision, 0,            &
                    MPI_Comm_world, MPIerror)
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


#ifndef IBM
    use ifport
#endif

!   Input variables.
    integer, intent(in) :: idefecttypes
    real(8), dimension (idefecttypes), intent(in) :: unitweight
    logical, intent(in) :: start

!   Local variables.
    integer :: ii
    integer, dimension (3) :: now
    real(8) :: aux, interval
#ifdef IBM
    integer, dimension (:), allocatable :: seed
    integer, external :: time
#endif

    if (start) then

#ifdef IBM
       now(1) = time()
       call random_seed (SIZE=now(2))
       allocate (seed(now(2)))
       seed = now(1) + now(2)
       call random_seed (PUT=seed)
       call random_number (aux)
       deallocate (seed)
#else
       call itime (now)
       aux = drand (3600*now(1)+60*now(2)+now(3))
#endif

    else

#ifdef IBM
       call random_number (aux)
#else
       aux = drand (0)
#endif
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


#ifndef IBM
    use ifport
#endif

!   Input variables.
    integer, intent(in) :: nunits
    logical, intent(in) :: start

!   Local variables.
    integer, dimension (3) :: now
    real(8) :: aux
#ifdef IBM
    integer, allocatable, dimension (:) :: seed
    integer, external :: time
#endif

    if (start) then

#ifdef IBM
       now(1) = time()
       call random_seed (SIZE=now(2))
       allocate (seed(now(2)))
       seed = now(1) + now(2)
       call random_seed (PUT=seed)
       call random_number (aux)
       deallocate (seed)
#else
       call itime (now)
       aux = rand (3600*now(1)+60*now(2)+now(3))
#endif

    else

#ifdef IBM
       call random_number (aux)
#else
       aux = drand (0)
#endif
       print*, "irandomize_index: aux 2", aux

    endif

    irandomize_index = NINT(aux*(nunits-1)+2)


  end function irandomize_index


!  *******************************************************************  !


END MODULE idsrdr_random
