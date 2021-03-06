!  *******************************************************************  !
!  I-Disorder Fortran Code 2007-2014                                    !
!                                                                       !
!  Written by Pedro Brandimarte (brandimarte@gmail.com),                !
!             Alberto Torres (alberto.trj@gmail.com) and                !
!             Alexandre Reily Rocha (reilya@ift.unesp.br).              !
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
!                         MODULE idsrdr_engrid                          !
!  *******************************************************************  !
!  Description: create and distribute to nodes the energy grid. If the  !
!  user chooses 'NumberTransmPoints' equal 0, then only Fermi energy    !
!  will be considered. Otherwise, other energy points will be           !
!  considered as defined by 'TEnergI' and 'TEnergF' at input file. In   !
!  this case, the energy points can be interpreted as a gate potential  !
!  that shifts the Fermi energy.                                        !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!                                                                       !
!  Modified by Alberto Torres, Jan 2014.                                !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: alberto.trj@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *******************************************************************  !

MODULE idsrdr_engrid

#ifdef MASTER_SLAVE
#include "master-slave.h"
#endif

!
!   Modules
!
  use parallel,        only: 
  use idsrdr_options,  only: 
  use idsrdr_leads,    only: 

  implicit none
  
  PUBLIC  :: engrid, freegrid, NTenerg_div, Ei
  PRIVATE :: energygrid

  integer :: NTenerg_div ! number of energy grid points per node

#ifdef MASTER_SLAVE
  integer, allocatable, dimension (:) :: MyEiRecord ! Keep track what
                                                    ! node does which Ei
#endif
  real(8), allocatable, dimension (:) :: Ei ! energy grid points


CONTAINS


!  *******************************************************************  !
!                                engrid                                 !
!  *******************************************************************  !
!  Description: subroutine to create and to distribute the energy grid  !
!  to the nodes.                                                        !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer Nodes               : Total number of nodes (MPI_Comm_size)  !
!  integer NTenerg             : Number of transmission energy points   !
!  real*8 TEnergI              : Initial transmission energy            !
!  real*8 TEnergF              : Final transmission energy              !
!  real*8 temp                 : Electronic temperature                 !
!  real*8 EfLead               : Lead Fermi energy                      !
!  ***************************** OUTPUT ******************************  !
!  integer NTenerg_div         : Number of energy grid points per node  !
!  real*8 Ei(NTenerg_div)      : Energy grid points                     !
!  *******************************************************************  !
  subroutine engrid

!
!   Modules
!
    use parallel,        only: IOnode, Nodes
    use idsrdr_options,  only: NTenerg, TEnergI, TEnergF, temp
    use idsrdr_leads,    only: EfLead

    if (IOnode) write (6,'(/,a)') 'engrid: Computing energy grid...'

#if defined MPI && !defined MASTER_SLAVE
    NTenerg_div = NTenerg / Nodes
    if (NTenerg_div == 0) NTenerg_div = 1
#else
    NTenerg_div = NTenerg
#endif

!   Allocate the energy grid points and weights arrays.
    allocate (Ei(NTenerg_div))
#ifdef MASTER_SLAVE
    allocate(MyEiRecord(NTenerg_div))
    MyEiRecord = NOT_ME
#endif

!   Compute the energy grid.
    if (NTenerg == 0) then
       TEnergI = EfLead
       TEnergF = EfLead
       Ei(NTenerg_div) = EfLead
    else
       call energygrid (NTenerg, NTenerg_div, TEnergI, TEnergF, Ei)
    endif
    if (IOnode) write(6,'(/,a)') 'engrid: done!'


  end subroutine engrid


!  *******************************************************************  !
!                              energygrid                               !
!  *******************************************************************  !
!  Description: computes the energy grid.                               !
!                                                                       !
!  Written by Alberto Torres, Jan 2014.                                 !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: alberto.trj@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:  January 2014                                      !
!  *********************** INPUT FROM MODULES ************************  !
!  integer Node                : Actual node (MPI_Comm_rank)            !
!  ****************************** INPUT ******************************  !
!  integer NTenerg             : Number of transmission energy points   !
!  integer NTenerg_div         : Number of energy grid points per node  !
!  real*8 TEnergI              : Initial transmission energy            !
!  real*8 TEnergF              : Final transmission energy              !
!  ***************************** OUTPUT ******************************  !
!  real*8 Ei(NTenerg_div)      : Energy grid points                     !
!  *******************************************************************  !
  subroutine energygrid (NTenerg, NTenerg_div, TEnergI, TEnergF, Ei)

!
!   Modules
!
    use parallel,        only: Node
    use idsrdr_options,  only: deltaEn

!   Input variables.
    integer, intent(in) :: NTenerg, NTenerg_div
    real(8), intent(in) :: TenergI, TenergF
    real(8), dimension (NTenerg_div), intent(out) :: Ei

!   Local variables.
    integer :: i

!   Initializations.  
    Ei = 0.d0
    deltaEn = (TenergF - TenergI) / DBLE(NTenerg)

#if defined MPI && !defined MASTER_SLAVE
!   Only MPI
    do i = 1, NTenerg_div
!      Energ = initial + node shift          + point shift
       Ei(i) = TenergI + Node*NTenerg_div*deltaEn + (i-1)*deltaEn
#else
!   MASTER_SLAVE and serial fall here
    do i = 1, NTenerg
       Ei(i) = TenergI + (i-1)*deltaEn
#endif
    end do


  end subroutine energygrid


!  *******************************************************************  !
!                               freegrid                                !
!  *******************************************************************  !
!  Description: free allocated vectors.                                 !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *******************************************************************  !
  subroutine freegrid


!   Free memory.
    deallocate (Ei)


  end subroutine freegrid


!  *******************************************************************  !


END MODULE idsrdr_engrid

