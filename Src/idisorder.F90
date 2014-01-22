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
!                           PROGRAM IDISORDER                           !
!  *******************************************************************  !
!  Description: electron transport in disordered systems with           !
!  electron-phonon interaction.                                         !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!                                                                       !
!  Modified by Alberto Torres, Jan 2014.                                !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: alberto.trj@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !
#ifdef MASTER_SLAVE
#include "master-slave.h"
#endif

PROGRAM IDISORDER

!
! Modules
!
  use parallel,        only: IOnode
#ifdef MASTER_SLAVE
  use parallel,        only: Node
  use master_slave,    only: Master_SetupLoop, Slave_AskWork
  use idsrdr_engrid,   only: MyEiRecord
#endif
  use idsrdr_init,     only: init
  use idsrdr_units,    only: makeunits
  use idsrdr_engrid,   only: engrid, NTenerg_div, Ei
  use idsrdr_hilbert,  only: hilbertkernel
  use idsrdr_spectral, only: spectralinit, spectral, writespectral
  use idsrdr_green,    only: greeninit, greenfunctions
  use idsrdr_options,  only: nspin
#ifdef MASTER_SLAVE
  use idsrdr_options,  only: NTenerg
#endif
  use idsrdr_leads,    only: leadsSelfEn
  use idsrdr_current,  only: current
  use idsrdr_end,      only: finalize

  implicit none

! Local variables.
  integer :: ienergy, ispin

! Proper initialization and reading of input options.
  call init

! Read and build disorder units.
  call makeunits

! Create the energy grid and distribute over the nodes.
  call engrid

! Compute interpolation kernel function.
  call hilbertkernel

! Initialize spectral function and DOS arrays.
  call spectralinit

! Initialize Green's functions structures.
  call greeninit

  if (IOnode) write (6,'(/,28("*"),a,28("*"),/)')                     &
       ' Transport Calculation '

#ifdef MASTER_SLAVE
  call Master_SetupLoop(NTenerg)
  do while (.true.)
     call Slave_AskWork(ienergy)
     write(*,'(a,i2,a,i4)') 'Node ', Node, ' received iE= ', ienergy 
     if (ienergy == ENDWORK_MSG) exit
     MyEiRecord(ienergy) = Node ! Keep track which energies I (the node) calculate
#else
  do ienergy = 1,NTenerg_div ! over energy grid
#endif
     do ispin = 1,nspin ! over spin components

!       Calculate lead's self-energies.
        call leadsSelfEn (Ei(ienergy), ispin)

!       Calculate required Green's functions.
        call greenfunctions (Ei(ienergy), ispin)

!       Compute spectral function and DOS for units with e-ph.
        call spectral (ienergy, ispin)

!       Calculate the current.
        call current (Ei(ienergy), ispin)

     enddo
  enddo

! Write spectral function and DOS to file.
  call writespectral

! Proper ending.
  call finalize


!  *******************************************************************  !


END PROGRAM IDISORDER
