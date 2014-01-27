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
!                           PROGRAM IDISORDER                           !
!  *******************************************************************  !
!  Description: electron transport in disordered systems with           !
!  electron-phonon interaction.                                         !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  Modified version:    January 2014 (CUDA implementation by Alberto    !
!                       Torres, Universidade de Sao Paulo,              !
!                       e-mail: alberto.trj@gmail.com)                  !
!  *******************************************************************  !

PROGRAM IDISORDER

!
! Modules
!
#ifdef MASTER_SLAVE
  use parallel,        only: IOnode, Node
#else
  use parallel,        only: IOnode
#endif
  use idsrdr_init,     only: init
  use idsrdr_units,    only: makeunits
#ifdef MASTER_SLAVE
  use idsrdr_engrid,   only: engrid, NTenerg_div, Ei, MyEiRecord
#else
  use idsrdr_engrid,   only: engrid, NTenerg_div, Ei
#endif
  use idsrdr_hilbert,  only: hilbertkernel
  use idsrdr_spectral, only: spectralinit, spectral
  use idsrdr_green,    only: greeninit, greenfunctions
#ifdef MASTER_SLAVE
  use idsrdr_options,  only: nspin, NIVP, VInitial, dV, NTenerg
#else
  use idsrdr_options,  only: nspin, NIVP, VInitial, dV
#endif
  use idsrdr_leads,    only: leadsSelfEn
  use idsrdr_current,  only: currentinit, current
  use idsrdr_out,      only: output
  use idsrdr_end,      only: finalize
#ifdef MASTER_SLAVE
  use master_slave,    only: Master_SetupLoop, Slave_AskWork
#endif

  implicit none

#ifdef MASTER_SLAVE
  include "master-slave.h"
#endif

! Local variables.
  integer :: ienergy, ispin, iv
  real(8) :: Vbias

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

! Initialize calculated current array.
  call currentinit

! Initialize Green's functions structures.
  call greeninit

  if (IOnode) write (6,'(/,28("*"),a,28("*"),/)')                       &
       ' Transport Calculation '

#ifdef MASTER_SLAVE
  call Master_SetupLoop (NTenerg)

  do while (.true.)

     call Slave_AskWork (ienergy)

     write (*,'(a,i2,a,i4)') 'Node ', Node, ' received iE = ', ienergy 

     if (ienergy == ENDWORK_MSG) exit

     MyEiRecord (ienergy) = Node ! Keep track which energies I
                                 ! (the node) calculates
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

        Vbias = VInitial

        do iv = 1,NIVP ! over bias points

!          Calculate dissipated power into phonon system.


!          Calculate the current.
           call current (Ei(ienergy), ienergy, ispin, iv, Vbias)

           Vbias = Vbias + dV

        enddo

     enddo
  enddo

! Write calculated values at output files.
  call output

! Proper ending.
  call finalize


!  *******************************************************************  !


END PROGRAM IDISORDER
