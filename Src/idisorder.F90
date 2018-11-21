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

#ifdef MASTER_SLAVE
#include "master-slave.h"
#endif

!
! Modules
!
#ifdef MASTER_SLAVE
  use master_slave,    only: Master_SetupLoop, Slave_AskWork
  use parallel,        only: IOnode, Node
  use idsrdr_options,  only: nspin, NIVP, VInitial, dV, NTenerg
  use idsrdr_engrid,   only: engrid, NTenerg_div, Ei, MyEiRecord
#else
  use parallel,        only: IOnode
  use idsrdr_options,  only: nspin, NIVP, VInitial, dV
  use idsrdr_engrid,   only: engrid, NTenerg_div, Ei
#endif
  use idsrdr_init,     only: init
  use idsrdr_units,    only: makeunits
  use idsrdr_hilbert,  only: hilbertkernel
  use idsrdr_arrays,   only: initarrays
  use idsrdr_leads,    only: leadsSelfEn
  use idsrdr_green,    only: greenfunctions
  use idsrdr_spectral, only: spectral
  use idsrdr_power,    only: powerTr, power
  use idsrdr_current,  only: currentTr, current
  use idsrdr_conduct,  only: conduct
  use idsrdr_out,      only: output
  use idsrdr_end,      only: finalize

  implicit none

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

! Allocate and initialize some arrays.
  call initarrays


  if (IOnode) write (6,'(/,28("*"),a,28("*"),/)')                       &
       ' Transport Calculation '


#ifdef MASTER_SLAVE
  call Master_SetupLoop (NTenerg)

  do while (.true.)

     call Slave_AskWork (ienergy)

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

        if (IOnode) then
           write (6,'(a)', advance='no')                                &
                '      computing dissipated power and current... '
           flush (6)
        endif

!       Calculate bias independent part of dissipated power.
        call powerTr (ispin)

!       Calculate bias independent part of current.
#ifdef DEBUG
        call currentTr (Ei(ienergy), ienergy, ispin)
#else
        call currentTr (ienergy, ispin)
#endif

        Vbias = VInitial

        do iv = 1,NIVP ! over bias points

!          Calculate dissipated power into phonon system.
           call power (ienergy, ispin, iv, Vbias)

!          Calculate the current.
           call current (Ei(ienergy), ienergy, ispin, iv, Vbias)

           Vbias = Vbias + dV

        enddo

        if (IOnode) write(6,'(a)') " ok!"

     enddo
  enddo

! Compute differential conductances ('dI/dV' and 'd2I/dV2').
  call conduct

! Write calculated values at output files.
  call output

! Proper ending.
  call finalize


!  *******************************************************************  !


END PROGRAM IDISORDER
