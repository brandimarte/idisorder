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
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !

PROGRAM IDISORDER

!
! Modules
!
  use idsrdr_init,     only: init
  use idsrdr_engrid,   only: engrid, NTenerg_div, Ei
  use idsrdr_units,    only: makeunits
  use idsrdr_options,  only: nspin
  use idsrdr_leads,    only: leadsSelfEn
  use idsrdr_green,    only: greeninit, greenfunctions, gfHead, gfTail
  use idsrdr_end,      only: finalize

  implicit none

! Local variables.
  integer :: ienergy, ispin
  integer, allocatable, dimension (:,:) :: INFO, NCHAN

  call init

  call engrid

  call makeunits

! Allocate arrays.
  allocate (INFO(NTenerg_div,nspin))
  allocate (NCHAN(NTenerg_div,nspin))

  call gfHead
  call greeninit
  do ienergy = 1,NTenerg_div ! over energy grid
     do ispin = 1,nspin ! over spin components

        call leadsSelfEn (Ei(ienergy), ispin, INFO(ienergy,ispin),      &
                          NCHAN(ienergy,ispin))

        call greenfunctions (Ei(ienergy), ispin)

     enddo
  enddo
  call gfTail

! Free memory.
  deallocate (INFO, NCHAN)

  call finalize


!  *******************************************************************  !


END PROGRAM IDISORDER
