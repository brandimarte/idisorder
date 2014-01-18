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
!                         MODULE idsrdr_distrib                         !
!  *******************************************************************  !
!  Description: compute Fermi-Dirac and Bose-Einstein distributions.    !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *******************************************************************  !

MODULE idsrdr_distrib

  implicit none

  PUBLIC :: FermiDirac, BoseEinstein


CONTAINS


!  *******************************************************************  !
!                              FermiDirac                               !
!  *******************************************************************  !
!  Description: calculates the Fermi-Dirac distribution for fermions    !
!  in equilibrium with the electrochemical potential 'mu'.              !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2013                                    !
!  ****************************** INPUT ******************************  !
!  real*8 E                   : Energy point                            !
!  real*8 mu                  : Electrochemical potential               !
!  real*8 kT                  : Electronic temperature (Ry) times       !
!                               Boltzmann constant                      !
!  *******************************************************************  !
  real(8) function FermiDirac (E, mu, kT)

!   Input variables.
    real(8), intent(in) :: mu, kT, E

    if ((E - mu) / kT .gt. 42.d0) then
       FermiDirac = 0.d0
    else
       FermiDirac = 1.d0 / (DEXP((E - mu) / kT) + 1.d0)
    endif


  end function FermiDirac


!  *******************************************************************  !
!                             BoseEinstein                              !
!  *******************************************************************  !
!  Description: calculates the Bose-Einstein distribution for phonons   !
!  with energy 'freq' in thermal equilibrium.                           !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2013                                    !
!  ****************************** INPUT ******************************  !
!  real*8 freq                : Vibrational mode frequency (energy)     !
!  real*8 kT                  : Electronic temperature (Ry) times       !
!                               Boltzmann constant                      !
!  *******************************************************************  !
  real(8) function BoseEinstein (freq, kT)

!   Input variables.
    real(8), intent(in) :: freq, kT

    if (freq/kT .gt. 42.d0) then
       BoseEinstein = 0.d0
    elseif (freq/kT .lt. -42.d0) then
       BoseEinstein = -1.d0
    elseif (freq/kT .lt. 1.d-15) then
       BoseEinstein = 1.d18
    else
       BoseEinstein = 1.d0 / (DEXP(freq/kT) - 1.d0)
    endif


  end function BoseEinstein


!  *******************************************************************  !


END MODULE idsrdr_distrib


