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
!  PS.: module extracted from i-smeagol.                                !
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

  PUBLIC :: DISTfermiDirac, DISTfermiDiracPM, DISTboseEinstein,         &
            DISTsnglFermiDirac, DISTsnglBoseEinstein,                   &
            DISTsnglZFermiDirac


CONTAINS


!  *******************************************************************  !
!                          DISTsnglZFermiDirac                          !
!  *******************************************************************  !
!  Description: Calculates the Fermi-Dirac distribution for fermions    !
!  in equilibrium with the electrochemical potential 'mu'.              !
!                                                                       !
!  Written by Pedro Brandimarte, Aug 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    August 2013                                     !
!  ****************************** INPUT ******************************  !
!  complex*8 E                : Energy point                            !
!  real*8 mu                  : Electrochemical potential               !
!  real*8 kT                  : Electronic temperature (Ry) times       !
!                               Boltzmann constant                      !
!  ***************************** OUTPUT ******************************  !
!  complex*8 dFermi           : Fermi-Dirac distribution                !
!  *******************************************************************  !
  subroutine DISTsnglZFermiDirac (E, mu, kT, dFermi)

!   Input variables.
    real(8) :: mu, kT
    complex(8) :: E, dFermi

    if ((DREAL(E) - mu) / kT .gt. 40.D0) then
       dFermi = 0.D0
    else
       dFermi = 1.D0 / (CDEXP(DCMPLX((DREAL(E) - mu)) / kT) + 1.D0)
    endif


  end subroutine DISTsnglZFermiDirac


!  *******************************************************************  !
!                          DISTsnglFermiDirac                           !
!  *******************************************************************  !
!  Description: Calculates the Fermi-Dirac distribution for fermions    !
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
!  ***************************** OUTPUT ******************************  !
!  real*8 dFermi              : Fermi-Dirac distribution                !
!  *******************************************************************  !
  subroutine DISTsnglFermiDirac (E, mu, kT, dFermi)

!   Input variables.
    real(8) :: mu, kT
    real(8) :: dFermi
    complex(8) :: E

    if ((DREAL(E) - mu) / kT .gt. 40.D0) then
       dFermi = 0.D0
    else
       dFermi = 1.D0 / (DEXP((DREAL(E) - mu) / kT) + 1.D0)
    endif


  end subroutine DISTsnglFermiDirac


!  *******************************************************************  !
!                            DISTfermiDirac                             !
!  *******************************************************************  !
!  Description: Calculates the Fermi-Dirac distribution for fermions    !
!  in equilibrium with the electrochemical potential 'mu'.              !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2013                                    !
!  ****************************** INPUT ******************************  !
!  integer Nenerg_div         : Local # of energy points                !
!  real*8 Ei(Nenerg_div)      : Energy points                           !
!  real*8 mu                  : Electrochemical potential               !
!  real*8 kT                  : Electronic temperature (Ry) times       !
!                               Boltzmann constant                      !
!  ***************************** OUTPUT ******************************  !
!  real*8 dFermi(Nenerg_div)  : Fermi-Dirac distribution                !
!  *******************************************************************  !
  subroutine DISTfermiDirac (Nenerg_div, Ei, mu, kT, dFermi)

!   Input variables.
    integer :: Nenerg_div
    real(8) :: mu, kT
    real(8) :: dFermi(Nenerg_div)
    complex(8) :: Ei(Nenerg_div)

!   Local variables.
    integer :: e

    do e = 1,Nenerg_div ! over energy grid
       call DISTsnglFermiDirac (Ei(e), mu, kT, dFermi(e))
    enddo


  end subroutine DISTfermiDirac


!  *******************************************************************  !
!                           DISTfermiDiracPM                            !
!  *******************************************************************  !
!  Description: Calculates the Fermi-Dirac distribution for fermions    !
!  in equilibrium with the electrochemical potential 'mu', at all       !
!  energy grid plus ('dFermip') and minus ('dFermim') the vibrational   !
!  mode energy ('freq').                                                !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2013                                    !
!  ****************************** INPUT ******************************  !
!  integer Nenerg_div         : Local # of energy points                !
!  complex*8 Ei(Nenerg_div)   : Energy points                           !
!  integer nModes             : # of vibrational modes                  !
!  real*8 freq(nModes)        : Vibrational mode frequencies (energy)   !
!  real*8 mu                  : Electrochemical potential               !
!  real*8 kT                  : Electronic temperature (Ry) times       !
!                               Boltzmann constant                      !
!  ***************************** OUTPUT ******************************  !
!  real*8 dFermip(Nenerg_div,nModes)   : Fermi-Dirac distribution       !
!                                        at 'E+freq'                    !
!  real*8 dFermim(Nenerg_div,nModes)   : Fermi-Dirac distribution       !
!                                        at 'E+freq'                    !
!  *******************************************************************  !
  subroutine DISTfermiDiracPM (Nenerg_div, Ei, nModes, freq,            &
                               mu, kT, dFermip, dFermim)

!   Input variables.
    integer :: Nenerg_div, nModes
    real(8) :: mu, kT
    real(8) :: freq(nModes)
    real(8) :: dFermip(Nenerg_div,nModes), dFermim(Nenerg_div,nModes)
    complex(8) :: Ei(Nenerg_div)

!   Local variables.
    integer :: lambda, k

    do lambda = 1,nModes ! over the modes
       do k = 1,Nenerg_div ! over energy grid
             call DISTsnglFermiDirac (Ei(k)+freq(lambda), mu, kT,       &
                  dFermip(k,lambda))
             call DISTsnglFermiDirac (Ei(k)-freq(lambda), mu, kT,       &
                  dFermim(k,lambda))
       enddo
    enddo


  end subroutine DISTfermiDiracPM


!  *******************************************************************  !
!                         DISTsnglBoseEinstein                          !
!  *******************************************************************  !
!  Description: Calculates the Bose-Einstein distribution for phonons   !
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
!  ***************************** OUTPUT ******************************  !
!  real*8 phonBEdist          : Bose-Einstein distribution              !
!  *******************************************************************  !
  subroutine DISTsnglBoseEinstein (freq, kT, phonBEdist)

!   Input variables.
    real(8) :: kT
    real(8) :: freq, phonBEdist

!!$    if (freq/kT .gt. 40.D0) then
!!$       phonBEdist = 0.D0
!!$    elseif (freq/kT .lt. 1.D-15) then
!!$       phonBEdist = 1.D18
!!$    else
       phonBEdist = 1.D0 / (DEXP(freq/kT) - 1.D0)
!!$    endif


     end subroutine DISTsnglBoseEinstein


!  *******************************************************************  !
!                           DISTboseEinstein                            !
!  *******************************************************************  !
!  Description: Calculates the Bose-Einstein distribution for phonons   !
!  with energy 'freq' in thermal equilibrium.                           !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2013                                    !
!  ****************************** INPUT ******************************  !
!  integer nModes             : # of vibrational modes                  !
!  real*8 freq(nModes)        : Vibrational mode frequencies (energy)   !
!  real*8 kT                  : Electronic temperature (Ry) times       !
!                               Boltzmann constant                      !
!  ***************************** OUTPUT ******************************  !
!  real*8 phonBEdist(nModes)  : Bose-Einstein distribution              !
!  *******************************************************************  !
  subroutine DISTboseEinstein (nModes, freq, kT, phonBEdist)

!   Input variables.
    integer :: nModes
    real(8) :: kT
    real(8) :: freq(nModes), phonBEdist(nModes)

!   Local variables.
    integer :: lambda

    do lambda = 1,nModes ! over the modes
       call DISTsnglBoseEinstein (freq(lambda), kT, phonBEdist(lambda))
    enddo


  end subroutine DISTboseEinstein


!  *******************************************************************  !


END MODULE idsrdr_distrib


