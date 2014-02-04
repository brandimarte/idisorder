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
!                          MODULE idsrdr_power                          !
!  *******************************************************************  !
!  Description: compute the power dissipated into the phonon system.    !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *******************************************************************  !

MODULE idsrdr_power

!
!   Modules
!
  use idsrdr_units,    only: 
  use idsrdr_options,  only: 
  use idsrdr_engrid,   only: 
  use idsrdr_leads,    only: 
  use idsrdr_ephcoupl, only: 
  use idsrdr_green,    only: 
  use idsrdr_io,       only: 
  use idsrdr_distrib,  only: 

  implicit none
  
  PUBLIC  :: powerinit, powerTr, power, freepower, phonon, phPwr, phOccup
  PRIVATE :: electronhole, effectiveEmission, occupPower,               &
             modePhon, modePwr, h

  TYPE phonon
     sequence
     real(8), pointer :: P(:,:,:,:)
  END TYPE phonon

! Dissipated power per phonon mode (for all units).
  TYPE(phonon), allocatable, dimension (:) :: phPwr

! Steady state phonon occupation for each mode (for all units).
  TYPE(phonon), allocatable, dimension (:) :: phOccup

! Type for storing calculated powers per mode contribution.
  TYPE modePhon
     real(8), pointer :: Pem(:) ! dissipated power by effective emission
     real(8), pointer :: Peh(:) ! dissipated power by e-h damping
  END TYPE modePhon

! Calculated powers per mode (for all units with e-ph).
  TYPE(modePhon), allocatable, dimension (:) :: modePwr

! Constants. (physical constants from CODATA 2013)
  real(8), parameter :: h   = 4.135667516D-15 / 13.60569253D0 ! Ry*s
!!$  real(8), parameter :: h   = 4.135667516D-15 ! eV*s
!!$  real(8), parameter :: Rhc = 13.60569253D0 ! eV


CONTAINS


!  *******************************************************************  !
!                               powerinit                               !
!  *******************************************************************  !
!  Description: allocate array of type 'phonon' for storing dissipated  !
!  power and the phonon occupation.                                     !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer nModes(neph)        : Number of vibrational modes            !
!  integer ephIdx(ntypeunits+2): Unit index (those with e-ph)           !
!  integer nunitseph           : Number of units with eph               !
!  integer eph_type(nunitseph) : Units types with eph                   !
!  integer nspin               : Number of spin components              !
!  integer NIVP                : Number of bias potential points        !
!  integer NTenerg_div         : Number of energy grid points per node  !
!  *******************************************************************  !
  subroutine powerinit

!
!   Modules
!
    use idsrdr_ephcoupl, only: nModes, ephIdx
    use idsrdr_units,    only: nunitseph, eph_type
    use idsrdr_options,  only: nspin, NIVP
    use idsrdr_engrid,   only: NTenerg_div

!   Local variables.
    integer :: u, idx

!   Allocate pointer arrays.
    allocate (phPwr(nunitseph))
    allocate (phOccup(nunitseph))
    allocate (modePwr(nunitseph))

    do u = 1,nunitseph ! over units with e-ph

       idx = ephIdx(eph_type(u))

!      Allocate memory.
       allocate (phPwr(u)%P(NTenerg_div,nspin,NIVP,nModes(idx)))
       allocate (phOccup(u)%P(NTenerg_div,nspin,NIVP,nModes(idx)))
       allocate (modePwr(u)%Peh(nModes(idx)))
       allocate (modePwr(u)%Pem(nModes(idx)))
       phPwr(u)%P = 0.d0
       phOccup(u)%P = 0.d0

    enddo


  end subroutine powerinit


!  *******************************************************************  !
!                                powerTr                                !
!  *******************************************************************  !
!  Description: interface subroutine for computing the bias             !
!  independent part of the phonon power.                                !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer NL                   : Number of left lead orbitals          !
!  integer NR                   : Number of right lead orbitals         !
!  complex(8) Sigma_L(NL,NL)    : Left-lead self-energy                 !
!  complex(8) Sigma_R(NR,NR)    : Right-lead self-energy                !
!  logical writeondisk          : Write GFs on disk?                    !
!  integer nunitseph            : Number of units with eph              !
!  integer eph_type(nunitseph)  : Units types with eph                  !
!  TYPE(green) Gr_nn(nunitseph)%G(unitdimensions,unitdimensions) :      !
!                                                  [complex] G^r_{n,n}  !
!  TYPE(green) Gr_1n(nunitseph)%G(NL,unitdimensions) : [complex]        !
!                                                      G^r_{1,n}        !
!  TYPE(green) Gr_Mn(nunitseph)%G(NR,unitdimensions) : [complex]        !
!                                                      G^r_{M,n}        !
!  ****************************** INPUT ******************************  !
!  integer ispin                : Spin component index                  !
!  *******************************************************************  !
  subroutine powerTr (ispin)

!
!   Modules
!
    use idsrdr_leads,    only: NL, NR, Sigma_L, Sigma_R
    use idsrdr_options,  only: writeondisk
    use idsrdr_units,    only: nunitseph, eph_type
    use idsrdr_ephcoupl, only: nModes, norbDyn, ephIdx
    use idsrdr_green,    only: Gr_nn, Gr_1n, Gr_Mn, greenload,          &
                               Gr_nn_disk, Gr_1n_disk, Gr_Mn_disk
    use idsrdr_io,       only: IOopenStream, IOcloseStream

!   Input variables.
    integer, intent(in) :: ispin

!   Local variables.
    integer :: ueph, idx
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8), allocatable, dimension (:,:) :: Gamma_L, Gamma_R,       &
                                                Gnn, G1n, GMn

!   Allocate coupling and auxiliary matrices.
    allocate (Gamma_L(NL,NL))
    allocate (Gamma_R(NR,NR))

!   Sets the lead's coupling matrices.
    Gamma_L = zi * (Sigma_L - DCONJG(Sigma_L))
    Gamma_R = zi * (Sigma_R - DCONJG(Sigma_R))

    IF (writeondisk) THEN

!      Open Green's functions files.
       call IOopenStream (Gr_nn_disk%fname, Gr_nn_disk%lun)
       call IOopenStream (Gr_1n_disk%fname, Gr_1n_disk%lun)
       call IOopenStream (Gr_Mn_disk%fname, Gr_Mn_disk%lun)

       do ueph = 1,nunitseph ! over unit with e-ph

          idx = ephIdx(eph_type(ueph))

!         Allocate auxiliary matrix.
          allocate (Gnn(norbDyn(idx),norbDyn(idx)))

          call greenload (Gr_nn_disk%lun, norbDyn(idx), norbDyn(idx),   &
                          1, norbDyn(idx), 1, norbDyn(idx), Gnn,        &
                          norbDyn(idx), norbDyn(idx))

!         Compute the electron-hole damping contribution.
          call electronhole (ispin, ueph, idx, nModes(idx),             &
                             norbDyn(idx), Gnn)

!         Free memory.
          deallocate (Gnn)

!         Allocate auxiliary matrices.
          allocate (G1n(NL,norbDyn(idx)))
          allocate (GMn(NR,norbDyn(idx)))

          call greenload (Gr_1n_disk%lun, NL, norbDyn(idx), 1, NL,      &
                          1, norbDyn(idx), G1n, NL, norbDyn(idx))
          call greenload (Gr_Mn_disk%lun, NR, norbDyn(idx), 1, NR,      &
                          1, norbDyn(idx), GMn, NR, norbDyn(idx))

!         Compute the nonequilibrium "effective emission" contribution.
          call effectiveEmission (ispin, ueph, idx, nModes(idx),        &
                                  NL, Gamma_L, NR, Gamma_R,             &
                                  norbDyn(idx), G1n, GMn)

!         Free memory.
          deallocate (G1n)
          deallocate (GMn)

       enddo

!      Close Green's functions files.
       call IOcloseStream (Gr_nn_disk%fname, Gr_nn_disk%lun)
       call IOcloseStream (Gr_1n_disk%fname, Gr_1n_disk%lun)
       call IOcloseStream (Gr_Mn_disk%fname, Gr_Mn_disk%lun)

    ELSE ! write on memory...

       do ueph = 1,nunitseph ! over unit with e-ph

          idx = ephIdx(eph_type(ueph))

!         Compute the electron-hole damping contribution.
          call electronhole (ispin, ueph, idx, nModes(idx),             &
                             norbDyn(idx), Gr_nn(ueph)%G)

          call effectiveEmission (ispin, ueph, idx, nModes(idx),        &
                                  NL, Gamma_L, NR, Gamma_R,             &
                                  norbDyn(idx), Gr_1n(ueph)%G,          &
                                  Gr_Mn(ueph)%G)

       enddo

    ENDIF ! IF (writeondisk)

!   Free memory.
    deallocate (Gamma_L)
    deallocate (Gamma_R)


  end subroutine powerTr


!  *******************************************************************  !
!                                 power                                 !
!  *******************************************************************  !
!  Description: interface subroutine for computing the phonon power     !
!  and occupation.                                                      !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer nunitseph            : Number of units with eph              !
!  integer eph_type(nunitseph)  : Units types with eph                  !
!  integer nModes               : Number of vibrational modes           !
!  integer ephIdx(ntypeunits+2) : Unit index (those with e-ph)          !
!  ****************************** INPUT ******************************  !
!  integer ienergy              : Energy grid index                     !
!  integer ispin                : Spin component index                  !
!  integer iv                   : Bias potential index                  !
!  real*8 Vbias                 : Bias potential                        !
!  *******************************************************************  !
  subroutine power (ienergy, ispin, iv, Vbias)

!
!   Modules
!
    use idsrdr_units,    only: nunitseph, eph_type
    use idsrdr_ephcoupl, only: nModes, ephIdx

!   Input variables.
    integer, intent(in) :: ienergy, ispin, iv
    real(8), intent(in) :: VBias

!   Local variables.
    integer :: ueph, idx

    do ueph = 1,nunitseph ! over unit with e-ph

       idx = ephIdx(eph_type(ueph))

!      Compute the modes occupation and dissipated power.
       call occupPower (ienergy, ispin, iv, ueph, idx,                  &
                        nModes(idx), Vbias)

    enddo


  end subroutine power


!  *******************************************************************  !
!                             electronhole                              !
!  *******************************************************************  !
!  Description: compute the first term of the dissipated power          !
!  expression at Lowest Order Approximation (second) which describes    !
!  the equilibrium energy exchange between the vibrational and          !
!  electronic degrees of freedom (electron-hole pair damping of the     !
!  vibrations).                                                         !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  TYPE(ephCplng) Meph(neph)%M(norbDyn,norbDyn,nspin,nModes) :          !
!                                   [complex*8] E-ph coupling matrices  !
!  ****************************** INPUT ******************************  !
!  integer ispin               : Spin component index                   !
!  integer ueph                : E-ph unit index                        !
!  integer idx                 : Unit type index                        !
!  integer nModes              : Number of vibrational modes            !
!  integer norbDyn             : Number of orbitals from dynamic atoms  !
!  complex*8 Gr_nn(norbDyn,norbDyn) : G^r_{n,n}                         !
!  *******************************************************************  !
  subroutine electronhole (ispin, ueph, idx, nModes, norbDyn, Gr_nn)

!
!   Modules
!
    use idsrdr_ephcoupl, only: Meph

!   Input variables.
    integer, intent(in) :: ispin, ueph, idx, nModes, norbDyn
    complex(8), dimension (norbDyn,norbDyn), intent(in) :: Gr_nn

!   Local variables.
    integer :: w, i
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8), allocatable, dimension(:,:) :: Aux1, Aux2
    external :: HI_zsymm, HI_zgemm

!   Allocate auxiliary matrices.
    allocate (Aux1(norbDyn,norbDyn))
    allocate (Aux2(norbDyn,norbDyn))

!   Spectral matrix (obs.: 'Gr_nn' is symmetric).
    Aux2 = zi * (Gr_nn - DCONJG(Gr_nn))

    do w = 1,nModes ! over phonon modes

!      -- 'Meph*A*Meph*A' --

!      ('Aux1 = Meph * Aux2')
       call HI_zsymm ('L', 'L', norbDyn, norbDyn, (1.d0,0.d0),          &
                      Meph(idx)%M(:,:,ispin,w), norbDyn, Aux2,          &
                      norbDyn, (0.d0,0.d0), Aux1, norbDyn)

!      ('Aux2 = Aux1 * Aux1')
       call HI_zgemm ('N', 'N', norbDyn, norbDyn, norbDyn,              &
                      (1.d0,0.d0), Aux1, norbDyn, Aux1,                 &
                      norbDyn, (0.d0,0.d0), Aux2, norbDyn)

!      Compute the trace.
       modePwr(ueph)%Peh(w) = 0.d0
       do i = 1,norbDyn
          modePwr(ueph)%Peh(w) = modePwr(ueph)%Peh(w) + DREAL(Aux2(i,i))
       enddo

    enddo ! do w = 1,nModes

!   Free memory.
    deallocate (Aux1)
    deallocate (Aux2)


  end subroutine electronhole


!  *******************************************************************  !
!                           effectiveEmission                           !
!  *******************************************************************  !
!  Description: compute the second term of the dissipated power         !
!  expression at Lowest Order Approximation (second) which is a         !
!  nonequilibrium term related an effective emession rate of            !
!  vibrational quanta under finite bias.                                !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  TYPE(ephCplng) Meph(neph)%M(norbDyn,norbDyn,nspin,nModes) :          !
!                                   [complex*8] E-ph coupling matrices  !
!  ****************************** INPUT ******************************  !
!  integer ispin              : Spin component index                    !
!  integer ueph               : E-ph unit index                         !
!  integer idx                : Unit type index                         !
!  integer nModes             : Number of vibrational modes             !
!  integer norbDyn            : Number of orbitals from dynamic atoms   !
!  complex*8 Gr_nn(norbDyn,norbDyn) : G^r_{n,n}                         !
!  integer NL                       : Number of left lead orbitals      !
!  complex*8 Gr_1n(NL,norbDyn)      : G^r_{1,n}                         !
!  integer NR                       : Number of right lead orbitals     !
!  complex*8 Gr_Mn(NR,norbDyn)      : G^r_{M,n}                         !
!  complex*8 Gamma_L(NL,NL)         : Left-lead coupling matrix         !
!  complex*8 Gamma_R(NR,NR)         : Right-lead coupling matrix        !
!  *******************************************************************  !
  subroutine effectiveEmission (ispin, ueph, idx, nModes,               &
                                NL, Gamma_L, NR, Gamma_R,               &
                                norbDyn, Gr_1n, Gr_Mn)

!
!   Modules
!
    use idsrdr_ephcoupl, only: Meph

!   Input variables.
    integer, intent(in) :: ispin, ueph, idx, nModes, NL, NR, norbDyn
    complex(8), dimension (NL,NL), intent(in) :: Gamma_L
    complex(8), dimension (NR,NR), intent(in) :: Gamma_R
    complex(8), dimension (NL,norbDyn), intent(in) :: Gr_1n
    complex(8), dimension (NR,norbDyn), intent(in) :: Gr_Mn

!   Local variables.
    integer :: w, i
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8), allocatable, dimension(:,:) :: Gr1nCJG, GrMnCJG, Aux1,  &
                                               Aux2, Aux3, Aux4, Aux5
    external :: HI_zsymm, HI_zhemm, HI_zgemm

!   Allocate auxiliary matrices.
    allocate (Gr1nCJG(NL,norbDyn))
    allocate (GrMnCJG(NR,norbDyn))
    allocate (Aux1(NL,norbDyn))
    allocate (Aux2(norbDyn,norbDyn))
    allocate (Aux3(norbDyn,norbDyn))
    allocate (Aux4(NR,norbDyn))
    allocate (Aux5(norbDyn,norbDyn))

!   Copy the complex conjugate of 'Gr_1n' and 'Gr_Mn'.
    Gr1nCJG = DCONJG(Gr_1n)
    GrMnCJG = DCONJG(Gr_Mn)

    do w = 1,nModes ! over phonon modes

!      -- 1st PART: 'Meph*G*Gamma_L*G^dagger' --

!      ('Aux1 = Gamma_L * Gr_1n^*')
       call HI_zhemm ('L', 'L', NL, norbDyn, (1.d0,0.d0), Gamma_L, NL,  &
                      Gr1nCJG, NL, (0.d0,0.d0), Aux1, NL)

!      ('Aux2 = Gr_1n^T * Aux1')
       call HI_zgemm ('T', 'N', norbDyn, norbDyn, NL, (1.d0,0.d0),      &
                      Gr_1n, NL, Aux1, NL, (0.d0,0.d0), Aux2, norbDyn)

!      ('Aux3 = Meph * Aux2')
       call HI_zsymm ('L', 'L', norbDyn, norbDyn, (1.d0,0.d0),          &
                      Meph(idx)%M(:,:,ispin,w), norbDyn, Aux2,          &
                      norbDyn, (0.d0,0.d0), Aux3, norbDyn)

!      -- 2nd PART: 'Meph*G*Gamma_R*G^dagger' --

!      ('Aux4 = Gamma_R * Gr_Mn^*')
       call HI_zhemm ('L', 'L', NR, norbDyn, (1.d0,0.d0), Gamma_R, NR,  &
                      GrMnCJG, NR, (0.d0,0.d0), Aux4, NR)

!      ('Aux2 = Gr_Mn^T * Aux4')
       call HI_zgemm ('T', 'N', norbDyn, norbDyn, NR, (1.d0,0.d0),      &
                      Gr_Mn, NR, Aux4, NR, (0.d0,0.d0), Aux2, norbDyn)

!      ('Aux5 = Meph * Aux2')
       call HI_zsymm ('L', 'L', norbDyn, norbDyn, (1.d0,0.d0),          &
                      Meph(idx)%M(:,:,ispin,w), norbDyn, Aux2,          &
                      norbDyn, (0.d0,0.d0), Aux5, norbDyn)

!      -- 1st PART * 2nd PART --

       call HI_zgemm ('N', 'N', norbDyn, norbDyn, norbDyn,              &
                      (1.d0,0.d0), Aux3, norbDyn, Aux5,                 &
                      norbDyn, (0.d0,0.d0), Aux2, norbDyn)

!      Compute the trace.
       modePwr(ueph)%Pem(w) = 0.d0
       do i = 1,norbDyn
          modePwr(ueph)%Pem(w) = modePwr(ueph)%Pem(w) + DREAL(Aux2(i,i))
       enddo

    enddo ! do w = 1,nModes

!   Free memory.
    deallocate (Gr1nCJG)
    deallocate (GrMnCJG)
    deallocate (Aux1)
    deallocate (Aux2)
    deallocate (Aux3)
    deallocate (Aux4)
    deallocate (Aux5)


  end subroutine effectiveEmission


!  *******************************************************************  !
!                              occupPower                               !
!  *******************************************************************  !
!  Description: compute the bias-dependent modes distribution and the   !
!  dissipated power by electrons to each phonon mode.                   !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  TYPE(ephFreq) freq(neph)%F(nModes) : [real*8] Vibrational mode's     !
!                                       frequencies                     !
!  real*8 temp                        : Electronic temperature          !
!  logical phonEqui                   : Phonons in equilibrium?         !
!  real*8 phonDamp                    : Damping parameter               !
!  ****************************** INPUT ******************************  !
!  integer ienergy            : Energy grid index                       !
!  integer ispin              : Spin component index                    !
!  integer iv                 : Bias potential index                    !
!  integer ueph               : Unit index                              !
!  integer idx                : Unit type index                         !
!  integer nModes             : Number of vibrational modes             !
!  real*8 Vbias               : Bias potential value                    !
!  *******************************************************************  !
  subroutine occupPower (ienergy, ispin, iv, ueph, idx, nModes, Vbias)

!
!   Modules
!
    use idsrdr_ephcoupl, only: freq
    use idsrdr_options,  only: temp
    use idsrdr_options,  only: phonEqui, phonDamp
    use idsrdr_distrib,  only: BoseEinstein

!   Input variables.
    integer, intent(in) :: ienergy, ispin, iv, ueph, idx, nModes
    real(8), intent(in) :: Vbias

!   Local variables.
    integer :: w
    real(8) :: coshVT, sinhVT, foo, bose, Pem, Peh

!   Compute bias and temperature dependent factors.
    coshVT = DCOSH (Vbias / temp)
    sinhVT = DSINH (Vbias / temp)

    do w = 1,nModes ! over phonon modes

!      Compute effective-emission pre-factor.
       foo = freq(idx)%F(w) * (coshVT - 1.d0) /                         &
            DTANH(freq(idx)%F(w)/(2.d0*temp))                           &
            - Vbias * sinhVT
       foo = foo / (h * (DCOSH(freq(idx)%F(w)/temp) - coshVT))
       Pem = modePwr(ueph)%Pem(w) * foo

!      Compute eletron-hole pre-factor.
       Peh = modePwr(ueph)%Peh(w) * freq(idx)%F(w) / h

!      Equilibrium distribution.
       bose = BoseEinstein (freq(idx)%F(w), temp)

!      Compute phonon occupation.
       if (phonEqui) then
          phOccup(ueph)%P(ienergy,ispin,iv,w) =  bose
       else
          phOccup(ueph)%P(ienergy,ispin,iv,w) =  bose + Pem /           &
               (freq(idx)%F(w) * (Peh + phonDamp))
       endif

!      Compute dissipated power.
       phPwr(ueph)%P(ienergy,ispin,iv,w) = freq(idx)%F(w) *             &
            ((bose - phOccup(ueph)%P(ienergy,ispin,iv,w)) * Peh + Pem)

    enddo


  end subroutine occupPower


!  *******************************************************************  !
!                               freepower                               !
!  *******************************************************************  !
!  Description: free allocated memory.                                  !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer nunitseph           : Number of units with eph               !
!  *******************************************************************  !
  subroutine freepower

!
!   Modules
!
    use idsrdr_units,    only: nunitseph

!   Local variables.
    integer :: u

!   First deallocates pointed arrays.
    do u = 1,nunitseph ! over units with e-ph
       deallocate (phPwr(u)%P)
       deallocate (phOccup(u)%P)
       deallocate (modePwr(u)%Peh)
       deallocate (modePwr(u)%Pem)
    enddo
    deallocate (phPwr)
    deallocate (phOccup)
    deallocate (modePwr)


  end subroutine freepower


!  *******************************************************************  !


END MODULE idsrdr_power
