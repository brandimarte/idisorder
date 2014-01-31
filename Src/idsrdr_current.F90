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
!                         MODULE idsrdr_current                         !
!  *******************************************************************  !
!  Description: compute the eletronic current.                          !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *******************************************************************  !

MODULE idsrdr_current

!
!   Modules
!
  use idsrdr_options,  only: 
  use idsrdr_engrid,   only: 
  use idsrdr_leads,    only: 
  use idsrdr_ephcoupl, only: 
  use idsrdr_green,    only: 
  use idsrdr_units,    only: 
  use idsrdr_recipes,  only: 
  use idsrdr_distrib,  only: 
  use idsrdr_hilbert,  only: 
  use idsrdr_io,       only: 
  use idsrdr_power,    only: 

  implicit none
  
  PUBLIC  :: currentinit, current, calcCurr, allcurr, freecurr,         &
             sumCalcCurr
  PRIVATE :: elastic, transmission, inelSymm, inelSymmDisk, asymmPre,   &
             inelAsymm, inelAsymmDisk, testInelSymm, testInelAsymm,     &
             kbTol, eoverh

! Type for storing calculated currents.
  TYPE calcCurr
     sequence
     real(8) :: el ! elastic part
     real(8) :: isymm ! inelastic symmetric
     real(8) :: iasymm ! inelastic asymmetric
  END TYPE calcCurr

! Calculated currents.
  TYPE(calcCurr), allocatable, dimension (:,:,:) :: allcurr

  real(8), parameter :: kbTol = 18.d0 ! tolerance value for temperature

! Constants. (physical constants from CODATA 2013)
!!$  real(8), parameter :: e   = 1.602176565D-19 ! C
!!$  real(8), parameter :: h   = 4.135667516D-15 ! eV*s
!!$  real(8), parameter :: Rhc = 13.60569253D0 ! eV
  real(8), parameter :: eoverh = 1.602176565D-4 * 13.60569253D0 /       &
                                 4.135667516D0

CONTAINS


!  *******************************************************************  !
!                              currentinit                              !
!  *******************************************************************  !
!  Description: allocate array of type 'calcCurr' for storing           !
!  calculated currents.                                                 !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer nspin               : Number of spin components              !
!  integer NIVP                : Number of bias potential points        !
!  integer NTenerg_div         : Number of energy grid points per node  !
!  *******************************************************************  !
  subroutine currentinit

!
!   Modules
!
    use idsrdr_options,  only: nspin, NIVP
    use idsrdr_engrid,   only: NTenerg_div

!   Allocate and initializes current array.
    allocate (allcurr(NTenerg_div,nspin,NIVP))
    allcurr%el = 0.d0
    allcurr%isymm = 0.d0
    allcurr%iasymm = 0.d0


  end subroutine currentinit


!  *******************************************************************  !
!                                current                                !
!  *******************************************************************  !
!  Description: interface subroutine for computing the current.         !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer NL                   : Number of left lead orbitals          !
!  integer NR                   : Number of right lead orbitals         !
!  complex(8) Sigma_L(NL,NL)    : Left-lead self-energy                 !
!  complex(8) Sigma_R(NR,NR)    : Right-lead self-energy                !
!  logical writeondisk          : Write GFs on disk?                    !
!  ****************************** INPUT ******************************  !
!  real*8 Ei                    : Energy grid point                     !
!  integer ienergy              : Energy grid index                     !
!  integer ispin                : Spin component index                  !
!  integer iv                   : Bias potential index                  !
!  real*8 Vbias                 : Bias potential value                  !
!  *******************************************************************  !
  subroutine current (Ei, ienergy, ispin, iv, Vbias)

!
!   Modules
!
    use idsrdr_leads,    only: NL, NR, Sigma_L, Sigma_R
    use idsrdr_options,  only: writeondisk

!   Input variables.
    integer, intent(in) :: ienergy, ispin, iv
    real(8), intent(in) :: Ei, Vbias

!   Local variables.
    real(8) :: Iel, Isymm, Iasymm
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8), allocatable, dimension (:,:) :: Gamma_L
    complex(8), allocatable, dimension (:,:) :: Gamma_R

!   Allocate coupling and auxiliary matrices.
    allocate (Gamma_L(NL,NL))
    allocate (Gamma_R(NR,NR))

!   Sets the lead's coupling matrices.
    Gamma_L = zi * (Sigma_L - DCONJG(Sigma_L))
    Gamma_R = zi * (Sigma_R - DCONJG(Sigma_R))

    IF (writeondisk) THEN

!      Compute elastic contribution.
       call elastic (Iel, NL, Gamma_L, NR, Gamma_R, Vbias)

!      Compute symmetric part of inelastic contribution.
#ifdef DEBUG
       call inelSymmDisk (Isymm, ienergy, ispin, iv, NL, Gamma_L,       &
                          NR, Gamma_R, Vbias, Ei)
#else
       call inelSymmDisk (Isymm, ienergy, ispin, iv, NL, Gamma_L,       &
                          NR, Gamma_R, Vbias)
#endif

!      Compute asymmetric part of inelastic contribution.
       call inelAsymmDisk (Iasymm, ispin, NL, Gamma_L,                  &
                           NR, Gamma_R, Vbias, Ei)

!      Store calculated currents.
       allcurr(ienergy,ispin,iv)%el = Iel
       allcurr(ienergy,ispin,iv)%isymm = Isymm
       allcurr(ienergy,ispin,iv)%iasymm = Iasymm

    ELSE ! write on memory...

!      Compute elastic contribution.
       call elastic (Iel, NL, Gamma_L, NR, Gamma_R, Vbias)

!      Compute symmetric part of inelastic contribution.
#ifdef DEBUG
       call inelSymm (Isymm, ienergy, ispin, iv, NL, Gamma_L,           &
                      NR, Gamma_R, Vbias, Ei)
#else
       call inelSymm (Isymm, ienergy, ispin, iv, NL, Gamma_L,           &
                      NR, Gamma_R, Vbias)
#endif

!      Compute asymmetric part of inelastic contribution.
       call inelAsymm (Iasymm, ispin, NL, Gamma_L,                      &
                       NR, Gamma_R, Vbias, Ei)

!      Store calculated currents.
       allcurr(ienergy,ispin,iv)%el = Iel
       allcurr(ienergy,ispin,iv)%isymm = Isymm
       allcurr(ienergy,ispin,iv)%iasymm = Iasymm

    ENDIF ! IF (writeondisk)

!   Free memory.
    deallocate (Gamma_L)
    deallocate (Gamma_R)


  end subroutine current


!  *******************************************************************  !
!                                elastic                                !
!  *******************************************************************  !
!  Description: compute the elastic component of the current.           !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  complex(8) Gr_1M(NL,NR)             : G^r_{1,M}                      !
!  ****************************** INPUT ******************************  !
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Gamma_L(NL,NL)           : Left-lead coupling matrix      !
!  complex(8) Gamma_R(NR,NR)           : Right-lead coupling matrix     !
!  real*8 Vbias                        : Bias potential                 !
!  ***************************** OUTPUT ******************************  !
!  real*8 Iel                          : Elastic current                !
!  *******************************************************************  !
  subroutine elastic (Iel, NL, Gamma_L, NR, Gamma_R, Vbias)

!
!   Modules
!
    use idsrdr_green,    only: Gr_1M

!   Input variables.
    integer, intent(in) :: NL, NR
    real(8), intent(in) :: Vbias
    real(8), intent(out) :: Iel
    complex(8), dimension (NL,NL), intent(in) :: Gamma_L
    complex(8), dimension (NR,NR), intent(in) :: Gamma_R

!   Calculates the transmission coefficient.
    call transmission (NL, Gamma_L, NR, Gamma_R, Gr_1M, Iel)

    Iel = eoverh * Vbias * Iel


  end subroutine elastic


!  *******************************************************************  !
!                             transmission                              !
!  *******************************************************************  !
!  Description: calculates the Landauer-Buttiker transmission           !
!  coefficient between probes I and J (for example, the electrodes L    !
!  and R).                                                              !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  ****************************** INPUT ******************************  !
!  integer NI                  : Number of basis orbitals from probe I  !
!  complex*8 Gamma_I(NI,NI)    : Low triangular I-probe coupling        !
!                                matrix (for a given spin)              !
!  integer NJ                  : Number of basis orbitals in the right  !
!                                lead, including spin components        !
!  complex*8 Gamma_J(NJ,NJ)    : Low triangular J-probe coupling        !
!                                matrix (for a given spin)              !
!  complex*8 Gr_ij(NI,NJ)      : I-J probes part of unperturbed         !
!                                retarded Green's function              !
!  ***************************** OUTPUT ******************************  !
!  complex*8 Tij               : Transmission coefficient               !
!  *******************************************************************  !
  subroutine transmission (NI, Gamma_I, NJ, Gamma_J, Gr_ij, Tij)


!   Input variables.
    integer, intent(in) :: NI, NJ
    real(8), intent(out) :: Tij
    complex(8), dimension (NI,NI), intent(in) :: Gamma_I
    complex(8), dimension (NJ,NJ), intent(in) :: Gamma_J
    complex(8), dimension (NI,NJ), intent(in) :: Gr_ij

!   Local variables.
    integer :: r, c
    complex(8), dimension (:,:), allocatable :: IJ, JI
    external :: HI_zhemm

!   Allocate auxiliary matrices.
    allocate (IJ(NI,NJ))
    allocate (JI(NJ,NI))

!   ('IJ = Gamma_I * Gr_ij')
    call HI_zhemm ('L', 'L', NI, NJ, (1.d0,0.d0), Gamma_I, NI,          &
                   Gr_ij, NI, (0.d0,0.d0), IJ, NI)

!   ('JI^dagger = Gr_ij * Gamma_J')
    call HI_zhemm ('R', 'L', NI, NJ, (1.d0,0.d0), Gamma_J, NJ,          &
                   Gr_ij, NI, (0.d0,0.d0), JI, NI)

!   Calculates the transmission coefficient.
    Tij = 0.d0
    do c = 1,NJ
       do r = 1,NI
          Tij = Tij + DREAL(IJ(r,c) * DCONJG(JI(r,c)))
       enddo
    enddo

!   Free memory.
    deallocate (IJ)
    deallocate (JI)


  end subroutine transmission


!  *******************************************************************  !
!                               inelSymm                                !
!  *******************************************************************  !
!  Description: compute the symmetric part of inelastic transmission.   !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer nModes(neph)        : Number of vibrational modes            !
!  integer norbDyn(neph)       : Number of orbitals from dynamic atoms  !
!  TYPE(ephCplng) Meph(neph)%M(norbDyn,norbDyn,nspin,nModes) :          !
!                                   [complex*8] E-ph coupling matrices  !
!  TYPE(ephFreq) freq(neph)%F(nModes) : [real*8] Vibrational mode's     !
!                                       frequencies                     !
!  integer ephIdx(ntypeunits+2)       : Unit index (those with e-ph)    !
!  integer nunitseph                  : Number of units with eph        !
!  integer eph_type(nunitseph)        : Units types with eph            !
!  TYPE(green) Gr_Mn(nunitseph)%G(NR,unitdimensions) : [complex]        !
!                                                      G^r_{M,n}        !
!  TYPE(green) Gr_1n(nunitseph)%G(NL,unitdimensions) : [complex]        !
!                                                      G^r_{1,n}        !
!  TYPE(green) Gr_nn(nunitseph)%G(unitdimensions,unitdimensions) :      !
!                                                  [complex] G^r_{n,n}  !
!  complex(8) Gr_1M(NL,NR)                           : G^r_{1,M}        !
!  real*8 temp                 : Electronic temperature                 !
!  TYPE(phonon) phOccup(nunitseph)%P(NTenerg_div,nspin,NIVP,nModes)     !
!                              : [real] Phonon occupation               !
!  ****************************** INPUT ******************************  !
!  integer ienergy                     : Energy grid index              !
!  integer ispin                       : Spin component index           !
!  integer iv                          : Bias potential index           !
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  complex*8 Gamma_L(NL,NL)            : Left-lead coupling matrix      !
!  complex*8 Gamma_R(NR,NR)            : Right-lead coupling matrix     !
!  real*8 Vbias                        : Bias potential                 !
!  real*8 Ei                           : [DEBUG] Energy grid point      !
!  ***************************** OUTPUT ******************************  !
!  real*8 Isymm                : Symmetric part of inelastic current    !
!  *******************************************************************  !
  subroutine inelSymm (Isymm, ienergy, ispin, iv, NL, Gamma_L,          &
#ifdef DEBUG
                       NR, Gamma_R, Vbias, Ei)
#else
                       NR, Gamma_R, Vbias)
#endif

!
!   Modules
!
    use idsrdr_ephcoupl, only: nModes, norbDyn, Meph, freq, ephIdx
    use idsrdr_units,    only: nunitseph, eph_type
    use idsrdr_green,    only: Gr_Mn, Gr_1n, Gr_nn, Gr_1M
    use idsrdr_options,  only: temp
    use idsrdr_distrib,  only: BoseEinstein
    use idsrdr_power,    only: phOccup

!   Input variables.
    integer, intent(in) :: ienergy, ispin, iv, NL, NR
    real(8), intent(out) :: Isymm
    real(8), intent(in) :: Vbias
#ifdef DEBUG
    real(8), intent(in) :: Ei
#endif
    complex(8), dimension (NL,NL), intent(in) :: Gamma_L
    complex(8), dimension (NR,NR), intent(in) :: Gamma_R

!   Local variables.
    integer :: j, w, i, idx
    real(8) :: foo
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8), allocatable, dimension(:,:) :: Aux1, Aux2, Aux3, Aux4,  &
                                               Aux5, Aux6, Aux7, GrCJG, A
    external :: HI_zsymm, HI_zhemm, HI_zgemm

!   Initialize variable.
    Isymm = 0.d0

    do j = 1,nunitseph ! over unit with e-ph

       idx = ephIdx(eph_type(j))

!      Allocate auxiliary matrices.
       allocate (GrCJG(NR,norbDyn(idx)))
       allocate (A(norbDyn(idx),norbDyn(idx)))
       allocate (Aux1(NR,norbDyn(idx)))
       allocate (Aux2(NR,norbDyn(idx)))
       allocate (Aux3(norbDyn(idx),norbDyn(idx)))
       allocate (Aux4(norbDyn(idx),norbDyn(idx)))
       allocate (Aux5(NL,norbDyn(idx)))
       allocate (Aux6(NL,norbDyn(idx)))
       allocate (Aux7(NR,NR))

!      Copy the complex conjugate of 'Gr_Mn'.
       GrCJG = DCONJG(Gr_Mn(j)%G)

!      Spectral matrix (obs.: 'Gr_nn' is symmetric).
       A = zi * (Gr_nn(j)%G - DCONJG(Gr_nn(j)%G))

       do w = 1,nModes(idx) ! over phonon modes

!         -- 1st PART: 'G*Meph*G*Gamma_R*G^dagger*Meph' --

!         ('Aux1 = Gamma_R * Gr_Mn^*')
          call HI_zhemm ('L', 'L', NR, norbDyn(idx), (1.d0,0.d0),       &
                         Gamma_R, NR, GrCJG, NR, (0.d0,0.d0), Aux1, NR)

!         ('Aux2 = Aux1 * Meph')
          call HI_zsymm ('R', 'L', NR, norbDyn(idx), (1.d0,0.d0),       &
                         Meph(idx)%M(:,:,ispin,w), norbDyn(idx),        &
                         Aux1, NR, (0.d0,0.d0), Aux2, NR)

!         ('Aux3 = Gr_Mn^T * Aux2')
          call HI_zgemm ('T', 'N', norbDyn(idx), norbDyn(idx), NR,      &
                         (1.d0,0.d0), Gr_Mn(j)%G, NR, Aux2, NR,         &
                         (0.d0,0.d0), Aux3, norbDyn(idx))

!         ('Aux4 = Meph * Aux3')
          call HI_zsymm ('L', 'L', norbDyn(idx), norbDyn(idx),          &
                         (1.d0,0.d0), Meph(idx)%M(:,:,ispin,w),         &
                         norbDyn(idx), Aux3, norbDyn(idx), (0.d0,0.d0), &
                         Aux4, norbDyn(idx))

!         ('Aux5 = Gr_1n * Aux4')
          call HI_zgemm ('N', 'N', NL, norbDyn(idx), norbDyn(idx),      &
                         (1.d0,0.d0), Gr_1n(j)%G, NL, Aux4,             &
                         norbDyn(idx), (0.d0,0.d0), Aux5, NL)

!         -- 2nd PART: 'Gamma_R*G^dagger*Meph*A*Meph' --

!         ('Aux1 = Aux2 * A') (where 'Aux2 = Gamma_R * Gr_Mn^* * Meph')
          call HI_zhemm ('R', 'L', NR, norbDyn(idx), (1.d0,0.d0), A,    &
                         norbDyn(idx), Aux2, NR, (0.d0,0.d0), Aux1, NR)

!         ('Aux2 = Aux1 * Meph')
          call HI_zsymm ('R', 'L', NR, norbDyn(idx), (1.d0,0.d0),       &
                         Meph(idx)%M(:,:,ispin,w), norbDyn(idx),        &
                         Aux1, NR, (0.d0,0.d0), Aux2, NR)

!         -- 3rd PART: 'G^dagger*Gamma_L*(1st PART + i/2*G*2nd PART)' --

!         ('Aux5 = i/2 * Gr_1M * Aux2 + Aux5')
          call HI_zgemm ('N', 'N', NL, norbDyn(idx), NR, (0.d0,0.5d0),  &
                         Gr_1M, NL, Aux2, NR, (1.d0,0.d0), Aux5, NL)

!         ('Aux6 = Gamma_L * Aux5')
          call HI_zhemm ('L', 'L', NL, norbDyn(idx), (1.d0,0.d0),       &
                         Gamma_L, NL, Aux5, NL, (0.d0,0.d0), Aux6, NL)

!         ('Aux4 = Gr_1n^dagger * Aux6')
          call HI_zgemm ('C', 'N', norbDyn(idx), norbDyn(idx), NL,      &
                         (1.d0,0.d0), Gr_1n(j)%G, NL, Aux6, NL,         &
                         (0.d0,0.d0), Aux4, norbDyn(idx))

!         -- 4rd PART: 'G^dagger*Gamma_L*G*(-i/2 * H.c.)' --

!         ('Aux6 = Gamma_L * Gr_1n')
          call HI_zhemm ('L', 'L', NL, norbDyn(idx), (1.d0,0.d0),       &
                         Gamma_L, NL, Gr_1n(j)%G, NL, (0.d0,0.d0),      &
                         Aux6, NL)

!         ('Aux1 = Gr_1M^dagger * Aux6')
          call HI_zgemm ('C', 'N', NR, norbDyn(idx), NL, (1.d0,0.d0),   &
                         Gr_1M, NL, Aux6, NL, (0.d0,0.d0), Aux1, NR)

!         ('Aux7 = Aux1 * (-i/2 * H.c.)') (where 'H.c. = Aux2^dagger')
          call HI_zgemm ('N', 'C', NR, NR, norbDyn(idx), (0.d0,-0.5d0), &
                         Aux1, NR, Aux2, NR, (0.d0,0.d0), Aux7, NR)

!         Compute the trace.
          do i = 1,norbDyn(idx)
             Isymm = Isymm + DREAL(Aux4(i,i))
          enddo
          do i = 1,NR
             Isymm = Isymm + DREAL(Aux7(i,i))
          enddo
          
#ifdef DEBUG
!         [test] Compute the matrices multiplication with full matrices.
          call testInelSymm (Ei, ispin, NL, Gamma_L,                    &
                             NR, Gamma_R, j, idx, norbDyn(idx),         &
                             Meph(idx)%M(:,:,ispin,w), Aux4, Aux7)
#endif

!         Compute symmetric pre-factor.
          foo = 2.d0 * Vbias * phOccup(j)%P(ienergy,ispin,iv,w)
          foo = foo + (freq(idx)%F(w) - Vbias)                          &
                * BoseEinstein (freq(idx)%F(w) - Vbias, temp)
          foo = foo - (freq(idx)%F(w) + Vbias)                          &
                * BoseEinstein (freq(idx)%F(w) + Vbias, temp)
          Isymm = eoverh * foo * Isymm

       enddo ! do w = 1,nModes(idx)

!      Free memory.
       deallocate (GrCJG)
       deallocate (A)
       deallocate (Aux1)
       deallocate (Aux2)
       deallocate (Aux3)
       deallocate (Aux4)
       deallocate (Aux5)
       deallocate (Aux6)
       deallocate (Aux7)

    enddo ! do j = 1,nunitseph


  end subroutine inelSymm


!  *******************************************************************  !
!                             inelSymmDisk                              !
!  *******************************************************************  !
!  Description: compute the symmetric part of inelastic transmission    !
!  (Green's functions are read from disk).                              !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer nModes(neph)        : Number of vibrational modes            !
!  integer norbDyn(neph)       : Number of orbitals from dynamic atoms  !
!  TYPE(ephCplng) Meph(neph)%M(norbDyn,norbDyn,nspin,nModes) :          !
!                                   [complex*8] E-ph coupling matrices  !
!  TYPE(ephFreq) freq(neph)%F(nModes) : [real*8] Vibrational mode's     !
!                                       frequencies                     !
!  integer ephIdx(ntypeunits+2)       : Unit index (those with e-ph)    !
!  integer nunitseph                  : Number of units with eph        !
!  integer eph_type(nunitseph)        : Units types with eph            !
!  TYPE(green) Gr_Mn(nunitseph)%G(NR,unitdimensions) : [complex]        !
!                                                      G^r_{M,n}        !
!  TYPE(green) Gr_1n(nunitseph)%G(NL,unitdimensions) : [complex]        !
!                                                      G^r_{1,n}        !
!  TYPE(green) Gr_nn(nunitseph)%G(unitdimensions,unitdimensions) :      !
!                                                  [complex] G^r_{n,n}  !
!  complex(8) Gr_1M(NL,NR)                           : G^r_{1,M}        !
!  real*8 temp                 : Electronic temperature                 !
!  TYPE(phonon) phOccup(nunitseph)%P(NTenerg_div,nspin,NIVP,nModes)     !
!                              : [real] Phonon occupation               !
!  ****************************** INPUT ******************************  !
!  integer ienergy                     : Energy grid index              !
!  integer ispin                       : Spin component index           !
!  integer iv                          : Bias potential index           !
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Gamma_L(NL,NL)           : Left-lead coupling matrix      !
!  complex(8) Gamma_R(NR,NR)           : Right-lead coupling matrix     !
!  real*8 Vbias                        : Bias potential                 !
!  real*8 Ei                           : [optional] Energy grid point   !
!  ***************************** OUTPUT ******************************  !
!  real*8 Isymm                : Symmetric part of inelastic current    !
!  *******************************************************************  !
  subroutine inelSymmDisk (Isymm, ienergy, ispin, iv, NL, Gamma_L,      &
#ifdef DEBUG
                           NR, Gamma_R, Vbias, Ei)
#else
                           NR, Gamma_R, Vbias)
#endif

!
!   Modules
!
    use idsrdr_ephcoupl, only: nModes, norbDyn, Meph, freq, ephIdx
    use idsrdr_units,    only: nunitseph, eph_type
    use idsrdr_green,    only: Gr_Mn_disk, Gr_1n_disk, Gr_nn_disk,      &
                               Gr_1M, greenload
    use idsrdr_options,  only: temp
    use idsrdr_distrib,  only: BoseEinstein
    use idsrdr_io,       only: IOopenStream, IOcloseStream
    use idsrdr_power,    only: phOccup

!   Input variables.
    integer, intent(in) :: ienergy, ispin, iv, NL, NR
    real(8), intent(out) :: Isymm
    real(8), intent(in) :: Vbias
#ifdef DEBUG
    real(8), optional, intent(in) :: Ei
#endif
    complex(8), dimension (NL,NL), intent(in) :: Gamma_L
    complex(8), dimension (NR,NR), intent(in) :: Gamma_R

!   Local variables.
    integer :: j, w, i, idx
    real(8) :: foo
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8), allocatable, dimension(:,:) :: Aux1, Aux2, Aux3, Aux4,  &
                                               Aux5, Aux6, Aux7, GrCJG, &
                                               A, auxGr_1n, auxGr_Mn
    external :: HI_zsymm, HI_zhemm, HI_zgemm

!   Initialize variable.
    Isymm = 0.d0

!   Open Green's functions files.
    call IOopenStream (Gr_nn_disk%fname, Gr_nn_disk%lun)
    call IOopenStream (Gr_1n_disk%fname, Gr_1n_disk%lun)
    call IOopenStream (Gr_Mn_disk%fname, Gr_Mn_disk%lun)

    do j = 1,nunitseph ! over unit with e-ph

       idx = ephIdx(eph_type(j))

!      Allocate auxiliary matrices.
       allocate (GrCJG(NR,norbDyn(idx)))
       allocate (A(norbDyn(idx),norbDyn(idx)))
       allocate (Aux1(NR,norbDyn(idx)))
       allocate (Aux2(NR,norbDyn(idx)))
       allocate (Aux3(norbDyn(idx),norbDyn(idx)))
       allocate (Aux4(norbDyn(idx),norbDyn(idx)))
       allocate (Aux5(NL,norbDyn(idx)))
       allocate (Aux6(NL,norbDyn(idx)))
       allocate (Aux7(NR,NR))

!      Spectral matrix (obs.: 'Gr_nn' is symmetric).
       call greenload (Gr_nn_disk%lun, norbDyn(idx), norbDyn(idx),      &
                       1, norbDyn(idx), 1, norbDyn(idx), Aux3,          &
                       norbDyn(idx), norbDyn(idx))
       A = zi * (Aux3 - DCONJG(Aux3))

!      Copy 'Gr_1n' from file to auxiliary matrix.
       allocate (auxGr_1n(NL,norbDyn(idx)))
       call greenload (Gr_1n_disk%lun, NL, norbDyn(idx), 1, NL,         &
                       1, norbDyn(idx), auxGr_1n, NL, norbDyn(idx))

!      Copy 'Gr_Mn' from file to auxiliary matrix.
       allocate (auxGr_Mn(NR,norbDyn(idx)))
       call greenload (Gr_Mn_disk%lun, NR, norbDyn(idx), 1, NR,         &
                       1, norbDyn(idx), auxGr_Mn, NR, norbDyn(idx))

!      Copy the complex conjugate of 'Gr_Mn'.
       GrCJG = DCONJG(auxGr_Mn)

       do w = 1,nModes(idx) ! over phonon modes

!         -- 1st PART: 'G*Meph*G*Gamma_R*G^dagger*Meph' --

!         ('Aux1 = Gamma_R * Gr_Mn^*')
          call HI_zhemm ('L', 'L', NR, norbDyn(idx), (1.d0,0.d0),       &
                         Gamma_R, NR, GrCJG, NR, (0.d0,0.d0), Aux1, NR)

!         ('Aux2 = Aux1 * Meph')
          call HI_zsymm ('R', 'L', NR, norbDyn(idx), (1.d0,0.d0),       &
                         Meph(idx)%M(:,:,ispin,w), norbDyn(idx),        &
                         Aux1, NR, (0.d0,0.d0), Aux2, NR)

!         ('Aux3 = Gr_Mn^T * Aux2')
          call HI_zgemm ('T', 'N', norbDyn(idx), norbDyn(idx), NR,      &
                         (1.d0,0.d0), auxGr_Mn, NR, Aux2, NR,           &
                         (0.d0,0.d0), Aux3, norbDyn(idx))

!         ('Aux4 = Meph * Aux3')
          call HI_zsymm ('L', 'L', norbDyn(idx), norbDyn(idx),          &
                         (1.d0,0.d0), Meph(idx)%M(:,:,ispin,w),         &
                         norbDyn(idx), Aux3, norbDyn(idx), (0.d0,0.d0), &
                         Aux4, norbDyn(idx))

!         ('Aux5 = Gr_1n * Aux4')
          call HI_zgemm ('N', 'N', NL, norbDyn(idx), norbDyn(idx),      &
                         (1.d0,0.d0), auxGr_1n, NL, Aux4,               &
                         norbDyn(idx), (0.d0,0.d0), Aux5, NL)

!         -- 2nd PART: 'Gamma_R*G^dagger*Meph*A*Meph' --

!         ('Aux1 = Aux2 * A') (where 'Aux2 = Gamma_R * Gr_Mn^* * Meph')
          call HI_zhemm ('R', 'L', NR, norbDyn(idx), (1.d0,0.d0), A,    &
                         norbDyn(idx), Aux2, NR, (0.d0,0.d0), Aux1, NR)

!         ('Aux2 = Aux1 * Meph')
          call HI_zsymm ('R', 'L', NR, norbDyn(idx), (1.d0,0.d0),       &
                         Meph(idx)%M(:,:,ispin,w), norbDyn(idx),        &
                         Aux1, NR, (0.d0,0.d0), Aux2, NR)

!         -- 3rd PART: 'G^dagger*Gamma_L*(1st PART + i/2*G*2nd PART)' --

!         ('Aux5 = i/2 * Gr_1M * Aux2 + Aux5')
          call HI_zgemm ('N', 'N', NL, norbDyn(idx), NR, (0.d0,0.5d0),  &
                         Gr_1M, NL, Aux2, NR, (1.d0,0.d0), Aux5, NL)

!         ('Aux6 = Gamma_L * Aux5')
          call HI_zhemm ('L', 'L', NL, norbDyn(idx), (1.d0,0.d0),       &
                         Gamma_L, NL, Aux5, NL, (0.d0,0.d0), Aux6, NL)

!         ('Aux4 = Gr_1n^dagger * Aux6')
          call HI_zgemm ('C', 'N', norbDyn(idx), norbDyn(idx), NL,      &
                         (1.d0,0.d0), auxGr_1n, NL, Aux6, NL,           &
                         (0.d0,0.d0), Aux4, norbDyn(idx))

!         -- 4rd PART: 'G^dagger*Gamma_L*G*(-i/2 * H.c.)' --

!         ('Aux6 = Gamma_L * Gr_1n')
          call HI_zhemm ('L', 'L', NL, norbDyn(idx), (1.d0,0.d0),       &
                         Gamma_L, NL, auxGr_1n, NL, (0.d0,0.d0),        &
                         Aux6, NL)

!         ('Aux1 = Gr_1M^dagger * Aux6')
          call HI_zgemm ('C', 'N', NR, norbDyn(idx), NL, (1.d0,0.d0),   &
                         Gr_1M, NL, Aux6, NL, (0.d0,0.d0), Aux1, NR)

!         ('Aux7 = Aux1 * (-i/2 * H.c.)') (where 'H.c. = Aux2^dagger')
          call HI_zgemm ('N', 'C', NR, NR, norbDyn(idx), (0.d0,-0.5d0), &
                         Aux1, NR, Aux2, NR, (0.d0,0.d0), Aux7, NR)

!         Compute the trace.
          do i = 1,norbDyn(idx)
             Isymm = Isymm + DREAL(Aux4(i,i))
          enddo
          do i = 1,NR
             Isymm = Isymm + DREAL(Aux7(i,i))
          enddo
          
#ifdef DEBUG
!         [test] Compute the matrices multiplication with full matrices.
          call testInelSymm (Ei, ispin, NL, Gamma_L,                    &
                             NR, Gamma_R, j, idx, norbDyn(idx),         &
                             Meph(idx)%M(:,:,ispin,w), Aux4, Aux7)
#endif

!         Compute symmetric pre-factor.
          foo = 2.d0 * Vbias * phOccup(j)%P(ienergy,ispin,iv,w)
          foo = foo + (freq(idx)%F(w) - Vbias)                          &
                * BoseEinstein (freq(idx)%F(w) - Vbias, temp)
          foo = foo - (freq(idx)%F(w) + Vbias)                          &
                * BoseEinstein (freq(idx)%F(w) + Vbias, temp)
          Isymm = eoverh * foo * Isymm

       enddo ! do w = 1,nModes(idx)

!      Free memory.
       deallocate (GrCJG)
       deallocate (A)
       deallocate (Aux1)
       deallocate (Aux2)
       deallocate (Aux3)
       deallocate (Aux4)
       deallocate (Aux5)
       deallocate (Aux6)
       deallocate (Aux7)
       deallocate (auxGr_1n)
       deallocate (auxGr_Mn)

    enddo ! do j = 1,nunitseph

!   Close Green's functions files.
    call IOcloseStream (Gr_nn_disk%fname, Gr_nn_disk%lun)
    call IOcloseStream (Gr_1n_disk%fname, Gr_1n_disk%lun)
    call IOcloseStream (Gr_Mn_disk%fname, Gr_Mn_disk%lun)


  end subroutine inelSymmDisk


!  *******************************************************************  !
!                               asymmPre                                !
!  *******************************************************************  !
!  Description: compute pre-factors from the asymmetric part of         !
!  inelastic current.                                                   !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer nAsymmPts           : Number of energy grid points           !
!                                for asymmetric term integral           !
!  real*8 temp                 : Electronic temperature                 !
!  ****************************** INPUT ******************************  !
!  real*8 Ei                           : Energy grid point              !
!  real*8 freq                         : Vibrational mode frequency     !
!  real*8 Vbias                        : Bias potential                 !
!  *******************************************************************  !
  real(8) function asymmPre (Ei, freq, VBias)

!
!   Modules
!
    use idsrdr_options,  only: nAsymmPts, temp
    use idsrdr_recipes,  only: RECPSsimpson
    use idsrdr_distrib,  only: FermiDirac
    use idsrdr_hilbert,  only: hilbert
    use idsrdr_leads,    only: EfLead

!   Input variables.
    real(8), intent(in) :: Ei, freq, Vbias

!   Local variables.
    integer :: k
    real(8) :: enI, enF
    real(8), allocatable, dimension (:) :: En ! energy grid points
    real(8), allocatable, dimension (:) :: We ! energy grid weights
    complex(8), allocatable, dimension (:) :: aux

!   Set lower and upper limit of energy integration.
    if (Vbias >= 0.d0) then
       enI = Ei - 2.d0*freq - kbTol*temp - Vbias
       enF = Ei + 2.d0*freq + kbTol*temp
    else
       enI = Ei - 2.d0*freq - kbTol*temp
       enF = Ei + 2.d0*freq + kbTol*temp - Vbias
    endif

!   Allocate the energy grid points and weights arrays.
    allocate (En(nAsymmPts), We(nAsymmPts))

!   Compute energy points and weights on an equidistant grid.
    call RECPSsimpson (enI, enF, nAsymmPts, En, We)

!   Allocate auxiliary array.
    allocate (aux(2*nAsymmPts))

!   ('aux = f(E-w) - f(E+w)')
    aux = 0.d0
    do k = 1,nAsymmPts
       aux(k) = FermiDirac (En(k)+freq, EfLead, temp)                   &
                - FermiDirac (En(k)-freq, EfLead, temp)
    enddo

!   Compute the pre-factor (with Hilbert transform).
    call hilbert (nAsymmPts, aux)

    asymmPre = 0.d0
    do k = 1,nAsymmPts
       asymmPre = asymmPre + We(k) * aux(k)                             &
                  * (FermiDirac (En(k), EfLead, temp)                   &
                  - FermiDirac (En(k)-Vbias, EfLead, temp))
    enddo
    asymmPre = eoverh * asymmPre / 2.d0

!   Free memory.
    deallocate (En, We)
    deallocate (aux)


  end function asymmPre


!  *******************************************************************  !
!                               inelAsymm                               !
!  *******************************************************************  !
!  Description: compute the asymmetric part of inelastic transmission.  !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer nModes(neph)        : Number of vibrational modes            !
!  integer norbDyn(neph)       : Number of orbitals from dynamic atoms  !
!  TYPE(ephCplng) Meph(neph)%M(norbDyn,norbDyn,nspin,nModes) :          !
!                                   [complex*8] E-ph coupling matrices  !
!  TYPE(ephFreq) freq(neph)%F(nModes) : [real*8] Vibrational mode's     !
!                                       frequencies                     !
!  integer ephIdx(ntypeunits+2)       : Unit index (those with e-ph)    !
!  integer nunitseph                  : Number of units with eph        !
!  integer eph_type(nunitseph)        : Units types with eph            !
!  TYPE(green) Gr_Mn(nunitseph)%G(NR,unitdimensions) : [complex]        !
!                                                      G^r_{M,n}        !
!  TYPE(green) Gr_1n(nunitseph)%G(NL,unitdimensions) : [complex]        !
!                                                      G^r_{1,n}        !
!  complex(8) Gr_1M(NL,NR)                           : G^r_{1,M}        !
!  ****************************** INPUT ******************************  !
!  integer ispin                       : Spin component index           !
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Gamma_L(NL,NL)           : Left-lead coupling matrix      !
!  complex(8) Gamma_R(NR,NR)           : Right-lead coupling matrix     !
!  real*8 Vbias                        : Bias potential                 !
!  real*8 Ei                           : Energy grid point              !
!  ***************************** OUTPUT ******************************  !
!  real*8 Iasymm               : Asymmetric part of inelastic current   !
!  *******************************************************************  !
  subroutine inelAsymm (Iasymm, ispin, NL, Gamma_L,                     &
                        NR, Gamma_R, Vbias, Ei)

!
!   Modules
!
    use idsrdr_ephcoupl, only: nModes, norbDyn, Meph, freq, ephIdx
    use idsrdr_units,    only: nunitseph, eph_type
    use idsrdr_green,    only: Gr_Mn, Gr_1n, Gr_1M

!   Input variables.
    integer, intent(in) :: ispin, NL, NR
    real(8), intent(in) :: Ei, Vbias
    real(8), intent(out) :: Iasymm
    complex(8), dimension (NL,NL), intent(in) :: Gamma_L
    complex(8), dimension (NR,NR), intent(in) :: Gamma_R

!   Local variables.
    integer :: j, w, i, idx
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8), allocatable, dimension(:,:) :: Aux1, Aux2, Aux3,        &
                                               Aux4, Aux5, Aux6,        &
                                               Gr_MnCJG, Gr_1nCJG
    external :: HI_zsymm, HI_zhemm, HI_zgemm

!   Initialize variable.
    Iasymm = 0.d0

    do j = 1,nunitseph ! over unit with e-ph

       idx = ephIdx(eph_type(j))

!      Allocate auxiliary matrices.
       allocate (Gr_MnCJG(NR,norbDyn(idx)))
       allocate (Gr_1nCJG(NL,norbDyn(idx)))
       allocate (Aux1(NR,norbDyn(idx)))
       allocate (Aux2(NR,norbDyn(idx)))
       allocate (Aux3(norbDyn(idx),norbDyn(idx)))
       allocate (Aux4(NL,norbDyn(idx)))
       allocate (Aux5(NL,norbDyn(idx)))
       allocate (Aux6(NR,NR))

!      Copy the complex conjugate of 'Gr_Mn' and 'Gr_1n'.
       Gr_MnCJG = DCONJG(Gr_Mn(j)%G)
       Gr_1nCJG = DCONJG(Gr_1n(j)%G)

       do w = 1,nModes(idx) ! over phonon modes

!         -- 1st PART: 'G*Gamma_R*G^dagger*Meph' --

!         ('Aux1 = Gamma_R * Gr_Mn^*')
          call HI_zhemm ('L', 'L', NR, norbDyn(idx), (1.d0,0.d0),       &
                         Gamma_R, NR, Gr_MnCJG, NR, (0.d0,0.d0),        &
                         Aux1, NR)

!         ('Aux2 = Aux1 * Meph')
          call HI_zsymm ('R', 'L', NR, norbDyn(idx), (1.d0,0.d0),       &
                         Meph(idx)%M(:,:,ispin,w), norbDyn(idx),        &
                         Aux1, NR, (0.d0,0.d0), Aux2, NR)

!         ('Aux3 = Gr_Mn^T * Aux2')
          call HI_zgemm ('T', 'N', norbDyn(idx), norbDyn(idx), NR,      &
                         (1.d0,0.d0), Gr_Mn(j)%G, NR, Aux2, NR,         &
                         (0.d0,0.d0), Aux3, norbDyn(idx))

!         -- 2nd PART: '-G*Gamma_L*G^dagger*Meph + 1st PART' --

!         ('Aux4 = Gr_1n^* * Meph')
          call HI_zsymm ('R', 'L', NL, norbDyn(idx), (1.d0,0.d0),       &
                         Meph(idx)%M(:,:,ispin,w), norbDyn(idx),        &
                         Gr_1nCJG, NL, (0.d0,0.d0), Aux4, NL)

!         ('Aux5 = Gamma_L * Aux4')
          call HI_zhemm ('L', 'L', NL, norbDyn(idx), (1.d0,0.d0),       &
                         Gamma_L, NL, Aux4, NL, (0.d0,0.d0), Aux5, NL)

!         ('Aux3 = - Gr_1n^T * Aux5 + Aux3')
          call HI_zgemm ('T', 'N', norbDyn(idx), norbDyn(idx), NL,      &
                         (-1.d0,0.d0), Gr_1n(j)%G, NL, Aux5, NL,        &
                         (1.d0,0.d0), Aux3, norbDyn(idx))

!         -- 3rd PART:
!                'G^dagger*Gamma_L*G*Gamma_R*G^dagger*Meph*(2nd Part)' --

!         ('Aux1 = Aux2 * Aux3') (where 'Aux2= Gamma_R * Gr_Mn^* * Meph')
          call HI_zgemm ('N', 'N', NR, norbDyn(idx), norbDyn(idx),      &
                         (1.d0,0.d0), Aux2, NR, Aux3,                   &
                         norbDyn(idx), (0.d0,0.d0), Aux1, NR)

!         ('Aux4 = Gr_1M * Aux1')
          call HI_zgemm ('N', 'N', NL, norbDyn(idx), NR, (1.d0,0.d0),   &
                         Gr_1M, NL, Aux1, NR, (0.d0,0.d0), Aux4, NL)

!         ('Aux5 = Gamma_L*Aux4')
          call HI_zhemm ('L', 'L', NL, norbDyn(idx), (1.d0,0.d0),       &
                         Gamma_L, NL, Aux4, NL, (0.d0,0.d0), Aux5, NL)

!         ('Aux3 = Gr_1n^dagger * Aux5')
          call HI_zgemm ('C', 'N', norbDyn(idx), norbDyn(idx), NL,      &
                         (1.d0,0.d0), Gr_1n(j)%G, NL, Aux5, NL,         &
                         (0.d0,0.d0), Aux3, norbDyn(idx))

!         -- 4rd PART: 'G^dagger*Gamma_L*G*(H.c.)' --

!         ('Aux4 = Gamma_L * Gr_1n')
          call HI_zhemm ('L', 'L', NL, norbDyn(idx), (1.d0,0.d0),       &
                         Gamma_L, NL, Gr_1n(j)%G, NL, (0.d0,0.d0),      &
                         Aux4, NL)

!         ('Aux2 = Gr_1M^dagger * Aux4')
          call HI_zgemm ('C', 'N', NR, norbDyn(idx), NL, (1.d0,0.d0),   &
                         Gr_1M, NL, Aux4, NL, (0.d0,0.d0), Aux2, NR)

!         ('Aux6 = Aux2 * (H.c.)') ('H.c. = Aux1^dagger')
          call HI_zgemm ('N', 'C', NR, NR, norbDyn(idx), (1.d0,0.0d0),  &
                         Aux2, NR, Aux1, NR, (0.d0,0.d0), Aux6, NR)

!         Compute the trace.
          do i = 1,norbDyn(idx)
             Iasymm = Iasymm + DREAL(Aux3(i,i))
          enddo
          do i = 1,NR
             Iasymm = Iasymm + DREAL(Aux6(i,i))
          enddo

#ifdef DEBUG
!         [test] Compute the matrices multiplication with full matrices.
          call testInelAsymm (Ei, ispin, NL, Gamma_L,                   &
                              NR, Gamma_R, j, idx, norbDyn(idx),        &
                              Meph(idx)%M(:,:,ispin,w), Aux3, Aux6)
#endif

!         Compute asymmetric pre-factor.
          Iasymm = asymmPre (Ei, freq(idx)%F(w), Vbias) * Iasymm
          
       enddo ! do w = 1,nModes(idx)

!      Free memory.
       deallocate (Gr_MnCJG)
       deallocate (Gr_1nCJG)
       deallocate (Aux1)
       deallocate (Aux2)
       deallocate (Aux3)
       deallocate (Aux4)
       deallocate (Aux5)
       deallocate (Aux6)

    enddo ! do j = 1,nunitseph


  end subroutine inelAsymm


!  *******************************************************************  !
!                             inelAsymmDisk                             !
!  *******************************************************************  !
!  Description: compute the asymmetric part of inelastic transmission.  !
!  (Green's functions are read from disk).                              !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer nModes(neph)        : Number of vibrational modes            !
!  integer norbDyn(neph)       : Number of orbitals from dynamic atoms  !
!  TYPE(ephCplng) Meph(neph)%M(norbDyn,norbDyn,nspin,nModes) :          !
!                                   [complex*8] E-ph coupling matrices  !
!  TYPE(ephFreq) freq(neph)%F(nModes) : [real*8] Vibrational mode's     !
!                                       frequencies                     !
!  integer ephIdx(ntypeunits+2)       : Unit index (those with e-ph)    !
!  integer nunitseph                  : Number of units with eph        !
!  integer eph_type(nunitseph)        : Units types with eph            !
!  TYPE(green) Gr_Mn(nunitseph)%G(NR,unitdimensions) : [complex]        !
!                                                      G^r_{M,n}        !
!  TYPE(green) Gr_1n(nunitseph)%G(NL,unitdimensions) : [complex]        !
!                                                      G^r_{1,n}        !
!  complex(8) Gr_1M(NL,NR)                           : G^r_{1,M}        !
!  ****************************** INPUT ******************************  !
!  integer ispin                       : Spin component index           !
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Gamma_L(NL,NL)           : Left-lead coupling matrix      !
!  complex(8) Gamma_R(NR,NR)           : Right-lead coupling matrix     !
!  real*8 Vbias                        : Bias potential                 !
!  real*8 Ei                           : Energy grid point              !
!  ***************************** OUTPUT ******************************  !
!  real*8 Iasymm               : Asymmetric part of inelastic current   !
!  *******************************************************************  !
  subroutine inelAsymmDisk (Iasymm, ispin, NL, Gamma_L,                 &
                            NR, Gamma_R, Vbias, Ei)

!
!   Modules
!
    use idsrdr_ephcoupl, only: nModes, norbDyn, Meph, freq, ephIdx
    use idsrdr_units,    only: nunitseph, eph_type
    use idsrdr_green,    only: Gr_Mn_disk, Gr_1n_disk, Gr_1M, greenload
    use idsrdr_io,       only: IOopenStream, IOcloseStream

!   Input variables.
    integer, intent(in) :: ispin, NL, NR
    real(8), intent(in) :: Ei, Vbias
    real(8), intent(out) :: Iasymm
    complex(8), dimension (NL,NL), intent(in) :: Gamma_L
    complex(8), dimension (NR,NR), intent(in) :: Gamma_R

!   Local variables.
    integer :: j, w, i, idx
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8), allocatable, dimension(:,:) :: Aux1, Aux2, Aux3,        &
                                               Aux4, Aux5, Aux6,        &
                                               Gr_MnCJG, Gr_1nCJG,      &
                                               auxGr_1n, auxGr_Mn
    external :: HI_zsymm, HI_zhemm, HI_zgemm

!   Initialize variable.
    Iasymm = 0.d0

!   Open Green's functions files.
    call IOopenStream (Gr_1n_disk%fname, Gr_1n_disk%lun)
    call IOopenStream (Gr_Mn_disk%fname, Gr_Mn_disk%lun)

    do j = 1,nunitseph ! over unit with e-ph

       idx = ephIdx(eph_type(j))

!      Allocate auxiliary matrices.
       allocate (Gr_MnCJG(NR,norbDyn(idx)))
       allocate (Gr_1nCJG(NL,norbDyn(idx)))
       allocate (Aux1(NR,norbDyn(idx)))
       allocate (Aux2(NR,norbDyn(idx)))
       allocate (Aux3(norbDyn(idx),norbDyn(idx)))
       allocate (Aux4(NL,norbDyn(idx)))
       allocate (Aux5(NL,norbDyn(idx)))
       allocate (Aux6(NR,NR))

!      Copy 'Gr_1n' from file to auxiliary matrix.
       allocate (auxGr_1n(NL,norbDyn(idx)))
       call greenload (Gr_1n_disk%lun, NL, norbDyn(idx), 1, NL,         &
                       1, norbDyn(idx), auxGr_1n, NL, norbDyn(idx))

!      Copy 'Gr_Mn' from file to auxiliary matrix.
       allocate (auxGr_Mn(NR,norbDyn(idx)))
       call greenload (Gr_Mn_disk%lun, NR, norbDyn(idx), 1, NR,         &
                       1, norbDyn(idx), auxGr_Mn, NR, norbDyn(idx))

!      Copy the complex conjugate of 'Gr_Mn' and 'Gr_1n'.
       Gr_MnCJG = DCONJG(auxGr_Mn)
       Gr_1nCJG = DCONJG(auxGr_1n)

       do w = 1,nModes(idx) ! over phonon modes

!         -- 1st PART: 'G*Gamma_R*G^dagger*Meph' --

!         ('Aux1 = Gamma_R * Gr_Mn^*')
          call HI_zhemm ('L', 'L', NR, norbDyn(idx), (1.d0,0.d0),       &
                         Gamma_R, NR, Gr_MnCJG, NR, (0.d0,0.d0),        &
                         Aux1, NR)

!         ('Aux2 = Aux1 * Meph')
          call HI_zsymm ('R', 'L', NR, norbDyn(idx), (1.d0,0.d0),       &
                         Meph(idx)%M(:,:,ispin,w), norbDyn(idx),        &
                         Aux1, NR, (0.d0,0.d0), Aux2, NR)

!         ('Aux3 = Gr_Mn^T * Aux2')
          call HI_zgemm ('T', 'N', norbDyn(idx), norbDyn(idx), NR,      &
                         (1.d0,0.d0), auxGr_Mn, NR, Aux2, NR,           &
                         (0.d0,0.d0), Aux3, norbDyn(idx))

!         -- 2nd PART: '-G*Gamma_L*G^dagger*Meph + 1st PART' --

!         ('Aux4 = Gr_1n^* * Meph')
          call HI_zsymm ('R', 'L', NL, norbDyn(idx), (1.d0,0.d0),       &
                         Meph(idx)%M(:,:,ispin,w), norbDyn(idx),        &
                         Gr_1nCJG, NL, (0.d0,0.d0), Aux4, NL)

!         ('Aux5 = Gamma_L * Aux4')
          call HI_zhemm ('L', 'L', NL, norbDyn(idx), (1.d0,0.d0),       &
                         Gamma_L, NL, Aux4, NL, (0.d0,0.d0), Aux5, NL)

!         ('Aux3 = - Gr_1n^T * Aux5 + Aux3')
          call HI_zgemm ('T', 'N', norbDyn(idx), norbDyn(idx), NL,      &
                         (-1.d0,0.d0), auxGr_1n, NL, Aux5, NL,          &
                         (1.d0,0.d0), Aux3, norbDyn(idx))

!         -- 3rd PART:
!                'G^dagger*Gamma_L*G*Gamma_R*G^dagger*Meph*(2nd Part)' --

!         ('Aux1 = Aux2 * Aux3') (where 'Aux2= Gamma_R * Gr_Mn^* * Meph')
          call HI_zgemm ('N', 'N', NR, norbDyn(idx), norbDyn(idx),      &
                         (1.d0,0.d0), Aux2, NR, Aux3,                   &
                         norbDyn(idx), (0.d0,0.d0), Aux1, NR)

!         ('Aux4 = Gr_1M * Aux1')
          call HI_zgemm ('N', 'N', NL, norbDyn(idx), NR, (1.d0,0.d0),   &
                         Gr_1M, NL, Aux1, NR, (0.d0,0.d0), Aux4, NL)

!         ('Aux5 = Gamma_L*Aux4')
          call HI_zhemm ('L', 'L', NL, norbDyn(idx), (1.d0,0.d0),       &
                         Gamma_L, NL, Aux4, NL, (0.d0,0.d0), Aux5, NL)

!         ('Aux3 = Gr_1n^dagger * Aux5')
          call HI_zgemm ('C', 'N', norbDyn(idx), norbDyn(idx), NL,      &
                         (1.d0,0.d0), auxGr_1n, NL, Aux5, NL,           &
                         (0.d0,0.d0), Aux3, norbDyn(idx))

!         -- 4rd PART: 'G^dagger*Gamma_L*G*(H.c.)' --

!         ('Aux4 = Gamma_L * Gr_1n')
          call HI_zhemm ('L', 'L', NL, norbDyn(idx), (1.d0,0.d0),       &
                         Gamma_L, NL, auxGr_1n, NL, (0.d0,0.d0),        &
                         Aux4, NL)

!         ('Aux2 = Gr_1M^dagger * Aux4')
          call HI_zgemm ('C', 'N', NR, norbDyn(idx), NL, (1.d0,0.d0),   &
                         Gr_1M, NL, Aux4, NL, (0.d0,0.d0), Aux2, NR)

!         ('Aux6 = Aux2 * (H.c.)') ('H.c. = Aux1^dagger')
          call HI_zgemm ('N', 'C', NR, NR, norbDyn(idx), (1.d0,0.0d0),  &
                         Aux2, NR, Aux1, NR, (0.d0,0.d0), Aux6, NR)

!         Compute the trace.
          do i = 1,norbDyn(idx)
             Iasymm = Iasymm + DREAL(Aux3(i,i))
          enddo
          do i = 1,NR
             Iasymm = Iasymm + DREAL(Aux6(i,i))
          enddo

#ifdef DEBUG
!         [test] Compute the matrices multiplication with full matrices.
          call testInelAsymm (Ei, ispin, NL, Gamma_L,                   &
                              NR, Gamma_R, j, idx, norbDyn(idx),        &
                              Meph(idx)%M(:,:,ispin,w), Aux3, Aux6)
#endif

!         Compute asymmetric pre-factor.
          Iasymm = asymmPre (Ei, freq(idx)%F(w), Vbias) * Iasymm
          
       enddo ! do w = 1,nModes(idx)

!      Free memory.
       deallocate (Gr_MnCJG)
       deallocate (Gr_1nCJG)
       deallocate (Aux1)
       deallocate (Aux2)
       deallocate (Aux3)
       deallocate (Aux4)
       deallocate (Aux5)
       deallocate (Aux6)
       deallocate (auxGr_1n)
       deallocate (auxGr_Mn)

    enddo ! do j = 1,nunitseph

!   Close Green's functions files.
    call IOcloseStream (Gr_1n_disk%fname, Gr_1n_disk%lun)
    call IOcloseStream (Gr_Mn_disk%fname, Gr_Mn_disk%lun)


  end subroutine inelAsymmDisk


!  *******************************************************************  !
!                             testInelSymm                              !
!  *******************************************************************  !
!  Description: build the full hamiltonian of the scattering region     !
!  and compute the symmetric part of inelastic transmission in order    !
!  compare with the transmission obtained with 'inelSymm' subroutine.   !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer ntypeunits                   : Number of unit types          !
!  integer nunits                       : Total number of units         !
!  integer unit_type(nunits+2)          : Units types                   !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  real*8 unitshift(ntypeunits+2)       : Units shift                   !
!  real*8 S1unit(unitdimensions(ntypeunits),                            !
!                unitdimensions(ntypeunits)) : Unit coupling overlap    !
!  real*8 H1unit(unitdimensions(ntypeunits),                            !
!                unitdimensions(ntypeunits),nspin) : Unit coupling      !
!                                                    Hamiltonian        !
!  TYPE(unitS) Sunits(ntypeunits+2)%S(unitdimensions,unitdimensions) :  !
!                                         [real*8] Units overlap        !
!  TYPE(unitH)                                                          !
!        Hunits(ntypeunits+2)%H(unitdimensions,unitdimensions,nspin) :  !
!                                         [real*8] Units hamiltonian    !
!  complex(8) Sigma_L(NL,NL)           : Left-lead self-energy          !
!  complex(8) Sigma_R(NR,NR)           : Right-lead self-energy         !
!  integer ephIdx(ntypeunits+2)        : Unit index (those with e-ph)   !
!  integer idxF(neph)                  : First dynamic atom orbital     !
!  integer idxL(neph)                  : Last dynamic atom orbital      !
!  ****************************** INPUT ******************************  !
!  real*8 Ei                           : Energy grid point              !
!  integer ispin                       : Spin component index           !
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Gamma_L(NL,NL)           : Left-lead coupling matrix      !
!  complex(8) Gamma_R(NR,NR)           : Right-lead coupling matrix     !
!  integer ephunit                     : E-ph unit index                !
!  integer ephtype                     : E-ph unit type index           !
!  integer norbDyn                     : # of dynamic atoms orbitals    !
!  complex(8) Meph(norbDyn,norbDyn)    : E-ph coupling matrix           !
!  complex(8) Symm(norbDyn,norbDyn)    : Calculated transmission        !
!  complex(8) Hc(NR,NR)                : Hermitian conjugated part      !
!  *******************************************************************  !
  subroutine testInelSymm (Ei, ispin, NL, Gamma_L, NR, Gamma_R,         &
                           ephunit, ephtype, norbDyn, Meph, Symm, Hc)

!
!   Modules
!
    use idsrdr_options,  only: ntypeunits, nunits
    use idsrdr_units,    only: unit_type, unitdimensions, unitshift,    &
                               S1unit, H1unit, Sunits, Hunits
    use idsrdr_leads,    only: Sigma_L, Sigma_R
    use idsrdr_ephcoupl, only: ephIdx, idxF, idxL

!   Input variables.
    integer, intent(in) :: ispin, NL, NR, ephunit, ephtype, norbDyn
    real(8), intent(in) :: Ei
    complex(8), dimension (NL,NL), intent(in) :: Gamma_L
    complex(8), dimension (NR,NR), intent(in) :: Gamma_R
    complex(8), dimension (norbDyn,norbDyn), intent(in) :: Meph
    complex(8), dimension (norbDyn,norbDyn), intent(in) :: Symm
    complex(8), dimension (NR,NR), intent(in) :: Hc

!   Local variables.
    integer :: i, j, k, utype, dim, dimTot, dimCpl, idxAnt, ueph
    real(8), allocatable, dimension (:,:) :: STot
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8) :: foo
    complex(8), allocatable, dimension (:,:) :: HTot, GrTot, MephTot,   &
                                                Gamma_LTot, Gamma_RTot, &
                                                Aux1, Aux2, Aux3
    external :: HI_zsymm, HI_zhemm, HI_zgemm, HI_zgeInvert

!   Get the total dimension.
    dimTot = unitdimensions(ntypeunits+1) + unitdimensions(ntypeunits+2)
    do I = 2,nunits+1
       utype = unit_type(I) ! current unit type
       dimTot = dimTot + unitdimensions(utype)
    enddo

!   Allocate and initialize matrices.
    allocate (STot(dimTot,dimTot))
    allocate (HTot(dimTot,dimTot))
    allocate (GrTot(dimTot,dimTot))
    STot = 0.d0
    HTot = 0.d0

!   Build total hamiltonian and overlap matrices.
    dimCpl = unitdimensions(ntypeunits) ! coupling unit dimensions
    utype = ntypeunits+1 ! current unit type (first unit)
    dim = unitdimensions(utype) ! current unit dimensions

    STot(1:dim,1:dim) = (Ei-unitshift(utype))*Sunits(utype)%S(:,:)
    HTot(1:dim,1:dim) = Hunits(utype)%H(:,:,ispin)
    STot(dim-dimCpl+1:dim,dim+1:dim+dimCpl) =                           &
         (Ei-unitshift(ntypeunits))*S1unit
    HTot(dim-dimCpl+1:dim,dim+1:dim+dimCpl) = H1unit(:,:,ispin)
    STot(dim+1:dim+dimCpl,dim-dimCpl+1:dim) =                           &
         (Ei-unitshift(ntypeunits))*TRANSPOSE(S1unit)
    HTot(dim+1:dim+dimCpl,dim-dimCpl+1:dim) =                           &
         TRANSPOSE(H1unit(:,:,ispin))

    idxAnt = dim + 1
    do k = 2,nunits+1
       utype = unit_type(k) ! current unit type
       dim = unitdimensions(utype) ! current unit dimensions

       STot(idxAnt:idxAnt+dim-1,idxAnt:idxAnt+dim-1) =                  &
            (Ei-unitshift(utype))*Sunits(utype)%S(:,:)
       HTot(idxAnt:idxAnt+dim-1,idxAnt:idxAnt+dim-1) =                  &
            Hunits(utype)%H(:,:,ispin)
       STot(idxAnt+dim-dimCpl:idxAnt+dim-1,                             &
            idxAnt+dim:idxAnt+dim+dimCpl-1) =                           &
            (Ei-unitshift(ntypeunits))*S1unit
       HTot(idxAnt+dim-dimCpl:idxAnt+dim-1,                             &
            idxAnt+dim:idxAnt+dim+dimCpl-1) = H1unit(:,:,ispin)
       STot(idxAnt+dim:idxAnt+dim+dimCpl-1,                             &
            idxAnt+dim-dimCpl:idxAnt+dim-1) =                           &
            (Ei-unitshift(ntypeunits))*TRANSPOSE(S1unit)
       HTot(idxAnt+dim:idxAnt+dim+dimCpl-1,                             &
            idxAnt+dim-dimCpl:idxAnt+dim-1) =                           &
            TRANSPOSE(H1unit(:,:,ispin))

       idxAnt = idxAnt + dim

    enddo

    utype = ntypeunits + 2 ! current unit type (last unit)
    dim = unitdimensions(utype) ! current unit dimensions

    STot(idxAnt:idxAnt+dim-1,idxAnt:idxAnt+dim-1) =                     &
         (Ei-unitshift(utype))*Sunits(utype)%S(:,:)
    HTot(idxAnt:idxAnt+dim-1,idxAnt:idxAnt+dim-1) =                     &
         Hunits(utype)%H(:,:,ispin)

!   Green's function.
    GrTot = STot - HTot
    GrTot(1:NL,1:NL) = GrTot(1:NL,1:NL) - Sigma_L
    GrTot(dimTot-NR+1:dimTot,dimTot-NR+1:dimTot) =                      &
         GrTot(dimTot-NR+1:dimTot,dimTot-NR+1:dimTot) - Sigma_R
    call HI_zgeInvert (Grtot, dimTot)

!   Free memory.
    deallocate (HTot)
    deallocate (STot)

!   Allocate full matrices.
    allocate (Gamma_LTot(dimTot,dimTot))
    allocate (Gamma_RTot(dimTot,dimTot))
    allocate (MephTot(dimTot,dimTot))

!   Assign the full lead's coupling matrices.
    Gamma_LTot = (0.d0,0.d0)
    Gamma_RTot = (0.d0,0.d0)
    do j = 1,NL
       do i = 1,NL
          Gamma_LTot(i,j) = Gamma_L(i,j)
       enddo
    enddo
    do j = 1,NR
       do i = 1,NR
          Gamma_RTot(dimTot-NR+i,dimTot-NR+j) = Gamma_R(i,j)
       enddo
    enddo

!   Assign the full e-ph coupling matrix.
    MephTot = (0.d0,0.d0)
    if (ephtype == ephIdx(ntypeunits+1)) then
       do j = idxF(ephtype),idxL(ephtype)
          do i = idxF(ephtype),idxL(ephtype)
             MephTot(i,j) = Meph(i-idxF(ephtype)+1,j-idxF(ephtype)+1)
          enddo
       enddo
    else if (ephtype == ephIdx(ntypeunits+2)) then
       dim = unitdimensions(ntypeunits+2)
       do j = idxF(ephtype),idxL(ephtype)
          do i = idxF(ephtype),idxL(ephtype)
             MephTot(dimTot-dim+i,dimTot-dim+j) =                       &
                  Meph(i-idxF(ephtype)+1,j-idxF(ephtype)+1)
          enddo
       enddo
    else
       dim = unitdimensions(ntypeunits+1)
       if (ephIdx(ntypeunits+1) /= 0) then ! first unit has eph?
          ueph = 2 ! eph units indexing
       else
          ueph = 1 ! eph units indexing
       endif
       do k = 2,nunits+1
          utype = unit_type(k) ! current unit type
          if (ephIdx(utype) /= 0) then ! unit with eph?
             if (ephunit == ueph) then ! the unit we are looking?
                do j = idxF(ephtype),idxL(ephtype)
                   do i = idxF(ephtype),idxL(ephtype)
                      MephTot(dim+i,dim+j) =                            &
                           Meph(i-idxF(ephtype)+1,j-idxF(ephtype)+1)
                   enddo
                enddo
                EXIT
             else
                ueph = ueph + 1
             endif
          endif
          dim = dim + unitdimensions(utype)
       enddo
    endif

!   Allocate auxiliaries matrices.
    allocate (Aux1(dimTot,dimTot))
    allocate (Aux2(dimTot,dimTot))
    allocate (Aux3(dimTot,dimTot))

!   -- 1st PART: 'i/2*G*(Gamma_R*G^dagger*Meph*A*Meph - H.c.' --

!   ('Aux1 = Gamma_RTot * GrTot^dagger')
    call HI_zgemm ('N', 'C', dimTot, dimTot, dimTot, (1.d0,0.d0),       &
                   Gamma_RTot, dimTot, GrTot, dimTot,                   &
                   (0.d0,0.d0), Aux1, dimTot)

!   ('Aux2 = Aux1 * MephTot')
    call HI_zsymm ('R', 'L', dimTot, dimTot, (1.d0,0.d0), MephTot,      &
                   dimTot, Aux1, dimTot, (0.d0,0.d0), Aux2, dimTot)

!   Spectral matrix (obs.: 'GrTot' is symmetric).
    Aux1 = zi * (GrTot - DCONJG(GrTot))

!   ('Aux3 = Aux2 * Aux1')
    call HI_zhemm ('R', 'L', dimTot, dimTot, (1.d0,0.d0), Aux1,         &
                   dimTot, Aux2, dimTot, (0.d0,0.d0), Aux3, dimTot)

!   ('Aux1 = Aux3 * MephTot')
    call HI_zsymm ('R', 'L', dimTot, dimTot, (1.d0,0.d0), MephTot,      &
                   dimTot, Aux3, dimTot, (0.d0,0.d0), Aux1, dimTot)

!   ('Aux3 = i/2 * GrTot * Aux1')
    call HI_zgemm ('N', 'N', dimTot, dimTot, dimTot, (0.d0,0.5d0),      &
                   GrTot, dimTot, Aux1, dimTot, (0.d0,0.d0),            &
                   Aux3, dimTot)

!   ('Aux3 = -i/2 * GrTot * Aux1^dagger + Aux3')
    call HI_zgemm ('N', 'C', dimTot, dimTot, dimTot, (0.d0,-0.5d0),     &
                   GrTot, dimTot, Aux1, dimTot, (1.d0,0.d0),            &
                   Aux3, dimTot)

!   -- 2nd PART: 'G*Meph*G*Gamma_R*G^dagger*Meph + 1st PART' --

!   ('Aux1 = GrTot * Aux2')
!                    (where 'Aux2 = Gamma_RTot * GrTot^dagger * MephTot')
    call HI_zgemm ('N', 'N', dimTot, dimTot, dimTot, (1.d0,0.d0),       &
                   GrTot, dimTot, Aux2, dimTot, (0.d0,0.d0),            &
                   Aux1, dimTot)

!   ('Aux2 = MephTot * Aux1')
    call HI_zsymm ('L', 'L', dimTot, dimTot, (1.d0,0.d0), MephTot,      &
                   dimTot, Aux1, dimTot, (0.d0,0.d0), Aux2, dimTot)

!   ('Aux3 = GrTot * Aux2 + Aux3')
    call HI_zgemm ('N', 'N', dimTot, dimTot, dimTot, (1.d0,0.d0),       &
                   GrTot, dimTot, Aux2, dimTot, (1.d0,0.d0),            &
                   Aux3, dimTot)

!   -- 3rd PART: 'G^dagger*Gamma_L*(2nd PART)' --

!   ('Aux2 = Gamma_LTot * Aux3')
    call HI_zhemm ('L', 'L', dimTot, dimTot, (1.d0,0.d0), Gamma_LTot,   &
                   dimTot, Aux3, dimTot, (0.d0,0.d0), Aux2, dimTot)

!   ('Aux1 = GrTot^dagger * Aux2')
    call HI_zgemm ('C', 'N', dimTot, dimTot, dimTot, (1.d0,0.d0),       &
                   GrTot, dimTot, Aux2, dimTot, (0.d0,0.d0),            &
                   Aux1, dimTot)

    foo = 0.d0
    do i = 1,dimTot
       foo = foo + Aux1(i,i)
    enddo

    if (ephtype == ephIdx(ntypeunits+1)) then
       do j = idxF(ephtype),idxL(ephtype)
          do i = idxF(ephtype),idxL(ephtype)
             write (2221,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Aux1(i,j)),                                     &
                  DREAL(Symm(i-idxF(ephtype)+1,j-idxF(ephtype)+1)),     &
                  DIMAG(Aux1(i,j)),                                     &
                  DIMAG(Symm(i-idxF(ephtype)+1,j-idxF(ephtype)+1))
             write (2222,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Aux1(i,j)) -                                    &
                  DREAL(Symm(i-idxF(ephtype)+1,j-idxF(ephtype)+1)),     &
                  DIMAG(Aux1(i,j)) -                                    &
                  DIMAG(Symm(i-idxF(ephtype)+1,j-idxF(ephtype)+1))
          enddo
       enddo
    else if (ephtype == ephIdx(ntypeunits+2)) then
! OBS.: Only work if 'norbDyn' is equal to 'NR'!!
       dim = unitdimensions(ntypeunits+2)
       do j = idxF(ephtype),idxL(ephtype)
          do i = idxF(ephtype),idxL(ephtype)
             write (2221,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Aux1(dimTot-dim+i,dimTot-dim+j)),               &
                  DREAL(Symm(i-idxF(ephtype)+1,j-idxF(ephtype)+1)) +    &
                  DREAL(Hc(i-idxF(ephtype)+1,j-idxF(ephtype)+1)),       &
                  DIMAG(Aux1(dimTot-dim+i,dimTot-dim+j)),               &
                  DIMAG(Symm(i-idxF(ephtype)+1,j-idxF(ephtype)+1)) +    &
                  DIMAG(Hc(i-idxF(ephtype)+1,j-idxF(ephtype)+1))
             write (2222,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Aux1(dimTot-dim+i,dimTot-dim+j)) -              &
                  DREAL(Symm(i-idxF(ephtype)+1,j-idxF(ephtype)+1)) -    &
                  DREAL(Hc(i-idxF(ephtype)+1,j-idxF(ephtype)+1)),       &
                  DIMAG(Aux1(dimTot-dim+i,dimTot-dim+j)) -              &
                  DIMAG(Symm(i-idxF(ephtype)+1,j-idxF(ephtype)+1)) -    &
                  DIMAG(Hc(i-idxF(ephtype)+1,j-idxF(ephtype)+1))
          enddo
       enddo
    else
       dim = unitdimensions(ntypeunits+1)
       if (ephIdx(ntypeunits+1) /= 0) then ! first unit has eph?
          ueph = 2 ! eph units indexing
       else
          ueph = 1 ! eph units indexing
       endif
       do k = 2,nunits+1
          utype = unit_type(k) ! current unit type
          if (ephIdx(utype) /= 0) then ! unit with eph?
             if (ephunit == ueph) then ! the unit we are looking?
                do j = idxF(ephtype),idxL(ephtype)
                   do i = idxF(ephtype),idxL(ephtype)
                      write (2221,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')  &
                           DREAL(Aux1(dim+i,dim+j)),                    &
                           DREAL(Symm(i-idxF(ephtype)+1,                &
                                      j-idxF(ephtype)+1)),              &
                           DIMAG(Aux1(dim+i,dim+j)),                    &
                           DIMAG(Symm(i-idxF(ephtype)+1,                &
                                      j-idxF(ephtype)+1))
                      write (2222,'(e17.8e3,e17.8e3)')                  &
                           DREAL(Aux1(dim+i,dim+j)) -                   &
                           DREAL(Symm(i-idxF(ephtype)+1,                &
                                      j-idxF(ephtype)+1)),              &
                           DIMAG(Aux1(dim+i,dim+j)) -                   &
                           DIMAG(Symm(i-idxF(ephtype)+1,                &
                                      j-idxF(ephtype)+1))
                   enddo
                enddo
                EXIT
             else
                ueph = ueph + 1
             endif
          endif
          dim = dim + unitdimensions(utype)
       enddo
    endif
          
!   Free memory.
    deallocate (GrTot)
    deallocate (Gamma_LTot)
    deallocate (Gamma_RTot)
    deallocate (MephTot)
    deallocate (Aux1)
    deallocate (Aux2)
    deallocate (Aux3)


  end subroutine testInelSymm


!  *******************************************************************  !
!                             testInelAsymm                             !
!  *******************************************************************  !
!  Description: build the full hamiltonian of the scattering region     !
!  and compute the asymmetric part of inelastic transmission in order   !
!  compare with the transmission obtained with 'inelSymm' subroutine.   !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer ntypeunits                   : Number of unit types          !
!  integer nunits                       : Total number of units         !
!  integer unit_type(nunits+2)          : Units types                   !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  real*8 unitshift(ntypeunits+2)       : Units shift                   !
!  real*8 S1unit(unitdimensions(ntypeunits),                            !
!                unitdimensions(ntypeunits)) : Unit coupling overlap    !
!  real*8 H1unit(unitdimensions(ntypeunits),                            !
!                unitdimensions(ntypeunits),nspin) : Unit coupling      !
!                                                    Hamiltonian        !
!  TYPE(unitS) Sunits(ntypeunits+2)%S(unitdimensions,unitdimensions) :  !
!                                         [real*8] Units overlap        !
!  TYPE(unitH)                                                          !
!        Hunits(ntypeunits+2)%H(unitdimensions,unitdimensions,nspin) :  !
!                                         [real*8] Units hamiltonian    !
!  complex(8) Sigma_L(NL,NL)           : Left-lead self-energy          !
!  complex(8) Sigma_R(NR,NR)           : Right-lead self-energy         !
!  integer ephIdx(ntypeunits+2)        : Unit index (those with e-ph)   !
!  integer idxF(neph)                  : First dynamic atom orbital     !
!  integer idxL(neph)                  : Last dynamic atom orbital      !
!  ****************************** INPUT ******************************  !
!  real*8 Ei                           : Energy grid point              !
!  integer ispin                       : Spin component index           !
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Gamma_L(NL,NL)           : Left-lead coupling matrix      !
!  complex(8) Gamma_R(NR,NR)           : Right-lead coupling matrix     !
!  integer ephunit                     : E-ph unit index                !
!  integer ephtype                     : E-ph unit type index           !
!  integer norbDyn                     : # of dynamic atoms orbitals    !
!  complex(8) Meph(norbDyn,norbDyn)    : E-ph coupling matrix           !
!  complex(8) Asymm(norbDyn,norbDyn)   : Calculated transmission        !
!  complex(8) Hc(NR,NR)                : Hermitian conjugated part      !
!  *******************************************************************  !
  subroutine testInelAsymm (Ei, ispin, NL, Gamma_L, NR, Gamma_R,        &
                            ephunit, ephtype, norbDyn, Meph, Asymm, Hc)

!
!   Modules
!
    use idsrdr_options,  only: ntypeunits, nunits
    use idsrdr_units,    only: unit_type, unitdimensions, unitshift,    &
                               S1unit, H1unit, Sunits, Hunits
    use idsrdr_leads,    only: Sigma_L, Sigma_R
    use idsrdr_ephcoupl, only: ephIdx, idxF, idxL

!   Input variables.
    integer, intent(in) :: ispin, NL, NR, ephunit, ephtype, norbDyn
    real(8), intent(in) :: Ei
    complex(8), dimension (NL,NL), intent(in) :: Gamma_L
    complex(8), dimension (NR,NR), intent(in) :: Gamma_R
    complex(8), dimension (norbDyn,norbDyn), intent(in) :: Meph
    complex(8), dimension (norbDyn,norbDyn), intent(in) :: Asymm
    complex(8), dimension (NR,NR), intent(in) :: Hc

!   Local variables.
    integer :: i, j, k, utype, dim, dimTot, dimCpl, idxAnt, ueph
    real(8), allocatable, dimension (:,:) :: STot
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8) :: foo
    complex(8), allocatable, dimension (:,:) :: HTot, GrTot, MephTot,   &
                                                Gamma_LTot, Gamma_RTot, &
                                                Aux1, Aux2, Aux3
    external :: HI_zsymm, HI_zhemm, HI_zgemm, HI_zgeInvert

!   Get the total dimension.
    dimTot = unitdimensions(ntypeunits+1) + unitdimensions(ntypeunits+2)
    do I = 2,nunits+1
       utype = unit_type(I) ! current unit type
       dimTot = dimTot + unitdimensions(utype)
    enddo

!   Allocate and initialize matrices.
    allocate (STot(dimTot,dimTot))
    allocate (HTot(dimTot,dimTot))
    allocate (GrTot(dimTot,dimTot))
    STot = 0.d0
    HTot = 0.d0

!   Build total hamiltonian and overlap matrices.
    dimCpl = unitdimensions(ntypeunits) ! coupling unit dimensions
    utype = ntypeunits+1 ! current unit type (first unit)
    dim = unitdimensions(utype) ! current unit dimensions

    STot(1:dim,1:dim) = (Ei-unitshift(utype))*Sunits(utype)%S(:,:)
    HTot(1:dim,1:dim) = Hunits(utype)%H(:,:,ispin)
    STot(dim-dimCpl+1:dim,dim+1:dim+dimCpl) =                           &
         (Ei-unitshift(ntypeunits))*S1unit
    HTot(dim-dimCpl+1:dim,dim+1:dim+dimCpl) = H1unit(:,:,ispin)
    STot(dim+1:dim+dimCpl,dim-dimCpl+1:dim) =                           &
         (Ei-unitshift(ntypeunits))*TRANSPOSE(S1unit)
    HTot(dim+1:dim+dimCpl,dim-dimCpl+1:dim) =                           &
         TRANSPOSE(H1unit(:,:,ispin))

    idxAnt = dim + 1
    do k = 2,nunits+1
       utype = unit_type(k) ! current unit type
       dim = unitdimensions(utype) ! current unit dimensions

       STot(idxAnt:idxAnt+dim-1,idxAnt:idxAnt+dim-1) =                  &
            (Ei-unitshift(utype))*Sunits(utype)%S(:,:)
       HTot(idxAnt:idxAnt+dim-1,idxAnt:idxAnt+dim-1) =                  &
            Hunits(utype)%H(:,:,ispin)
       STot(idxAnt+dim-dimCpl:idxAnt+dim-1,                             &
            idxAnt+dim:idxAnt+dim+dimCpl-1) =                           &
            (Ei-unitshift(ntypeunits))*S1unit
       HTot(idxAnt+dim-dimCpl:idxAnt+dim-1,                             &
            idxAnt+dim:idxAnt+dim+dimCpl-1) = H1unit(:,:,ispin)
       STot(idxAnt+dim:idxAnt+dim+dimCpl-1,                             &
            idxAnt+dim-dimCpl:idxAnt+dim-1) =                           &
            (Ei-unitshift(ntypeunits))*TRANSPOSE(S1unit)
       HTot(idxAnt+dim:idxAnt+dim+dimCpl-1,                             &
            idxAnt+dim-dimCpl:idxAnt+dim-1) =                           &
            TRANSPOSE(H1unit(:,:,ispin))

       idxAnt = idxAnt + dim

    enddo

    utype = ntypeunits + 2 ! current unit type (last unit)
    dim = unitdimensions(utype) ! current unit dimensions

    STot(idxAnt:idxAnt+dim-1,idxAnt:idxAnt+dim-1) =                     &
         (Ei-unitshift(utype))*Sunits(utype)%S(:,:)
    HTot(idxAnt:idxAnt+dim-1,idxAnt:idxAnt+dim-1) =                     &
         Hunits(utype)%H(:,:,ispin)

!   Green's function.
    GrTot = STot - HTot
    GrTot(1:NL,1:NL) = GrTot(1:NL,1:NL) - Sigma_L
    GrTot(dimTot-NR+1:dimTot,dimTot-NR+1:dimTot) =                      &
         GrTot(dimTot-NR+1:dimTot,dimTot-NR+1:dimTot) - Sigma_R
    call HI_zgeInvert (Grtot, dimTot)

!   Free memory.
    deallocate (HTot)
    deallocate (STot)

!   Allocate full matrices.
    allocate (Gamma_LTot(dimTot,dimTot))
    allocate (Gamma_RTot(dimTot,dimTot))
    allocate (MephTot(dimTot,dimTot))

!   Assign the full lead's coupling matrices.
    Gamma_LTot = (0.d0,0.d0)
    Gamma_RTot = (0.d0,0.d0)
    do j = 1,NL
       do i = 1,NL
          Gamma_LTot(i,j) = Gamma_L(i,j)
       enddo
    enddo
    do j = 1,NR
       do i = 1,NR
          Gamma_RTot(dimTot-NR+i,dimTot-NR+j) = Gamma_R(i,j)
       enddo
    enddo

!   Assign the full e-ph coupling matrix.
    MephTot = (0.d0,0.d0)
    if (ephtype == ephIdx(ntypeunits+1)) then
       do j = idxF(ephtype),idxL(ephtype)
          do i = idxF(ephtype),idxL(ephtype)
             MephTot(i,j) = Meph(i-idxF(ephtype)+1,j-idxF(ephtype)+1)
          enddo
       enddo
    else if (ephtype == ephIdx(ntypeunits+2)) then
       dim = unitdimensions(ntypeunits+2)
       do j = idxF(ephtype),idxL(ephtype)
          do i = idxF(ephtype),idxL(ephtype)
             MephTot(dimTot-dim+i,dimTot-dim+j) =                       &
                  Meph(i-idxF(ephtype)+1,j-idxF(ephtype)+1)
          enddo
       enddo
    else
       dim = unitdimensions(ntypeunits+1)
       if (ephIdx(ntypeunits+1) /= 0) then ! first unit has eph?
          ueph = 2 ! eph units indexing
       else
          ueph = 1 ! eph units indexing
       endif
       do k = 2,nunits+1
          utype = unit_type(k) ! current unit type
          if (ephIdx(utype) /= 0) then ! unit with eph?
             if (ephunit == ueph) then ! the unit we are looking?
                do j = idxF(ephtype),idxL(ephtype)
                   do i = idxF(ephtype),idxL(ephtype)
                      MephTot(dim+i,dim+j) =                            &
                           Meph(i-idxF(ephtype)+1,j-idxF(ephtype)+1)
                   enddo
                enddo
                EXIT
             else
                ueph = ueph + 1
             endif
          endif
          dim = dim + unitdimensions(utype)
       enddo
    endif

!   Allocate auxiliaries matrices.
    allocate (Aux1(dimTot,dimTot))
    allocate (Aux2(dimTot,dimTot))
    allocate (Aux3(dimTot,dimTot))

!   -- 1st PART:
!            'Gamma_R*G^dagger*Meph*G*(Gamma_R-Gamma_L)*G^dagger*Meph' --

!   ('Aux1 = Gamma_RTot * GrTot^dagger')
    call HI_zgemm ('N', 'C', dimTot, dimTot, dimTot, (1.d0,0.d0),       &
                   Gamma_RTot, dimTot, GrTot, dimTot,                   &
                   (0.d0,0.d0), Aux1, dimTot)

!   ('Aux2 = Aux1 * MephTot')
    call HI_zsymm ('R', 'L', dimTot, dimTot, (1.d0,0.d0), MephTot,      &
                   dimTot, Aux1, dimTot, (0.d0,0.d0), Aux2, dimTot)

!   ('Aux1 = Aux2 * GrTot')
    call HI_zgemm ('N', 'N', dimTot, dimTot, dimTot, (1.d0,0.d0), Aux2, &
                   dimTot, GrTot, dimTot, (0.d0,0.d0), Aux1, dimTot)

!   ('Aux3 = Gamma_RTot - Gamma_LTot')
    Aux3 = Gamma_RTot - Gamma_LTot

!   ('Aux2 = Aux1 * Aux3')
    call HI_zhemm ('R', 'L', dimTot, dimTot, (1.d0,0.d0), Aux3,         &
                   dimTot, Aux1, dimTot, (0.d0,0.d0), Aux2, dimTot)

!   ('Aux3 = Aux2 * GrTot^dagger')
    call HI_zgemm ('N', 'C', dimTot, dimTot, dimTot, (1.d0,0.d0), Aux2, &
                   dimTot, GrTot, dimTot, (0.d0,0.d0), Aux3, dimTot)

!   ('Aux1 = Aux3 * MephTot')
    call HI_zsymm ('R', 'L', dimTot, dimTot, (1.d0,0.d0), MephTot,      &
                   dimTot, Aux3, dimTot, (0.d0,0.d0), Aux1, dimTot)

!   -- 2nd PART: 'G*(1st PART + H.c.)' --

!   ('Aux3 = GrTot * Aux1')
    call HI_zgemm ('N', 'N', dimTot, dimTot, dimTot, (1.d0,0.d0),       &
                   GrTot, dimTot, Aux1, dimTot, (0.d0,0.d0),            &
                   Aux3, dimTot)

!   ('Aux3 = GrTot * Aux1^dagger + Aux3')
    call HI_zgemm ('N', 'C', dimTot, dimTot, dimTot, (1.d0,0.d0),       &
                   GrTot, dimTot, Aux1, dimTot, (1.d0,0.d0),            &
                   Aux3, dimTot)

!   -- 3rd PART: 'G^dagger*Gamma_L*(2nd PART)' --

!   ('Aux2 = Gamma_LTot * Aux3')
    call HI_zhemm ('L', 'L', dimTot, dimTot, (1.d0,0.d0), Gamma_LTot,   &
                   dimTot, Aux3, dimTot, (0.d0,0.d0), Aux2, dimTot)

!   ('Aux1 = GrTot^dagger * Aux2')
    call HI_zgemm ('C', 'N', dimTot, dimTot, dimTot, (1.d0,0.d0),       &
                   GrTot, dimTot, Aux2, dimTot, (0.d0,0.d0),            &
                   Aux1, dimTot)

    foo = 0.d0
    do i = 1,dimTot
       foo = foo + Aux1(i,i)
    enddo

    if (ephtype == ephIdx(ntypeunits+1)) then
       do j = idxF(ephtype),idxL(ephtype)
          do i = idxF(ephtype),idxL(ephtype)
             write (3331,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Aux1(i,j)),                                     &
                  DREAL(Asymm(i-idxF(ephtype)+1,j-idxF(ephtype)+1)),    &
                  DIMAG(Aux1(i,j)),                                     &
                  DIMAG(Asymm(i-idxF(ephtype)+1,j-idxF(ephtype)+1))
             write (3332,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Aux1(i,j)) -                                    &
                  DREAL(Asymm(i-idxF(ephtype)+1,j-idxF(ephtype)+1)),    &
                  DIMAG(Aux1(i,j)) -                                    &
                  DIMAG(Asymm(i-idxF(ephtype)+1,j-idxF(ephtype)+1))
          enddo
       enddo
    else if (ephtype == ephIdx(ntypeunits+2)) then
! OBS.: Only work if 'norbDyn' is equal to 'NR'!!
       dim = unitdimensions(ntypeunits+2)
       do j = idxF(ephtype),idxL(ephtype)
          do i = idxF(ephtype),idxL(ephtype)
             write (3331,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Aux1(dimTot-dim+i,dimTot-dim+j)),               &
                  DREAL(Asymm(i-idxF(ephtype)+1,j-idxF(ephtype)+1)) +   &
                  DREAL(Hc(i-idxF(ephtype)+1,j-idxF(ephtype)+1)),       &
                  DIMAG(Aux1(dimTot-dim+i,dimTot-dim+j)),               &
                  DIMAG(Asymm(i-idxF(ephtype)+1,j-idxF(ephtype)+1)) +   &
                  DIMAG(Hc(i-idxF(ephtype)+1,j-idxF(ephtype)+1))
             write (3332,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Aux1(dimTot-dim+i,dimTot-dim+j)) -              &
                  DREAL(Asymm(i-idxF(ephtype)+1,j-idxF(ephtype)+1)) -   &
                  DREAL(Hc(i-idxF(ephtype)+1,j-idxF(ephtype)+1)),       &
                  DIMAG(Aux1(dimTot-dim+i,dimTot-dim+j)) -              &
                  DIMAG(Asymm(i-idxF(ephtype)+1,j-idxF(ephtype)+1)) -   &
                  DIMAG(Hc(i-idxF(ephtype)+1,j-idxF(ephtype)+1))
          enddo
       enddo
    else
       dim = unitdimensions(ntypeunits+1)
       if (ephIdx(ntypeunits+1) /= 0) then ! first unit has eph?
          ueph = 2 ! eph units indexing
       else
          ueph = 1 ! eph units indexing
       endif
       do k = 2,nunits+1
          utype = unit_type(k) ! current unit type
          if (ephIdx(utype) /= 0) then ! unit with eph?
             if (ephunit == ueph) then ! the unit we are looking?
                do j = idxF(ephtype),idxL(ephtype)
                   do i = idxF(ephtype),idxL(ephtype)
                      write (3331,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')  &
                           DREAL(Aux1(dim+i,dim+j)),                    &
                           DREAL(Asymm(i-idxF(ephtype)+1,               &
                                       j-idxF(ephtype)+1)),             &
                           DIMAG(Aux1(dim+i,dim+j)),                    &
                           DIMAG(Asymm(i-idxF(ephtype)+1,               &
                                       j-idxF(ephtype)+1))
                      write (3332,'(e17.8e3,e17.8e3)')                  &
                           DREAL(Aux1(dim+i,dim+j)) -                   &
                           DREAL(Asymm(i-idxF(ephtype)+1,               &
                                       j-idxF(ephtype)+1)),             &
                           DIMAG(Aux1(dim+i,dim+j)) -                   &
                           DIMAG(Asymm(i-idxF(ephtype)+1,               &
                                       j-idxF(ephtype)+1))
                   enddo
                enddo
                EXIT
             else
                ueph = ueph + 1
             endif
          endif
          dim = dim + unitdimensions(utype)
       enddo
    endif
          
!   Free memory.
    deallocate (GrTot)
    deallocate (Gamma_LTot)
    deallocate (Gamma_RTot)
    deallocate (MephTot)
    deallocate (Aux1)
    deallocate (Aux2)
    deallocate (Aux3)


  end subroutine testInelAsymm


!  *******************************************************************  !
!                               freecurr                                !
!  *******************************************************************  !
!  Description: free allocated memory.                                  !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *******************************************************************  !
  subroutine freecurr


!   Free memory.
    deallocate (allcurr)


  end subroutine freecurr


!  *******************************************************************  !
!                              sumCalcCurr                              !
!  *******************************************************************  !
!  Description: function for suming two items of type 'calcCurr' (to    !
!  be used at MPI collective operations like 'MPI_Reduce').             !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  ****************************** INPUT ******************************  !
!  TYPE(calcCurr) inCurr(:,:,:)  : 1st item to be added                 !
!  integer len                   : Total length of the items            !
!                                  (1st dim x 2nd dim x 3rd dim)        !
!  integer type                  : Created MPI data type                !
!  ************************** INPUT/OUTPUT ***************************  !
!  TYPE(calcCurr) outCurr(:,:,:) : 2nd item to be added                 !
!  *******************************************************************  !
  subroutine sumCalcCurr (inCurr, outCurr, len, type)

!   Input variables.
    integer, intent(in) :: len, type
    TYPE(calcCurr), intent(in) :: inCurr(len)
    TYPE(calcCurr), intent(inout) :: outCurr(len)

!   Local variables.
    integer :: i

    i = type ! I know that it seems useless...
    do i = 1,len
       outCurr(i)%el = inCurr(i)%el + outCurr(i)%el
       outCurr(i)%isymm = inCurr(i)%isymm + outCurr(i)%isymm
       outCurr(i)%iasymm = inCurr(i)%iasymm + outCurr(i)%iasymm
    enddo


  end subroutine sumCalcCurr


!  *******************************************************************  !


END MODULE idsrdr_current
