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
!                         MODULE idsrdr_current                         !
!  *******************************************************************  !
!  Description: compute the eletronic current (actually the             !
!  transmission coeficient at zero bias).                               !
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
  use parallel,        only: 
  use idsrdr_leads,    only: 
  use idsrdr_ephcoupl, only: 
  use idsrdr_green,    only: 
  use idsrdr_options,  only: 
  use idsrdr_units,    only: 
  use idsrdr_check,    only: 

  implicit none
  
  PUBLIC  :: current
  PRIVATE :: elastic, transmission, inelSymm, inelAsymm, writeTransm,   &
             testInelSymm, testInelAsymm


CONTAINS


!  *******************************************************************  !
!                                current                                !
!  *******************************************************************  !
!  Description: main subroutine for computing the current (actually     !
!  the transmission coeficient at zero bias).                           !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Sigma_L(NL,NL)           : Left-lead self-energy          !
!  complex(8) Sigma_R(NR,NR)           : Right-lead self-energy         !
!  ****************************** INPUT ******************************  !
!  real*8 Ei                           : Energy grid point              !
!  integer ispin                       : Spin component index           !
!  *******************************************************************  !
  subroutine current (Ei, ispin)

!
!   Modules
!
    use idsrdr_leads,    only: NL, NR, Sigma_L, Sigma_R

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    real(8) :: Tel, Tsymm, Tasymm
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8), allocatable, dimension (:,:) :: Gamma_L
    complex(8), allocatable, dimension (:,:) :: Gamma_R

!   Allocate coupling and auxiliary matrices.
    allocate (Gamma_L(NL,NL))
    allocate (Gamma_R(NR,NR))

!   Sets the lead's coupling matrices.
    Gamma_L = zi * (Sigma_L - DCONJG(Sigma_L))
    Gamma_R = zi * (Sigma_R - DCONJG(Sigma_R))

!   Compute elastic contribution.
    call elastic (Tel, NL, Gamma_L, NR, Gamma_R)

!   Compute symmetric part of inelastic contribution.
!   OBS.: change the commmented lines for testing.
    call inelSymm (Tsymm, ispin, NL, Gamma_L, NR, Gamma_R)
!!$    call inelSymm (Tsymm, ispin, NL, Gamma_L, NR, Gamma_R, Ei)

!   Compute asymmetric part of inelastic contribution.
!!$    call inelAsymm (Tasymm, ispin, NL, Gamma_L, NR, Gamma_R)
    call inelAsymm (Tasymm, ispin, NL, Gamma_L, NR, Gamma_R, Ei)

!   Write transmissions to outputfiles.
    call writeTransm (Ei, Tel, Tsymm, Tasymm)

!   Free memory.
    deallocate (Gamma_L)
    deallocate (Gamma_R)


  end subroutine current


!  *******************************************************************  !
!                                elastic                                !
!  *******************************************************************  !
!  Description: compute the elastic transmission.                       !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode                      : True if it is the I/O node     !
!  complex(8) Gr_1M(NL,NR)             : G^r_{1,M}                      !
!  ****************************** INPUT ******************************  !
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Gamma_L(NL,NL)           : Left-lead coupling matrix      !
!  complex(8) Gamma_R(NR,NR)           : Right-lead coupling matrix     !
!  ***************************** OUTPUT ******************************  !
!  real*8 Tel                          : Elastic conductance            !
!  *******************************************************************  !
  subroutine elastic (Tel, NL, Gamma_L, NR, Gamma_R)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_green,    only: Gr_1M

!   Input variables.
    integer, intent(in) :: NL, NR
    real(8), intent(out) :: Tel
    complex(8), dimension (NL,NL), intent(in) :: Gamma_L
    complex(8), dimension (NR,NR), intent(in) :: Gamma_R

!   Local variables.

    if (IOnode) write (6,'(a)', advance='no')                           &
         '      computing elastic current... '

!   Calculates the transmission coefficient.
    call transmission (NL, Gamma_L, NR, Gamma_R, Gr_1M, Tel)

    if (IOnode) write(6,'(a)') " ok!"


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
    external :: zhemm

!   Allocate auxiliary matrices.
    allocate (IJ(NI,NJ))
    allocate (JI(NJ,NI))

!   ('IJ = Gamma_I * Gr_ij')
    call zhemm ('L', 'L', NI, NJ, (1.d0,0.d0), Gamma_I, NI,             &
                Gr_ij, NI, (0.d0,0.d0), IJ, NI)

!   ('JI^dagger = Gr_ij * Gamma_J')
    call zhemm ('R', 'L', NI, NJ, (1.d0,0.d0), Gamma_J, NJ,             &
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
!  logical IOnode              : True if it is the I/O node             !
!  integer neph                : Number of units with e-ph interaction  !
!  integer nModes(neph)        : Number of vibrational modes            !
!  integer norbDyn(neph)       : Number of orbitals from dynamic atoms  !
!  TYPE(ephCplng) Meph(neph)%M(norbDyn,norbDyn,nspin,nModes) :          !
!                                   [complex*8] E-ph coupling matrices  !
!  TYPE(green) Gr_Mn(neph)%G(NR,unitdimensions) :  [complex] G^r_{M,n}  !
!  TYPE(green) Gr_1n(neph)%G(NL,unitdimensions) :  [complex] G^r_{1,n}  !
!  TYPE(green) Gr_nn(neph)%G(unitdimensions,unitdimensions) :           !
!                                                  [complex] G^r_{n,n}  !
!  complex(8) Gr_1M(NL,NR)                      : G^r_{1,M}             !
!  ****************************** INPUT ******************************  !
!  integer ispin                       : Spin component index           !
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Gamma_L(NL,NL)           : Left-lead coupling matrix      !
!  complex(8) Gamma_R(NR,NR)           : Right-lead coupling matrix     !
!  real*8 Ei                           : [optional] Energy grid point   !
!  ***************************** OUTPUT ******************************  !
!  real*8 Tsymm            : Symmetric part of inelastic transmission   !
!  *******************************************************************  !
  subroutine inelSymm (Tsymm, ispin, NL, Gamma_L, NR, Gamma_R, Ei)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_ephcoupl, only: neph, nModes, norbDyn, Meph
    use idsrdr_green,    only: Gr_Mn, Gr_1n, Gr_nn, Gr_1M

!   Input variables.
    integer, intent(in) :: ispin, NL, NR
    real(8), intent(out) :: Tsymm
    real(8), optional, intent(in) :: Ei
    complex(8), dimension (NL,NL), intent(in) :: Gamma_L
    complex(8), dimension (NR,NR), intent(in) :: Gamma_R

!   Local variables.
    integer :: j, w, i
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8), allocatable, dimension(:,:) :: Aux1, Aux2, Aux3, Aux4,  &
                                               Aux5, Aux6, Aux7, GrCJG, A
    external :: zsymm, zhemm, zgemm

    if (IOnode) write (6,'(a)', advance='no')                           &
         '      computing inelastic transmission (symmetric part)... '

!   Initialize variable.
    Tsymm = 0.d0

    do j = 1,neph ! over unit with e-ph

!      Allocate auxiliary matrices.
       allocate (GrCJG(NR,norbDyn(j)))
       allocate (A(norbDyn(j),norbDyn(j)))
       allocate (Aux1(NR,norbDyn(j)))
       allocate (Aux2(NR,norbDyn(j)))
       allocate (Aux3(norbDyn(j),norbDyn(j)))
       allocate (Aux4(norbDyn(j),norbDyn(j)))
       allocate (Aux5(NL,norbDyn(j)))
       allocate (Aux6(NL,norbDyn(j)))
       allocate (Aux7(NR,NR))

!      Copy the complex conjugate of 'Gr_Mn'.
       GrCJG = DCONJG(Gr_Mn(j)%G)

!      Spectral matrix (obs.: 'Gr_nn' is symmetric).
       A = zi * (Gr_nn(j)%G - DCONJG(Gr_nn(j)%G))

       do w = 1,nModes(j) ! over phonon modes

!         -- 1st PART: 'G*Meph*G*Gamma_R*G^dagger*Meph' --

!         ('Aux1 = Gamma_R * Gr_Mn^*')
          call zhemm ('L', 'L', NR, norbDyn(j), (1.d0,0.d0),            &
                      Gamma_R, NR, GrCJG, NR, (0.d0,0.d0), Aux1, NR)

!         ('Aux2 = Aux1 * Meph')
          call zsymm ('R', 'L', NR, norbDyn(j), (1.d0,0.d0),            &
                      Meph(j)%M(:,:,ispin,w), norbDyn(j),               &
                      Aux1, NR, (0.d0,0.d0), Aux2, NR)

!         ('Aux3 = Gr_Mn^T * Aux2')
          call zgemm ('T', 'N', norbDyn(j), norbDyn(j), NR,             &
                      (1.d0,0.d0), Gr_Mn(j)%G, NR, Aux2, NR,            &
                      (0.d0,0.d0), Aux3, norbDyn(j))

!         ('Aux4 = Meph * Aux3')
          call zsymm ('L', 'L', norbDyn(j), norbDyn(j), (1.d0,0.d0),    &
                      Meph(j)%M(:,:,ispin,w), norbDyn(j),               &
                      Aux3, norbDyn(j), (0.d0,0.d0), Aux4, norbDyn(j))

!         ('Aux5 = Gr_1n * Aux4')
          call zgemm ('N', 'N', NL, norbDyn(j), norbDyn(j),             &
                      (1.d0,0.d0), Gr_1n(j)%G, NL, Aux4,                &
                      norbDyn(j), (0.d0,0.d0), Aux5, NL)

!         -- 2nd PART: 'Gamma_R*G^dagger*Meph*A*Meph' --

!         ('Aux1 = Aux2 * A') (where 'Aux2 = Gamma_R * Gr_Mn^* * Meph')
          call zhemm ('R', 'L', NR, norbDyn(j), (1.d0,0.d0),            &
                      A, norbDyn(j), Aux2, NR, (0.d0,0.d0), Aux1, NR)

!         ('Aux2 = Aux1 * Meph')
          call zsymm ('R', 'L', NR, norbDyn(j), (1.d0,0.d0),            &
                      Meph(j)%M(:,:,ispin,w), norbDyn(j),               &
                      Aux1, NR, (0.d0,0.d0), Aux2, NR)

!         -- 3rd PART: 'G^dagger*Gamma_L*(1st PART + i/2*G*2nd PART)' --

!         ('Aux5 = i/2 * Gr_1M * Aux2 + Aux5')
          call zgemm ('N', 'N', NL, norbDyn(j), NR, (0.d0,0.5d0),       &
                      Gr_1M, NL, Aux2, NR, (1.d0,0.d0), Aux5, NL)

!         ('Aux6 = Gamma_L * Aux5')
          call zhemm ('L', 'L', NL, norbDyn(j), (1.d0,0.d0),            &
                      Gamma_L, NL, Aux5, NL, (0.d0,0.d0), Aux6, NL)

!         ('Aux4 = Gr_1n^dagger * Aux6')
          call zgemm ('C', 'N', norbDyn(j), norbDyn(j), NL,             &
                      (1.d0,0.d0), Gr_1n(j)%G, NL, Aux6, NL,            &
                      (0.d0,0.d0), Aux4, norbDyn(j))

!         -- 4rd PART: 'G^dagger*Gamma_L*G*(-i/2 * H.c.)' --

!         ('Aux6 = Gamma_L * Gr_1n')
          call zhemm ('L', 'L', NL, norbDyn(j), (1.d0,0.d0), Gamma_L,   &
                      NL, Gr_1n(j)%G, NL, (0.d0,0.d0), Aux6, NL)

!         ('Aux1 = Gr_1M^dagger * Aux6')
          call zgemm ('C', 'N', NR, norbDyn(j), NL, (1.d0,0.d0),        &
                      Gr_1M, NL, Aux6, NL, (0.d0,0.d0), Aux1, NR)

!         ('Aux7 = Aux1 * (-i/2 * H.c.)') (where 'H.c. = Aux2^dagger')
          call zgemm ('N', 'C', NR, NR, norbDyn(j), (0.d0,-0.5d0),      &
                      Aux1, NR, Aux2, NR, (0.d0,0.d0), Aux7, NR)

!         Compute the trace.
          do i = 1,norbDyn(j)
             Tsymm = Tsymm + DREAL(Aux4(i,i))
          enddo
          do i = 1,NR
             Tsymm = Tsymm + DREAL(Aux7(i,i))
          enddo
          
!         [test] Compute the matrices multiplication with full matrices.
          if (present(Ei)) then
             call testInelSymm (Ei, ispin, NL, Gamma_L, NR, Gamma_R,    &
                                j, norbDyn(j), Meph(j)%M(:,:,ispin,w),  &
                                Aux4, Aux7)
          endif

       enddo ! do w = 1,nModes(j)

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

    enddo ! do j = 1,neph

    if (IOnode) write(6,'(a)') " ok!"


  end subroutine inelSymm


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
!  logical IOnode              : True if it is the I/O node             !
!  integer neph                : Number of units with e-ph interaction  !
!  integer nModes(neph)        : Number of vibrational modes            !
!  integer norbDyn(neph)       : Number of orbitals from dynamic atoms  !
!  TYPE(ephCplng) Meph(neph)%M(norbDyn,norbDyn,nspin,nModes) :          !
!                                   [complex*8] E-ph coupling matrices  !
!  TYPE(green) Gr_Mn(neph)%G(NR,unitdimensions) :  [complex] G^r_{M,n}  !
!  TYPE(green) Gr_1n(neph)%G(NL,unitdimensions) :  [complex] G^r_{1,n}  !
!  complex(8) Gr_1M(NL,NR)                      : G^r_{1,M}             !
!  ****************************** INPUT ******************************  !
!  integer ispin                       : Spin component index           !
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Gamma_L(NL,NL)           : Left-lead coupling matrix      !
!  complex(8) Gamma_R(NR,NR)           : Right-lead coupling matrix     !
!  real*8 Ei                           : [optional] Energy grid point   !
!  ***************************** OUTPUT ******************************  !
!  real*8 Tasymm           : Asymmetric part of inelastic transmission  !
!  *******************************************************************  !
  subroutine inelAsymm (Tasymm, ispin, NL, Gamma_L, NR, Gamma_R, Ei)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_ephcoupl, only: neph, nModes, norbDyn, Meph
    use idsrdr_green,    only: Gr_Mn, Gr_1n, Gr_1M

!   Input variables.
    integer, intent(in) :: ispin, NL, NR
    real(8), intent(out) :: Tasymm
    real(8), optional, intent(in) :: Ei
    complex(8), dimension (NL,NL), intent(in) :: Gamma_L
    complex(8), dimension (NR,NR), intent(in) :: Gamma_R

!   Local variables.
    integer :: j, w, i
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8), allocatable, dimension(:,:) :: Aux1, Aux2, Aux3,        &
                                               Aux4, Aux5, Aux6,        &
                                               Gr_MnCJG, Gr_1nCJG
    external :: zsymm, zhemm, zgemm

    if (IOnode) write (6,'(a)', advance='no')                           &
         '      computing inelastic transmission (asymmetric part)... '

!   Initialize variable.
    Tasymm = 0.d0

    do j = 1,neph ! over unit with e-ph

!      Allocate auxiliary matrices.
       allocate (Gr_MnCJG(NR,norbDyn(j)))
       allocate (Gr_1nCJG(NL,norbDyn(j)))
       allocate (Aux1(NR,norbDyn(j)))
       allocate (Aux2(NR,norbDyn(j)))
       allocate (Aux3(norbDyn(j),norbDyn(j)))
       allocate (Aux4(NL,norbDyn(j)))
       allocate (Aux5(NL,norbDyn(j)))
       allocate (Aux6(NR,NR))

!      Copy the complex conjugate of 'Gr_Mn' and 'Gr_1n'.
       Gr_MnCJG = DCONJG(Gr_Mn(j)%G)
       Gr_1nCJG = DCONJG(Gr_1n(j)%G)

       do w = 1,nModes(j) ! over phonon modes

!         -- 1st PART: 'G*Gamma_R*G^dagger*Meph' --

!         ('Aux1 = Gamma_R * Gr_Mn^*')
          call zhemm ('L', 'L', NR, norbDyn(j), (1.d0,0.d0),            &
                      Gamma_R, NR, Gr_MnCJG, NR, (0.d0,0.d0), Aux1, NR)

!         ('Aux2 = Aux1 * Meph')
          call zsymm ('R', 'L', NR, norbDyn(j), (1.d0,0.d0),            &
                      Meph(j)%M(:,:,ispin,w), norbDyn(j),               &
                      Aux1, NR, (0.d0,0.d0), Aux2, NR)

!         ('Aux3 = Gr_Mn^T * Aux2')
          call zgemm ('T', 'N', norbDyn(j), norbDyn(j), NR,             &
                      (1.d0,0.d0), Gr_Mn(j)%G, NR, Aux2, NR,            &
                      (0.d0,0.d0), Aux3, norbDyn(j))

!         -- 2nd PART: '-G*Gamma_L*G^dagger*Meph + 1st PART' --

!         ('Aux4 = Gr_1n^* * Meph')
          call zsymm ('R', 'L', NL, norbDyn(j), (1.d0,0.d0),            &
                      Meph(j)%M(:,:,ispin,w), norbDyn(j),               &
                      Gr_1nCJG, NL, (0.d0,0.d0), Aux4, NL)

!         ('Aux5 = Gamma_L * Aux4')
          call zhemm ('L', 'L', NL, norbDyn(j), (1.d0,0.d0),            &
                      Gamma_L, NL, Aux4, NL, (0.d0,0.d0), Aux5, NL)

!         ('Aux3 = - Gr_1n^T * Aux5 + Aux3')
          call zgemm ('T', 'N', norbDyn(j), norbDyn(j), NL,             &
                      (-1.d0,0.d0), Gr_1n(j)%G, NL, Aux5, NL,           &
                      (1.d0,0.d0), Aux3, norbDyn(j))

!         -- 3rd PART:
!                'G^dagger*Gamma_L*G*Gamma_R*G^dagger*Meph*(2nd Part)' --

!         ('Aux1 = Aux2 * Aux3') (where 'Aux2= Gamma_R * Gr_Mn^* * Meph')
          call zgemm ('N', 'N', NR, norbDyn(j), norbDyn(j),             &
                      (1.d0,0.d0), Aux2, NR, Aux3,                      &
                      norbDyn(j), (0.d0,0.d0), Aux1, NR)

!         ('Aux4 = Gr_1M * Aux1')
          call zgemm ('N', 'N', NL, norbDyn(j), NR, (1.d0,0.d0),        &
                      Gr_1M, NL, Aux1, NR, (0.d0,0.d0), Aux4, NL)

!         ('Aux5 = Gamma_L*Aux4')
          call zhemm ('L', 'L', NL, norbDyn(j), (1.d0,0.d0), Gamma_L,   &
                      NL, Aux4, NL, (0.d0,0.d0), Aux5, NL)

!         ('Aux3 = Gr_1n^dagger * Aux5')
          call zgemm ('C', 'N', norbDyn(j), norbDyn(j), NL,             &
                      (1.d0,0.d0), Gr_1n(j)%G, NL, Aux5, NL,            &
                      (0.d0,0.d0), Aux3, norbDyn(j))

!         -- 4rd PART: 'G^dagger*Gamma_L*G*(H.c.)' --

!         ('Aux4 = Gamma_L * Gr_1n')
          call zhemm ('L', 'L', NL, norbDyn(j), (1.d0,0.d0), Gamma_L,   &
                      NL, Gr_1n(j)%G, NL, (0.d0,0.d0), Aux4, NL)

!         ('Aux2 = Gr_1M^dagger * Aux4')
          call zgemm ('C', 'N', NR, norbDyn(j), NL, (1.d0,0.d0),        &
                      Gr_1M, NL, Aux4, NL, (0.d0,0.d0), Aux2, NR)

!         ('Aux6 = Aux2 * (H.c.)') ('H.c. = Aux1^dagger')
          call zgemm ('N', 'C', NR, NR, norbDyn(j), (1.d0,0.0d0),      &
                      Aux2, NR, Aux1, NR, (0.d0,0.d0), Aux6, NR)

!         Compute the trace.
          do i = 1,norbDyn(j)
             Tasymm = Tasymm + DREAL(Aux3(i,i))
          enddo
          do i = 1,NR
             Tasymm = Tasymm + DREAL(Aux6(i,i))
          enddo

!         [test] Compute the matrices multiplication with full matrices.
          if (present(Ei)) then
             call testInelAsymm (Ei, ispin, NL, Gamma_L, NR, Gamma_R,   &
                                 j, norbDyn(j), Meph(j)%M(:,:,ispin,w), &
                                 Aux3, Aux6)
          endif
          
       enddo ! do w = 1,nModes(j)

!      Free memory.
       deallocate (Gr_MnCJG)
       deallocate (Gr_1nCJG)
       deallocate (Aux1)
       deallocate (Aux2)
       deallocate (Aux3)
       deallocate (Aux4)
       deallocate (Aux5)
       deallocate (Aux6)

    enddo ! do j = 1,neph

    if (IOnode) write(6,'(a)') " ok!"


  end subroutine inelAsymm


!  *******************************************************************  !
!                              writeTransm                              !
!  *******************************************************************  !
!  Description: write the transmissions (components and total) to       !
!  output file.                                                         !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  ****************************** INPUT ******************************  !
!  real*8 Ei               : Energy grid point                          !
!  real*8 Tel              : Elastic transmission                       !
!  real*8 Tsymm            : Symmetric part of inelastic transmission   !
!  real*8 Tasymm           : Asymmetric part of inelastic transmission  !
!  *******************************************************************  !
  subroutine writeTransm (Ei, Tel, Tsymm, Tasymm)

!
!   Modules
!
    use parallel,        only: IOnode

    include "mpif.h"

!   Input variables.
    real(8), intent(in) :: Ei, Tel, Tsymm, Tasymm

!   Local variables.
    real(8) :: TelTot, TsymmTot, TasymmTot
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

!   Sum the computed transmissions from all nodes.
    call MPI_Reduce (Tel, TelTot, 1, MPI_Double_Precision,              &
                     MPI_Sum, 0, MPI_Comm_World, MPIerror)
    call MPI_Reduce (Tsymm, TsymmTot, 1, MPI_Double_Precision,          &
                     MPI_Sum, 0, MPI_Comm_World, MPIerror)
    call MPI_Reduce (Tasymm, TasymmTot, 1, MPI_Double_Precision,        &
                     MPI_Sum, 0, MPI_Comm_World, MPIerror)

    if (IOnode) then

       write (1102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')                 &
              Ei, Tel, -Tsymm+Tasymm, Tel-Tsymm+Tasymm
       

    endif


  end subroutine writeTransm


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
!  integer norbDyn                     : # of dynamic atoms orbitals    !
!  complex(8) Meph(norbDyn,norbDyn)    : E-ph coupling matrix           !
!  complex(8) Symm(norbDyn,norbDyn)    : Calculated transmission        !
!  complex(8) Hc(NR,NR)                : Hermitian conjugated part      !
!  *******************************************************************  !
  subroutine testInelSymm (Ei, ispin, NL, Gamma_L, NR, Gamma_R,         &
                           ephunit, norbDyn, Meph, Symm, Hc)

!
!   Modules
!
    use idsrdr_options,  only: ntypeunits, nunits
    use idsrdr_units,    only: unit_type, unitdimensions, unitshift,    &
                               S1unit, H1unit, Sunits, Hunits
    use idsrdr_leads,    only: Sigma_L, Sigma_R
    use idsrdr_ephcoupl, only: ephIdx, idxF, idxL
    use idsrdr_check,    only: CHECKzgetrf, CHECKzgetri

!   Input variables.
    integer, intent(in) :: ispin, NL, NR, ephunit, norbDyn
    real(8), intent(in) :: Ei
    complex(8), dimension (NL,NL), intent(in) :: Gamma_L
    complex(8), dimension (NR,NR), intent(in) :: Gamma_R
    complex(8), dimension (norbDyn,norbDyn), intent(in) :: Meph
    complex(8), dimension (norbDyn,norbDyn), intent(in) :: Symm
    complex(8), dimension (NR,NR), intent(in) :: Hc

!   Local variables.
    integer :: i, j, k, utype, dim, dimTot, dimCpl, idxAnt
    integer, allocatable, dimension (:) :: ipiv
    real(8), allocatable, dimension (:,:) :: STot
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8), allocatable, dimension (:,:) :: HTot, GrTot, MephTot,   &
                                                Gamma_LTot, Gamma_RTot, &
                                                Aux1, Aux2, Aux3
    external :: zsymm, zhemm, zgemm

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
    allocate (ipiv(dimTot))
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
    call CHECKzgetrf (dimTot, GrTot, ipiv)
    call CHECKzgetri (dimTot, GrTot, ipiv)

!   Free memory.
    deallocate (HTot)
    deallocate (STot)
    deallocate (ipiv)

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
    if (ephunit == ephIdx(ntypeunits+1)) then
       do j = idxF(ephunit),idxL(ephunit)
          do i = idxF(ephunit),idxL(ephunit)
             MephTot(i,j) = Meph(i-idxF(ephunit)+1,j-idxF(ephunit)+1)
          enddo
       enddo
    else if (ephunit == ephIdx(ntypeunits+2)) then
       dim = unitdimensions(ntypeunits+2)
       do j = idxF(ephunit),idxL(ephunit)
          do i = idxF(ephunit),idxL(ephunit)
             MephTot(dimTot-dim+i,dimTot-dim+j) =                       &
                  Meph(i-idxF(ephunit)+1,j-idxF(ephunit)+1)
          enddo
       enddo
    else
       dim = unitdimensions(ntypeunits+1)
       do k = 2,nunits+1
          utype = unit_type(k) ! current unit type
          if (ephunit == ephIdx(utype)) then
             do j = idxF(ephunit),idxL(ephunit)
                do i = idxF(ephunit),idxL(ephunit)
                   MephTot(dim+i,dim+j) =                               &
                        Meph(i-idxF(ephunit)+1,j-idxF(ephunit)+1)
                enddo
             enddo
             EXIT
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
    call zgemm ('N', 'C', dimTot, dimTot, dimTot, (1.d0,0.d0),          &
                Gamma_RTot, dimTot, GrTot, dimTot,                      &
                (0.d0,0.d0), Aux1, dimTot)

!   ('Aux2 = Aux1 * MephTot')
    call zsymm ('R', 'L', dimTot, dimTot, (1.d0,0.d0), MephTot,         &
                dimTot, Aux1, dimTot, (0.d0,0.d0), Aux2, dimTot)

!   Spectral matrix (obs.: 'GrTot' is symmetric).
    Aux1 = zi * (GrTot - DCONJG(GrTot))

!   ('Aux3 = Aux2 * Aux1')
    call zhemm ('R', 'L', dimTot, dimTot, (1.d0,0.d0), Aux1,            &
                dimTot, Aux2, dimTot, (0.d0,0.d0), Aux3, dimTot)

!   ('Aux1 = Aux3 * MephTot')
    call zsymm ('R', 'L', dimTot, dimTot, (1.d0,0.d0), MephTot,         &
                dimTot, Aux3, dimTot, (0.d0,0.d0), Aux1, dimTot)

!   ('Aux3 = i/2 * GrTot * Aux1')
    call zgemm ('N', 'N', dimTot, dimTot, dimTot, (0.d0,0.5d0), GrTot,  &
                dimTot, Aux1, dimTot, (0.d0,0.d0), Aux3, dimTot)

!   ('Aux3 = -i/2 * GrTot * Aux1^dagger + Aux3')
    call zgemm ('N', 'C', dimTot, dimTot, dimTot, (0.d0,-0.5d0), GrTot, &
                dimTot, Aux1, dimTot, (1.d0,0.d0), Aux3, dimTot)

!   -- 2nd PART: 'G*Meph*G*Gamma_R*G^dagger*Meph + 1st PART' --

!   ('Aux1 = GrTot * Aux2')
!                    (where 'Aux2 = Gamma_RTot * GrTot^dagger * MephTot')
    call zgemm ('N', 'N', dimTot, dimTot, dimTot, (1.d0,0.d0), GrTot,   &
                dimTot, Aux2, dimTot, (0.d0,0.d0), Aux1, dimTot)

!   ('Aux2 = MephTot * Aux1')
    call zsymm ('L', 'L', dimTot, dimTot, (1.d0,0.d0), MephTot,         &
                dimTot, Aux1, dimTot, (0.d0,0.d0), Aux2, dimTot)

!   ('Aux3 = GrTot * Aux2 + Aux3')
    call zgemm ('N', 'N', dimTot, dimTot, dimTot, (1.d0,0.d0), GrTot,   &
                dimTot, Aux2, dimTot, (1.d0,0.d0), Aux3, dimTot)

!   -- 3rd PART: 'G^dagger*Gamma_L*(2nd PART)' --

!   ('Aux2 = Gamma_LTot * Aux3')
    call zhemm ('L', 'L', dimTot, dimTot, (1.d0,0.d0), Gamma_LTot,      &
                dimTot, Aux3, dimTot, (0.d0,0.d0), Aux2, dimTot)

!   ('Aux1 = GrTot^dagger * Aux2')
    call zgemm ('C', 'N', dimTot, dimTot, dimTot, (1.d0,0.d0),          &
                GrTot, dimTot, Aux2, dimTot, (0.d0,0.d0), Aux1, dimTot)

    if (ephunit == ephIdx(ntypeunits+1)) then
       do j = idxF(ephunit),idxL(ephunit)
          do i = idxF(ephunit),idxL(ephunit)
             write (2221,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Aux1(i,j)),                                     &
                  DREAL(Symm(i-idxF(ephunit)+1,j-idxF(ephunit)+1)),     &
                  DIMAG(Aux1(i,j)),                                     &
                  DIMAG(Symm(i-idxF(ephunit)+1,j-idxF(ephunit)+1))
             write (2222,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Aux1(i,j)) -                                    &
                  DREAL(Symm(i-idxF(ephunit)+1,j-idxF(ephunit)+1)),     &
                  DIMAG(Aux1(i,j)) -                                    &
                  DIMAG(Symm(i-idxF(ephunit)+1,j-idxF(ephunit)+1))
          enddo
       enddo
    else if (ephunit == ephIdx(ntypeunits+2)) then
! OBS.: Only work if 'norbDyn' is equal to 'NR'!!
       dim = unitdimensions(ntypeunits+2)
       do j = idxF(ephunit),idxL(ephunit)
          do i = idxF(ephunit),idxL(ephunit)
             write (2221,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Aux1(dimTot-dim+i,dimTot-dim+j)),               &
                  DREAL(Symm(i-idxF(ephunit)+1,j-idxF(ephunit)+1)) +    &
                  DREAL(Hc(i-idxF(ephunit)+1,j-idxF(ephunit)+1)),       &
                  DIMAG(Aux1(dimTot-dim+i,dimTot-dim+j)),               &
                  DIMAG(Symm(i-idxF(ephunit)+1,j-idxF(ephunit)+1)) +    &
                  DIMAG(Hc(i-idxF(ephunit)+1,j-idxF(ephunit)+1))
             write (2222,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Aux1(dimTot-dim+i,dimTot-dim+j)) -              &
                  DREAL(Symm(i-idxF(ephunit)+1,j-idxF(ephunit)+1)) -    &
                  DREAL(Hc(i-idxF(ephunit)+1,j-idxF(ephunit)+1)),       &
                  DIMAG(Aux1(dimTot-dim+i,dimTot-dim+j)) -              &
                  DIMAG(Symm(i-idxF(ephunit)+1,j-idxF(ephunit)+1)) -    &
                  DIMAG(Hc(i-idxF(ephunit)+1,j-idxF(ephunit)+1))
          enddo
       enddo
    else
       dim = unitdimensions(ntypeunits+1)
       do k = 2,nunits+1
          utype = unit_type(k) ! current unit type
          if (ephunit == ephIdx(utype)) then
             do j = idxF(ephunit),idxL(ephunit)
                do i = idxF(ephunit),idxL(ephunit)
                   write (2221,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')     &
                      DREAL(Aux1(dim+i,dim+j)),                         &
                      DREAL(Symm(i-idxF(ephunit)+1,j-idxF(ephunit)+1)), &
                      DIMAG(Aux1(dim+i,dim+j)),                         &
                      DIMAG(Symm(i-idxF(ephunit)+1,j-idxF(ephunit)+1))
                   write (2222,'(e17.8e3,e17.8e3)')                     &
                      DREAL(Aux1(dim+i,dim+j)) -                        &
                      DREAL(Symm(i-idxF(ephunit)+1,j-idxF(ephunit)+1)), &
                      DIMAG(Aux1(dim+i,dim+j)) -                        &
                      DIMAG(Symm(i-idxF(ephunit)+1,j-idxF(ephunit)+1))
                enddo
             enddo
             EXIT
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
!  integer norbDyn                     : # of dynamic atoms orbitals    !
!  complex(8) Meph(norbDyn,norbDyn)    : E-ph coupling matrix           !
!  complex(8) Asymm(norbDyn,norbDyn)   : Calculated transmission        !
!  complex(8) Hc(NR,NR)                : Hermitian conjugated part      !
!  *******************************************************************  !
  subroutine testInelAsymm (Ei, ispin, NL, Gamma_L, NR, Gamma_R,        &
                            ephunit, norbDyn, Meph, Asymm, Hc)

!
!   Modules
!
    use idsrdr_options,  only: ntypeunits, nunits
    use idsrdr_units,    only: unit_type, unitdimensions, unitshift,    &
                               S1unit, H1unit, Sunits, Hunits
    use idsrdr_leads,    only: Sigma_L, Sigma_R
    use idsrdr_ephcoupl, only: ephIdx, idxF, idxL
    use idsrdr_check,    only: CHECKzgetrf, CHECKzgetri

!   Input variables.
    integer, intent(in) :: ispin, NL, NR, ephunit, norbDyn
    real(8), intent(in) :: Ei
    complex(8), dimension (NL,NL), intent(in) :: Gamma_L
    complex(8), dimension (NR,NR), intent(in) :: Gamma_R
    complex(8), dimension (norbDyn,norbDyn), intent(in) :: Meph
    complex(8), dimension (norbDyn,norbDyn), intent(in) :: Asymm
    complex(8), dimension (NR,NR), intent(in) :: Hc

!   Local variables.
    integer :: i, j, k, utype, dim, dimTot, dimCpl, idxAnt
    integer, allocatable, dimension (:) :: ipiv
    real(8), allocatable, dimension (:,:) :: STot
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8), allocatable, dimension (:,:) :: HTot, GrTot, MephTot,   &
                                                Gamma_LTot, Gamma_RTot, &
                                                Aux1, Aux2, Aux3
    external :: zsymm, zhemm, zgemm

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
    allocate (ipiv(dimTot))
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
    call CHECKzgetrf (dimTot, GrTot, ipiv)
    call CHECKzgetri (dimTot, GrTot, ipiv)

!   Free memory.
    deallocate (HTot)
    deallocate (STot)
    deallocate (ipiv)

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
    if (ephunit == ephIdx(ntypeunits+1)) then
       do j = idxF(ephunit),idxL(ephunit)
          do i = idxF(ephunit),idxL(ephunit)
             MephTot(i,j) = Meph(i-idxF(ephunit)+1,j-idxF(ephunit)+1)
          enddo
       enddo
    else if (ephunit == ephIdx(ntypeunits+2)) then
       dim = unitdimensions(ntypeunits+2)
       do j = idxF(ephunit),idxL(ephunit)
          do i = idxF(ephunit),idxL(ephunit)
             MephTot(dimTot-dim+i,dimTot-dim+j) =                       &
                  Meph(i-idxF(ephunit)+1,j-idxF(ephunit)+1)
          enddo
       enddo
    else
       dim = unitdimensions(ntypeunits+1)
       do k = 2,nunits+1
          utype = unit_type(k) ! current unit type
          if (ephunit == ephIdx(utype)) then
             do j = idxF(ephunit),idxL(ephunit)
                do i = idxF(ephunit),idxL(ephunit)
                   MephTot(dim+i,dim+j) =                               &
                        Meph(i-idxF(ephunit)+1,j-idxF(ephunit)+1)
                enddo
             enddo
             EXIT
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
    call zgemm ('N', 'C', dimTot, dimTot, dimTot, (1.d0,0.d0),          &
                Gamma_RTot, dimTot, GrTot, dimTot,                      &
                (0.d0,0.d0), Aux1, dimTot)

!   ('Aux2 = Aux1 * MephTot')
    call zsymm ('R', 'L', dimTot, dimTot, (1.d0,0.d0), MephTot,         &
                dimTot, Aux1, dimTot, (0.d0,0.d0), Aux2, dimTot)

!   ('Aux1 = Aux2 * GrTot')
    call zgemm ('N', 'N', dimTot, dimTot, dimTot, (1.d0,0.d0), Aux2,    &
                dimTot, GrTot, dimTot, (0.d0,0.d0), Aux1, dimTot)

!   ('Aux3 = Gamma_RTot - Gamma_LTot')
    Aux3 = Gamma_RTot - Gamma_LTot

!   ('Aux2 = Aux1 * Aux3')
    call zhemm ('R', 'L', dimTot, dimTot, (1.d0,0.d0), Aux3,            &
                dimTot, Aux1, dimTot, (0.d0,0.d0), Aux2, dimTot)

!   ('Aux3 = Aux2 * GrTot^dagger')
    call zgemm ('N', 'C', dimTot, dimTot, dimTot, (1.d0,0.d0), Aux2,    &
                dimTot, GrTot, dimTot, (0.d0,0.d0), Aux3, dimTot)

!   ('Aux1 = Aux3 * MephTot')
    call zsymm ('R', 'L', dimTot, dimTot, (1.d0,0.d0), MephTot,         &
                dimTot, Aux3, dimTot, (0.d0,0.d0), Aux1, dimTot)

!   -- 2nd PART: 'G*(1st PART + H.c.)' --

!   ('Aux3 = GrTot * Aux1')
    call zgemm ('N', 'N', dimTot, dimTot, dimTot, (1.d0,0.d0), GrTot,   &
                dimTot, Aux1, dimTot, (0.d0,0.d0), Aux3, dimTot)

!   ('Aux3 = GrTot * Aux1^dagger + Aux3')
    call zgemm ('N', 'C', dimTot, dimTot, dimTot, (1.d0,0.d0), GrTot,   &
                dimTot, Aux1, dimTot, (1.d0,0.d0), Aux3, dimTot)

!   -- 3rd PART: 'G^dagger*Gamma_L*(2nd PART)' --

!   ('Aux2 = Gamma_LTot * Aux3')
    call zhemm ('L', 'L', dimTot, dimTot, (1.d0,0.d0), Gamma_LTot,      &
                dimTot, Aux3, dimTot, (0.d0,0.d0), Aux2, dimTot)

!   ('Aux1 = GrTot^dagger * Aux2')
    call zgemm ('C', 'N', dimTot, dimTot, dimTot, (1.d0,0.d0),          &
                GrTot, dimTot, Aux2, dimTot, (0.d0,0.d0), Aux1, dimTot)

    if (ephunit == ephIdx(ntypeunits+1)) then
       do j = idxF(ephunit),idxL(ephunit)
          do i = idxF(ephunit),idxL(ephunit)
             write (3331,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Aux1(i,j)),                                     &
                  DREAL(Asymm(i-idxF(ephunit)+1,j-idxF(ephunit)+1)),    &
                  DIMAG(Aux1(i,j)),                                     &
                  DIMAG(Asymm(i-idxF(ephunit)+1,j-idxF(ephunit)+1))
             write (3332,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Aux1(i,j)) -                                    &
                  DREAL(Asymm(i-idxF(ephunit)+1,j-idxF(ephunit)+1)),    &
                  DIMAG(Aux1(i,j)) -                                    &
                  DIMAG(Asymm(i-idxF(ephunit)+1,j-idxF(ephunit)+1))
          enddo
       enddo
    else if (ephunit == ephIdx(ntypeunits+2)) then
! OBS.: Only work if 'norbDyn' is equal to 'NR'!!
       dim = unitdimensions(ntypeunits+2)
       do j = idxF(ephunit),idxL(ephunit)
          do i = idxF(ephunit),idxL(ephunit)
             write (3331,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Aux1(dimTot-dim+i,dimTot-dim+j)),               &
                  DREAL(Asymm(i-idxF(ephunit)+1,j-idxF(ephunit)+1)) +   &
                  DREAL(Hc(i-idxF(ephunit)+1,j-idxF(ephunit)+1)),       &
                  DIMAG(Aux1(dimTot-dim+i,dimTot-dim+j)),               &
                  DIMAG(Asymm(i-idxF(ephunit)+1,j-idxF(ephunit)+1)) +   &
                  DIMAG(Hc(i-idxF(ephunit)+1,j-idxF(ephunit)+1))
             write (3332,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Aux1(dimTot-dim+i,dimTot-dim+j)) -              &
                  DREAL(Asymm(i-idxF(ephunit)+1,j-idxF(ephunit)+1)) -   &
                  DREAL(Hc(i-idxF(ephunit)+1,j-idxF(ephunit)+1)),       &
                  DIMAG(Aux1(dimTot-dim+i,dimTot-dim+j)) -              &
                  DIMAG(Asymm(i-idxF(ephunit)+1,j-idxF(ephunit)+1)) -   &
                  DIMAG(Hc(i-idxF(ephunit)+1,j-idxF(ephunit)+1))
          enddo
       enddo
    else
       dim = unitdimensions(ntypeunits+1)
       do k = 2,nunits+1
          utype = unit_type(k) ! current unit type
          if (ephunit == ephIdx(utype)) then
             do j = idxF(ephunit),idxL(ephunit)
                do i = idxF(ephunit),idxL(ephunit)
                   write (3331,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')     &
                      DREAL(Aux1(dim+i,dim+j)),                         &
                      DREAL(Asymm(i-idxF(ephunit)+1,j-idxF(ephunit)+1)),&
                      DIMAG(Aux1(dim+i,dim+j)),                         &
                      DIMAG(Asymm(i-idxF(ephunit)+1,j-idxF(ephunit)+1))
                   write (3332,'(e17.8e3,e17.8e3)')                     &
                      DREAL(Aux1(dim+i,dim+j)) -                        &
                      DREAL(Asymm(i-idxF(ephunit)+1,j-idxF(ephunit)+1)),&
                      DIMAG(Aux1(dim+i,dim+j)) -                        &
                      DIMAG(Asymm(i-idxF(ephunit)+1,j-idxF(ephunit)+1))
                enddo
             enddo
             EXIT
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


END MODULE idsrdr_current
