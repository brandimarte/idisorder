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
  use idsrdr_green,    only: 

  implicit none
  
  PUBLIC  :: current
  PRIVATE :: elastic, transmission, inelSymm


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
    integer :: i, j
    real(8) :: Tel, Tsym, Tasym
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8), allocatable, dimension (:,:) :: Gamma_L
    complex(8), allocatable, dimension (:,:) :: Gamma_R

!   Allocate coupling and auxiliary matrices.
    allocate (Gamma_L(NL,NL))
    allocate (Gamma_R(NR,NR))

! ATENCAO: POSSO USAR QUE Sigma_L/R SAO SIMETRICAS AQUI!!!

!   Sets the lead's coupling matricesx (triangular inferior part).
    if (NL >= NR) then
       Do j = 1,NR
          do i = j,NR
             Gamma_L(i,j) = zi * (Sigma_L(i,j) - DCONJG(Sigma_L(j,i)))
             Gamma_R(i,j) = zi * (Sigma_R(i,j) - DCONJG(Sigma_R(j,i)))
          enddo
          do i = NR+1,NL
             Gamma_L(i,j) = zi * (Sigma_L(i,j) - DCONJG(Sigma_L(j,i)))
          enddo
       Enddo
       Do j = NR+1,NL
          do i = j,NL
             Gamma_L(i,j) = zi * (Sigma_L(i,j) - DCONJG(Sigma_L(j,i)))
          enddo
       Enddo
    else ! NL < NR
       Do j = 1,NL
          do i = j,NL
             Gamma_L(i,j) = zi * (Sigma_L(i,j) - DCONJG(Sigma_L(j,i)))
             Gamma_R(i,j) = zi * (Sigma_R(i,j) - DCONJG(Sigma_R(j,i)))
          enddo
          do i = NL+1,NR
             Gamma_R(i,j) = zi * (Sigma_R(i,j) - DCONJG(Sigma_R(j,i)))
          enddo
       Enddo
       Do j = NL+1,NR
          do i = j,NR
             Gamma_R(i,j) = zi * (Sigma_R(i,j) - DCONJG(Sigma_R(j,i)))
          enddo
       Enddo
    endif

!   Compute elastic contribution.
    call elastic (Tel, NL, Gamma_L, NR, Gamma_R)

!   Compute symmetric part of inelastic contribution.
    call inelSymm (Tsym, ispin, NL, Gamma_L, NR, Gamma_R)

!   Compute asymmetric part of inelastic contribution.
    call inelAsymm (Tasym, ispin, NL, Gamma_L, NR, Gamma_R)

!   MPI_Reduce


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
!  ***************************** OUTPUT ******************************  !
!  real*8 Tsym             : Symmetric part of inelastic transmission   !
!  *******************************************************************  !
  subroutine inelSymm (Tsym, ispin, NL, Gamma_L, NR, Gamma_R)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_ephcoupl, only: neph, nModes, norbDyn, Meph
    use idsrdr_green,    only: Gr_Mn, Gr_1n, Gr_nn, Gr_1M

!   Input variables.
    integer, intent(in) :: ispin, NL, NR
    real(8), intent(out) :: Tsym
    complex(8), dimension (NL,NL), intent(in) :: Gamma_L
    complex(8), dimension (NR,NR), intent(in) :: Gamma_R

!   Local variables.
    integer :: j, w, i
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8), dimension(:,:), allocatable :: Aux1, Aux2, Aux3, Aux4,  &
                                               Aux5, Aux6, Aux7, GrCJG, A
    external :: zsymm, zhemm, zgemm

    if (IOnode) write (6,'(a)', advance='no')                           &
         '      computing inelastic current (symmetric part)... '

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

!         -- 2nd PART: 'G*Gamma_R*G^dagger*Meph*A*Meph' --

!         ('Aux1 = Aux2 * A') (where 'Aux2 = Gamma_R * Gr_Mn^* * Meph')
          call zhemm ('R', 'L', NR, norbDyn(j), (1.d0,0.d0),            &
                      A, norbDyn(j), Aux2, NR, (0.d0,0.d0), Aux1, NR)

!         ('Aux2 = Aux1 * Meph')
          call zsymm ('R', 'L', NR, norbDyn(j), (1.d0,0.d0),            &
                      Meph(j)%M(:,:,ispin,w), norbDyn(j),               &
                      Aux1, NR, (0.d0,0.d0), Aux2, NR)

!         -- 3rd PART: 'G^dagger*Gamma_L*(1st + 2nd PARTS)' --

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
          Tsym = 0.d0
          do i = 1,norbDyn(j)
             Tsym = Tsym + DREAL(Aux4(i,i))
          enddo
          do i = 1,NR
             Tsym = Tsym + DREAL(Aux7(i,i))
          enddo
          
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
!  ***************************** OUTPUT ******************************  !
!  real*8 Tasym            : Asymmetric part of inelastic transmission  !
!  *******************************************************************  !
  subroutine inelAsymm (Tasym, ispin, NL, Gamma_L, NR, Gamma_R)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_ephcoupl, only: neph, nModes, norbDyn, Meph
    use idsrdr_green,    only: Gr_Mn, Gr_1n, Gr_1M

!   Input variables.
    integer, intent(in) :: ispin, NL, NR
    real(8), intent(out) :: Tasym
    complex(8), dimension (NL,NL), intent(in) :: Gamma_L
    complex(8), dimension (NR,NR), intent(in) :: Gamma_R

!   Local variables.
    integer :: j, w, i !, k, l
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8), dimension(:,:), allocatable :: Aux1, Aux2, Aux3,        &
                                               Aux4, Aux5, Aux6,        &
                                               Gr_MnCJG, Gr_1nCJG
    external :: zsymm, zhemm, zgemm

    if (IOnode) write (6,'(a)', advance='no')                           &
         '      computing inelastic current (asymmetric part)... '

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

!         ('Aux4 = Gr_Mn^* * Meph')
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

!         ('Aux3 = Gr_1n^T * Aux5')
          call zgemm ('T', 'N', norbDyn(j), norbDyn(j), NL,             &
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
          Tasym = 0.d0
          do i = 1,norbDyn(j)
             Tasym = Tasym + DREAL(Aux3(i,i))
          enddo
          do i = 1,NR
             Tasym = Tasym + DREAL(Aux6(i,i))
          enddo
          
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


END MODULE idsrdr_current
