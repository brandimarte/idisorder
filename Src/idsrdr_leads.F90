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
!                          MODULE idsrdr_leads                          !
!  *******************************************************************  !
!  Description: read the dimensions of electrodes, the principal layer  !
!  (PL) hamiltonians and overlaps, and the coupling hamiltonians and    !
!  overlaps between PLs.                                                !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *******************************************************************  !

MODULE idsrdr_leads

!
!   Modules
!
  use parallel,        only: 
  use idsrdr_options,  only: 
  use idsrdr_check,    only: 
  use idsrdr_io,       only: 
  use idsrdr_string,   only: 

  implicit none
  
  PUBLIC  :: NL, NR, numberL, numberR, EfLead, S0_L, S0_R, S1_L, S1_R,  &
             H0_L, H0_R, H1_L, H1_R, VHL, VHR, VSL, VSR, QL, QR,        &
             Sigma_L, Sigma_R, side_rankL, side_rankR,                  &
             readleads, leadsSelfEn, freeleads
  PRIVATE :: selfenergy

  integer :: NL ! Number of left lead orbitals
  integer :: NR ! Number of right lead orbitals
  integer, allocatable, dimension (:) :: numberL !
  integer, allocatable, dimension (:) :: numberR !

  real(8) :: EfLead ! Lead Fermi energy

  complex(8), allocatable, dimension (:,:) :: S0_L ! left PL overlap
  complex(8), allocatable, dimension (:,:) :: S0_R ! right PL overlap
  complex(8), allocatable, dimension (:,:) :: S1_L ! left coupl. overlap
  complex(8), allocatable, dimension (:,:) :: S1_R ! right coupl. overlap
  complex(8), allocatable, dimension (:,:,:) :: H0_L ! left PL
  complex(8), allocatable, dimension (:,:,:) :: H0_R ! right PL
  complex(8), allocatable, dimension (:,:,:) :: H1_L ! left coupl. 
  complex(8), allocatable, dimension (:,:,:) :: H1_R ! right coupl.
  complex(8), allocatable, dimension (:,:,:) :: VHL !
  complex(8), allocatable, dimension (:,:,:) :: VHR !
  complex(8), allocatable, dimension (:,:,:) :: VSL !
  complex(8), allocatable, dimension (:,:,:) :: VSR !
  complex(8), allocatable, dimension (:,:,:) :: QL !
  complex(8), allocatable, dimension (:,:,:) :: QR !
  complex(8), allocatable, dimension (:,:) :: Sigma_L ! Left self-energy
  complex(8), allocatable, dimension (:,:) :: Sigma_R ! Right self-energy

  character, dimension (:), allocatable :: side_rankL ! 
  character, dimension (:), allocatable :: side_rankR ! 


CONTAINS


!  *******************************************************************  !
!                               readleads                               !
!  *******************************************************************  !
!  Description: subroutine to read the dimensions of electrodes and     !
!  principal layer (PL) hamiltonians from 'bulklft.DAT' and             !
!  'bulkrgt.DAT' files. Also build and store the overlap matrix, the    !
!  coupling hamiltonians and overlaps between PLs.                      !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode               : True if it is the I/O node            !
!  integer nspin                : Number of spin components             !
!  real*8 temp                  : Electronic temperature                !
!  character(60) directory      : Working directory                     !
!  logical tightbinding         : Tight-binding calculation?            !
!  real*8 TBeFermi              : Tight-binding Fermi energy            !
!  ****************************** INPUT ******************************  !
!  integer nsc(2)               : Number of unit cells along parallel   !
!                                 directions                            !
!  ***************************** OUTPUT ******************************  !
!  integer NL                     : Number of left lead orbitals        !
!  integer NR                     : Number of right lead orbitals       !
!  real*8 EfLead                  : Lead Fermi energy                   !
!  integer numberL(nspin)         :                                     !
!  integer numberR(nspin)         :                                     !
!  complex(8) S0_L(NL,NL)         : Left lead PL overlap                !
!  complex(8) S0_R(NR,NR)         : Right lead PL overlap               !
!  complex(8) S1_L(NL,NL)         : Left lead coupling overlap          !
!  complex(8) S1_R(NR,NR)         : Right lead coupling overlap         !
!  complex(8) H0_L(NL,NL,nspin)   : Left lead PL hamiltonian            !
!  complex(8) H0_R(NR,NR,nspin)   : Right lead PL hamiltonian           !
!  complex(8) H1_L(NL,NL,nspin)   : Left lead coupling                  !
!  complex(8) H1_R(NR,NR,nspin)   : Right lead PL coupling              !
!  complex(8) VHL(NL,NL,nspin)    :                                     !
!  complex(8) VHR(NR,NR,nspin)    :                                     !
!  complex(8) VSL(NL,NL,nspin)    :                                     !
!  complex(8) VSR(NR,NR,nspin)    :                                     !
!  complex(8) QL(NL,NL,nspin)     :                                     !
!  complex(8) QR(NR,NR,nspin)     :                                     !
!  character(1) side_rankL(nspin) :                                     !
!  character(1) side_rankR(nspin) :                                     !
!  *******************************************************************  !
  subroutine readleads (nsc)

!
!   Modules
!
#ifdef MPI
    use parallel,        only: IOnode, MPI_Comm_MyWorld
#else
    use parallel,        only: IOnode
#endif
    use idsrdr_options,  only: nspin, temp, directory,                  &
                               tightbinding, TBeFermi
    use idsrdr_io,       only: IOassign, IOclose
    use idsrdr_string,   only: STRpaste

#ifdef MPI
    include "mpif.h"
#endif

!   Input variables.
    integer, intent(in) :: nsc(2)

!   Local variables.
    integer :: iu, nspinu, maxnh
    real(8) :: EfLeadR
    character (len=30) :: slabeli
    character(len=72) :: file
    external :: zhsunits, ranksvd
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

    if (IOnode) then

       write (6,'(/,30("*"),a,31("*"))')                                &
            ' Leads input data '

       call IOassign (iu)
       call STRpaste (directory, 'bulklft.DAT', file)
       open (iu, file=file, status='old')
       read (iu,*) slabeli, NL, nspinu, maxnh, EfLead
       write (6,2)                                                      &
            'readleads: Left lead system label                  ' //    &
            '      =  ', slabeli
       write (6,4)                                                      &
            'readleads: Number of left lead orbitals            ' //    &
            '      =', NL
       if (nspinu /= nspin) then
          write (6,'(a)') 'WARNING: spin components from left lead ' // &
               'differs from input option!'
       endif
       call IOclose (iu)
       call IOassign (iu)
       call STRpaste (directory, 'bulkrgt.DAT', file)
       open (iu, file=file, status='old')
       read (iu,*) slabeli, NR, nspinu, maxnh, EfLeadR
       if (nspinu /= nspin) then
          write (6,'(a)') 'WARNING: spin components from left lead ' // &
               'differs from input option! (using the input option...)'
       endif
       write (6,2)                                                      &
            'readleads: Right lead system label                 ' //    &
            '      =  ', slabeli
       write (6,4)                                                      &
            'readleads: Number of right lead orbitals           ' //    &
            '      =', NR
       if (EfLead /= EfLeadR) then
          write (6,'(a)') 'WARNING: the Fermi energy from leads ' //    &
               'differs! (using the left lead Fermi energy...)'
       endif
       if (tightbinding) then
          EfLead = TBeFermi
       endif
       write (6,6)                                                      &
            'readleads: Lead Fermi energy                       ' //    &
            '      =', EfLead, ' Ry'
       call IOclose (iu)

       write (6,'(2a,/)') 'readleads: ', repeat('*', 68)

    endif

#ifdef MPI
    call MPI_Bcast (NL, 1, MPI_Integer, 0, MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (NR, 1, MPI_Integer, 0, MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (EfLead, 1 ,MPI_Double_Precision, 0,                 &
                    MPI_Comm_MyWorld, MPIerror)
#endif

!   Allocate hamiltonians and overlap matrices.
    allocate (S0_L(NL,NL), S1_L(NL,NL))
    allocate (H0_L(NL,NL,nspin), H1_L(NL,NL,nspin))
    allocate (S0_R(NR,NR), S1_R(NR,NR))
    allocate (H0_R(NR,NR,nspin), H1_R(NR,NR,nspin))

!   Reads and calculates the Hamiltonians
!   and overlaps matrices of the leads.
    if (IOnode) then
       If (tightbinding) Then

!         Creates overlap and hamiltonians for tight-binding calculation.
          call TBreadleads

       Else
          call STRpaste (directory, 'bulklft.DAT', file)
          call zhsunits (nspin, nspin, NL, nsc, 1, 0, 0, .true.,        &
                         1, temp, H0_L, H1_L, S0_L, S1_L, file)
          call STRpaste (directory, 'bulkrgt.DAT', file)
          call zhsunits (nspin, nspin, NR, nsc, 1, 0, 0, .true.,        &
                         1, temp, H0_R, H1_R, S0_R, S1_R, file)
       EndIf
    endif

#ifdef MPI
    call MPI_Bcast (H0_L, NL*NL*nspin, MPI_Double_Complex, 0,           &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (S0_L, NL*NL, MPI_Double_Complex, 0,                 &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (H1_L, NL*NL*nspin, MPI_Double_Complex, 0,           &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (S1_L, NL*NL, MPI_Double_Complex, 0,                 &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (H0_R, NR*NR*nspin, MPI_Double_Complex, 0,           &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (S0_R, NR*NR, MPI_Double_Complex, 0,                 &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (H1_R, NR*NR*nspin, MPI_Double_Complex, 0,           &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (S1_R, NR*NR, MPI_Double_Complex, 0,                 &
                    MPI_Comm_MyWorld, MPIerror)
#endif


!   Allocate matrices and ranks for leads decimation.
    allocate (numberL(nspin), numberR(nspin))
    allocate (side_rankL(nspin), side_rankR(nspin))
    allocate (QL(NL,NL,nspin), VHL(NL,NL,nspin), VSL(NL,NL,nspin))
    allocate (QR(NR,NR,nspin), VHR(NR,NR,nspin), VSR(NR,NR,nspin))

!   Determine the direction of decimation and
!   perform a Singular Value Decomposiation (SVD).
    if (IOnode) write (6,'(/,a)')                                       &
         'ranksvd: Computing leads decimation...'
    call ranksvd (NL, nspin, H1_L, S1_L, side_rankL,                    &
                  numberL, VHL, VSL, QL)
    call ranksvd (NR, nspin, H1_R, S1_R, side_rankR,                    &
                  numberR, VHR, VSR, QR)
    if (IOnode) write(6,'(/,a)') 'ranksvd: done!'

!   Allocate lead self-energies matrices.
    allocate (Sigma_L(NL,NL), Sigma_R(NR,NR))

2   format(a,a)
4   format(a,i7)
6   format(a,f12.4,a)


  end subroutine readleads


!  *******************************************************************  !
!                              TBreadleads                              !
!  *******************************************************************  !
!  Description: Creates the lead's overlap and hamiltonian matrices     !
!  for tight-binding calculation.                                       !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer nspin                  : Number of spin components           !
!  real*8 TBenerg                 : Tight-binding site energy           !
!  real*8 TBcoupl                 : Tight-binding site couplings        !
!  ***************************** OUTPUT ******************************  !
!  complex(8) S0_L(NL,NL)         : Left lead PL overlap                !
!  complex(8) S0_R(NR,NR)         : Right lead PL overlap               !
!  complex(8) S1_L(NL,NL)         : Left lead coupling overlap          !
!  complex(8) S1_R(NR,NR)         : Right lead coupling overlap         !
!  complex(8) H0_L(NL,NL,nspin)   : Left lead PL hamiltonian            !
!  complex(8) H0_R(NR,NR,nspin)   : Right lead PL hamiltonian           !
!  complex(8) H1_L(NL,NL,nspin)   : Left lead coupling                  !
!  complex(8) H1_R(NR,NR,nspin)   : Right lead PL coupling              !
!  *******************************************************************  !
  subroutine TBreadleads

!
!   Modules
!
    use idsrdr_options,  only: nspin, TBenerg, TBcoupl

!   Local variables.
    integer :: i, j, s

    S0_L = 0.d0
    S1_L = 0.d0
    H0_L = 0.d0
    H1_L = 0.d0
    S0_R = 0.d0
    S1_R = 0.d0
    H0_R = 0.d0
    H1_R = 0.d0

    do i = 1,NL-1
       H0_L(i,i,1) = DCMPLX(TBenerg)
       H0_L(i,i+1,1) = DCMPLX(TBcoupl)
       H0_L(i+1,i,1) = DCMPLX(TBcoupl)
       S0_L(i,i) = (1.d0,0.d0)
    enddo
    H0_L(i,i,1) = DCMPLX(TBenerg)
    S0_L(i,i) = (1.d0,0.d0)
    H1_L(NL,1,1) = DCMPLX(TBcoupl)

    do i = 1,NR-1
       H0_R(i,i,1) = DCMPLX(TBenerg)
       H0_R(i,i+1,1) = DCMPLX(TBcoupl)
       H0_R(i+1,i,1) = DCMPLX(TBcoupl)
       S0_R(i,i) = (1.d0,0.d0)
    enddo
    H0_R(i,i,1) = DCMPLX(TBenerg)
    S0_R(i,i) = (1.d0,0.d0)
    H1_R(NR,1,1) = DCMPLX(TBcoupl)

    do s = 2,nspin ! over other spin components

       do i = 1,NL-1
          H0_L(i,i,s) = DCMPLX(TBenerg)
          H0_L(i,i+1,s) = DCMPLX(TBcoupl)
          H0_L(i+1,i,s) = DCMPLX(TBcoupl)
       enddo
       H0_L(i,i,s) = DCMPLX(TBenerg)
       H1_L(NL,1,s) = DCMPLX(TBcoupl)

       do i = 1,NR-1
          H0_R(i,i,s) = DCMPLX(TBenerg)
          H0_R(i,i+1,s) = DCMPLX(TBcoupl)
          H0_R(i+1,i,s) = DCMPLX(TBcoupl)
       enddo
       H0_R(i,i,s) = DCMPLX(TBenerg)
       H1_R(NR,1,s) = DCMPLX(TBcoupl)

    enddo

!   Write matrices on screen.
    do s = 1,nspin
       write (6, '(a)') 'H0_L'
       do i = 1,NL
          do j = 1,NL-1
             write (6, '(e17.8e3)', advance='no') DREAL(H0_L(i,j,s))
          enddo
          write (6, '(e17.8e3)') DREAL(H0_L(i,j,s))
       enddo

       write (6, '(a)') 'H1_L'
       do i = 1,NL
          do j = 1,NL-1
             write (6, '(e17.8e3)', advance='no') DREAL(H1_L(i,j,s))
          enddo
          write (6, '(e17.8e3)') DREAL(H1_L(i,j,s))
       enddo

       write (6, '(a)') 'S0_L'
       do i = 1,NL
          do j = 1,NL-1
             write (6, '(e17.8e3)', advance='no') DREAL(S0_L(i,j))
          enddo
          write (6, '(e17.8e3)') DREAL(S0_L(i,j))
       enddo

       write (6, '(a)') 'S1_L'
       do i = 1,NL
          do j = 1,NL-1
             write (6, '(e17.8e3)', advance='no') DREAL(S1_L(i,j))
          enddo
          write (6, '(e17.8e3)') DREAL(S1_L(i,j))
       enddo

       write (6, '(a)') 'H0_R'
       do i = 1,NR
          do j = 1,NR-1
             write (6, '(e17.8e3)', advance='no') DREAL(H0_R(i,j,s))
          enddo
          write (6, '(e17.8e3)') DREAL(H0_R(i,j,s))
       enddo

       write (6, '(a)') 'H1_R'
       do i = 1,NR
          do j = 1,NR-1
             write (6, '(e17.8e3)', advance='no') DREAL(H1_R(i,j,s))
          enddo
          write (6, '(e17.8e3)') DREAL(H1_R(i,j,s))
       enddo

       write (6, '(a)') 'S0_R'
       do i = 1,NR
          do j = 1,NR-1
             write (6, '(e17.8e3)', advance='no') DREAL(S0_R(i,j))
          enddo
          write (6, '(e17.8e3)') DREAL(S0_R(i,j))
       enddo

       write (6, '(a)') 'S1_R'
       do i = 1,NR
          do j = 1,NR-1
             write (6, '(e17.8e3)', advance='no') DREAL(S1_R(i,j))
          enddo
          write (6, '(e17.8e3)') DREAL(S1_R(i,j))
       enddo

    enddo


  end subroutine TBreadleads


!  *******************************************************************  !
!                              leadsSelfEn                              !
!  *******************************************************************  !
!  Description: computes the leads self-energies for a given energy.    !
!                                                                       !
!  Use 'SELFENERGY' subroutine of A.R.Rocha (Smeagol - 2003).           !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 Ei                            : Energy grid point             !
!  ***************************** OUTPUT ******************************  !
!  complex(8) Sigma_L(NL,NL)      : Left-lead self-energy               !
!  complex(8) Sigma_R(NR,NR)      : Right-lead self-energy              !
!  *******************************************************************  !
  subroutine leadsSelfEn (Ei, ispin)

!
!   Modules
!
    use parallel,        only: IOnode

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    integer :: INFO, NCHAN

    if (IOnode) then
       write (6,'(/,a,f8.5,a,i1,a,/)') '[energy(Ry),spin] = [',         &
            Ei, ',', ispin, ']'
       write (6,'(a)', advance='no')                                    &
            '      computing leads self-energies... '
    endif

    call selfenergy ('L', side_rankL(ispin), numberL(ispin), NL,        &
                     DCMPLX(Ei), H0_L(:,:,ispin), VHL(:,:,ispin),       &
                     S0_L, VSL(:,:,ispin), Sigma_L, QL(:,:,ispin),      &
                     INFO, NCHAN)

    if (INFO == 0) then
       call selfenergy ('R', side_rankR(ispin), numberR(ispin), NR,     &
                        DCMPLX(Ei), H0_R(:,:,ispin), VHR(:,:,ispin),    &
                        S0_R, VSR(:,:,ispin), Sigma_R, QR(:,:,ispin),   &
                        INFO, NCHAN)
    endif

    if (INFO == 0 .and. IOnode) then
       write(6,'(a)') ' ok!'
       write(6,'(a,i3)') '        -> number of open channels = ', NCHAN
    endif


  end subroutine leadsSelfEn


!  *******************************************************************  !
!                              selfenergy                               !
!  *******************************************************************  !
!  Description: calculates the self-energies of the right and left      !
!  hand-side leads.                                                     !
!                                                                       !
!  From SMEAGOL code - 2003.                                            !
!                                                                       !
!  Written by Alexandre Reily Rocha, Jun 2003.                          !
!  Computational Spintronics Group                                      !
!  Trinity College Dublin                                               !
!  e-mail: reilya@ift.unesp.br                                          !
!  ***************************** HISTORY *****************************  !
!  Original version:    June 2003                                       !
!  ****************************** INPUT ******************************  !
!  character(len=1) SIDE      : 'L' (left) or 'R' (right) depending on  !
!                               the self-energy calculated              !
!  character(len=1) side_rank : '0' or 'N' or 'C' depending on the      !
!                               direction of decimation                 !
!  integer number             : Number of decimated states              !
!  integer N2                 : Number of basis orbitals on the lead    !
!  complex*8 Ei               : Energy                                  !
!  complex*8 H0(NL,NL)        : Hamiltonian of the lead lead            !
!  complex*8 H1(N2,N2)        : Coupling Matrix                         !
!  complex*8 S0(N2,N2)        : On-site Overlap Matrix                  !
!  complex*8 S1(N2,N2)        : Neighbour Overlap Matrix                !
!  complex*8 Q(N2,N2)         : Decimation matrix                       !
!  ***************************** OUTPUT ******************************  !
!  complex*8 Sigma            : Calculated self-energy                  !
!  integer INFO               : If INFO=0, self-energy calculated       !
!                               with success, NOT 0 otherwise           !
!  integer nrchan             : Number of open channels                 !
!  *******************************************************************  !
  subroutine selfenergy (SIDE, side_rank, number, N2, Ei, H0, H1,       &
                         S0, S1, Sigma, Q, INFO, nrchan)

!
!   Modules
!
    use idsrdr_check,    only: CHECKzgetrf, CHECKzgetri

!   Input variables.
    integer, intent(in) :: number, N2
    integer, intent(out) :: INFO, nrchan
    complex(8), intent(in) :: Ei
    complex(8), dimension (N2,N2), intent(in) :: H0, H1, S0, S1, Q
    complex(8), dimension (N2,N2), intent(out) :: Sigma
    character(len=1), intent(in) :: SIDE, side_rank

!   Local variables.
    integer :: N3, nlchan, I, J
    integer, allocatable, dimension (:) :: IPIV, IPIV2
    complex(8), parameter :: smear = (0.d0,0.d0)
    complex(8), allocatable, dimension (:,:) :: T1_aux, T1_dag, T1,     &
                                                H0_aux, H1_aux,         &
                                                H1_dag_aux, h0corr,     &
                                                Gr_2, Gr_square
    complex(8), allocatable, dimension (:,:) :: Gr
    complex(8), allocatable, dimension (:,:) :: Gr_1
    complex(8), allocatable, dimension (:,:) :: foo23
    complex(8), allocatable, dimension (:,:) :: foo32, aux32
    external :: zgemm, DECIMATE_LEADS, LEADS

!   Allocate matrices and arrays.
    allocate (T1_aux(N2,N2), H0_aux(N2,N2), H1_aux(N2,N2),              &
              H1_dag_aux(N2,N2), h0corr(N2,N2))
    allocate (Gr(2*(N2-number),2*(N2-number)))
    allocate (Gr_1(N2-number,N2-number))
    allocate (IPIV(n2-number))
    allocate (IPIV2(n2))

!   Initialize variables.
    N3 = N2 - number
    H0_aux = H0 - Ei*S0
    Gr = (0.d0,0.d0)
    nrchan = 0
    nlchan = 0

!   ('T1_aux = Q^dagger * H0_aux')
    call zgemm ('C', 'N', N2, N2, N2, (1.d0,0.d0), Q, N2,               &
                H0_aux, N2, (0.d0,0.d0), T1_aux, N2)
!   ('H0_aux = T1_aux * Q')
    call zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), T1_aux, N2,          &
                Q, N2, (0.d0,0.d0), H0_aux, N2)

    IF (side_rank /= '0') THEN

!      Compute 'H1_dag_aux' according to the direction of decimation.
       If (side_rank == 'N') Then

          T1_aux = H1 - Ei*S1

!         ('H1_aux = Q^dagger * T1_aux')
          call zgemm ('C', 'N', N2, N2, N2, (1.d0,0.d0), Q, N2,         &
                      T1_aux, N2, (0.d0,0.d0), H1_aux, N2)

!         ('T1_aux = H1^dagger - Ei*S1^dagger')
          do I = 1,N2
             do J = 1,N2
                T1_aux(I,J) = DCONJG(H1(J,I)) - Ei*DCONJG(S1(J,I))
             enddo
          enddo

!         ('H1_dag_aux = T1_aux * Q')
          call zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), T1_aux, N2,    &
                      Q, N2, (0.d0,0.d0), H1_dag_aux, N2)

       Else

          T1_aux = H1 - Ei*S1

!         ('H1_aux = T1_aux * Q')
          call zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), T1_aux, N2,    &
                      Q, N2, (0.d0,0.d0), H1_aux, N2)

!         ('T1_aux = H1^dagger - Ei*S1^dagger')
          do I = 1,N2
             do J = 1,N2
                T1_aux(I,J) = DCONJG(H1(J,I)) - Ei*DCONJG(S1(J,I))
             enddo
          enddo

!         ('H1_dag_aux = Q^dagger * T1_aux')
          call zgemm ('C', 'N', N2, N2, N2, (1.d0,0.d0), Q, N2,         &
                      T1_aux, N2, (0.d0,0.d0), H1_dag_aux, N2)

       EndIf

!      Calculates the decimated Hamiltonian and the correction.
       call DECIMATE_LEADS (side_rank, SIDE, N2, number, H0_aux,        &
                            H1_aux, H1_dag_aux, h0corr)

    ELSE

       H0_aux = H0 - Ei*S0
       H1_aux = H1 - Ei*S1
       h0corr = 0.d0
       do I = 1,N2
          do J = 1,N2
             H1_dag_aux(I,J) = DCONJG(H1(J,I)) - Ei*DCONJG(S1(J,I))
          enddo
       enddo

    ENDIF ! IF (side_rank /= '0')

!   Calculates the surface Green's function.
    call LEADS (SIDE, N3, 2*N3, 4*N3, Ei,                               &
                H0_aux(number+1:N2,number+1:N2),                        &
                H1_aux(number+1:N2,number+1:N2),                        &
                H1_dag_aux(number+1:N2,number+1:N2),                    &
                Gr, nrchan, nlchan, INFO)

!   Free memory.
    deallocate (H0_aux)

!   Allocate memory.
    allocate (T1(N2,N2), T1_dag(N2,N2))
    allocate (foo23(N2,N2-number))
    allocate (foo32(N2-number,N2), aux32(N2-number,N2))

    IF (info == 0) THEN
       If (side == 'L') Then
          if (side_rank == 'N') then

!            Allocate matrix.
             allocate (Gr_2(N2,N2))

             T1(1:N3,:) = H1_aux(number+1:N2,:)
             T1_dag(:,1:N3) = H1_dag_aux(:,number+1:N2)
             Gr_1 = Gr(1:N3,1:N3)

!            ('T1_aux = Gr_1 * T1')
             foo32 = T1(1:N3,:)
             call zgemm ('N', 'N', N3, N2, N3,( 1.d0,0.d0), Gr_1, N3,   &
                         foo32, N3, (0.d0,0.d0), aux32, N3)

!            ('Sigma = T1_dag * T1_aux')
             foo23 = T1_dag(:,1:N3)
             call zgemm ('N', 'N', N2, N2, N3, (1.d0,0.d0), foo23, N2,  &
                         aux32, N3, (0.d0,0.d0), Sigma, N2)

!            (Gr_2 = (-h0corr - Sigma)^-1')
             Gr_2 = -h0corr - Sigma
             call CHECKzgetrf (N2, Gr_2, IPIV2)
             call CHECKzgetri (N2, Gr_2, IPIV2)

             if (DIMAG(smear) > 1.d-7) then ! PB: it will never happen
                                            ! since 'smear'is set to 0?

!               Allocate matrix.
                allocate (Gr_square(N2,N2))

!               ('Gr_square = Gr_2 * Gr_2')
                call zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), Gr_2,    &
                            N2, Gr_2, N2, (0.d0,0.d0), Gr_square, N2)

                Gr_2 = Gr_2 - smear*Gr_square

!               Free memory.
                deallocate (Gr_square)

             endif

             T1 = H1 - Ei*S1

!            ('T1_aux = T1 * Q^dagger')
             call zgemm ('N', 'C', N2, N2, N2, (1.d0,0.d0), T1, N2,     &
                         Q, N2, (0.d0,0.d0), T1_aux, N2)
!            ('T1 = Q^dagger * T1_aux')
             call zgemm ('C', 'N', N2, N2, N2, (1.d0,0.d0), Q, N2,      &
                         T1_aux, N2, (0.d0,0.d0), T1, N2)     

!            ('T1_dag = H1^dagger - Ei*S1^dagger')
             do i = 1,N2
                do j = 1,N2
                   T1_dag(I,J) = DCONJG(H1(J,I)) - Ei*DCONJG(S1(J,I))
                enddo
             enddo

!            ('T1_aux = T1_dag * Q')
             call zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), T1_dag, N2, &
                         Q, N2, (0.d0,0.d0), T1_aux, N2)
!            ('T1_dag = Q * T1_aux')
             call zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), Q, N2,      &
                         T1_aux, N2, (0.d0,0.d0), T1_dag, N2)     

!            ('T1_aux = Gr_2 * T1')
             call zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), Gr_2, N2,   &
                         T1, N2, (0.d0,0.d0), T1_aux, N2)

!            ('Sigma = T1_dag * T1_aux')
             call zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), T1_dag, N2, &
                         T1_aux, N2, (0.d0,0.d0), Sigma, N2)

!            Free memory.
             deallocate (Gr_2)

          else ! side_rank /= 'N'

             T1(1:N3,:) = H1(number+1:N2,:) - Ei*S1(number+1:N2,:)

!            ('T1_dag = H1^dagger - Ei*S1^dagger')
             do I = 1,N2
                do J = 1,N3
                   T1_dag(I,J) = DCONJG(H1(number+J,I))                 &
                        - Ei*DCONJG(S1(number+J,I))
                enddo
             enddo
        
             Gr_1 = Gr(1:N3,1:N3)

             if (side_rank /= '0') then

!               (Gr_1 = (Gr(1:N3,1:N3))^-1')
                call CHECKzgetrf (N3, Gr_1, IPIV)
                call CHECKzgetri (N3, Gr_1, IPIV)

!               (Gr_1 = (Gr_1 - h0corr)^-1')
                Gr_1 = Gr_1 - h0corr(number+1:N2,number+1:N2)
                call CHECKzgetrf (N3, Gr_1, IPIV)
                call CHECKzgetri (N3, Gr_1, IPIV)

             endif

             if (DIMAG(smear) > 1.d-7) then ! PB: it will never happen
                                            ! since 'smear'is set to 0?

!               Allocate matrix.
                allocate (Gr_square(N3,N3))

!               ('Gr_square = Gr_1 * Gr_1')
                call zgemm ('N', 'N', N3, N3, N3, (1.d0,0.d0),          &
                            Gr_1, N3, Gr_1,N3, (0.d0,0.d0),             &
                            Gr_square, N3)

                Gr_1 = Gr_1 - smear*Gr_square

!               Free memory.
                deallocate (Gr_square)

             endif

!            ('T1_aux = Gr_1 * T1')
             foo32 = T1(1:N3,:)
             call zgemm ('N', 'N', N3, N2, N3, (1.d0,0.d0), Gr_1, N3,   &
                         foo32, N3, (0.d0,0.d0), aux32, N3)

!            ('Sigma = T1_dag * T1_aux')
             foo23 = T1_dag(:,1:N3)
             call zgemm ('N', 'N', N2, N2, N3, (1.d0,0.d0), foo23, N2,  &
                         aux32, N3, (0.d0,0.d0), Sigma, N2)
 
          endif ! if (side_rank == 'N')
        
       Else ! side == 'R'

          if (side_rank == 'C') then

!            Allocate matrix.
             allocate (Gr_2(N2,N2))

             T1(:,1:N3) = H1_aux(:,number+1:N2)
             T1_dag(1:N3,:) = H1_dag_aux(number+1:N2,:)
             Gr_1 = Gr(N3+1:2*N3,N3+1:2*N3)
          
!            ('T1_aux = Gr_1 * T1_dag')
             foo32 = T1_dag(1:N3,:)
             call zgemm ('N', 'N', N3, N2, N3,( 1.d0,0.d0), Gr_1, N3,   &
                         foo32, N3, (0.d0,0.d0), aux32, N3)

!            ('Sigma = T1 * T1_aux')
             foo23 = T1(:,1:N3)
             call zgemm ('N', 'N', N2, N2, N3, (1.d0,0.d0), foo23, N2,  &
                         aux32, N3, (0.d0,0.d0), Sigma, N2)

!            (Gr_2 = (-h0corr - Sigma)^-1')
             Gr_2 = -h0corr - Sigma
             call CHECKzgetrf (N2, Gr_2, IPIV2)
             call CHECKzgetri (N2, Gr_2, IPIV2)

             if (DIMAG(smear) > 1.d-7) then ! PB: it will never happen
                                            ! since 'smear'is set to 0?

!               Allocate matrix.
                allocate (Gr_square(N2,N2))

!               ('Gr_square = Gr_2 * Gr_2')
                call zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), Gr_2,    &
                            N2, Gr_2, N2, (0.d0,0.d0), Gr_square, N2)

                Gr_2 = Gr_2 - smear*Gr_square

!               Free memory.
                deallocate (Gr_square)

             endif
                                               
             T1 = H1 - Ei*S1

!            ('T1_aux = T1 * Q')
             call zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), T1, N2,     &
                         Q, N2, (0.d0,0.d0), T1_aux, N2)
!            ('T1 = Q * T1_aux')
             call zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), Q, N2,      &
                         T1_aux, N2, (0.d0,0.d0), T1, N2)     

!            ('T1_dag = H1^dagger - Ei*S1^dagger')
             do I = 1,N2
                do J = 1,N2
                   T1_dag(I,J) =  H1(J,I) - Ei*S1(J,I)
                enddo
             enddo

!            ('T1_aux = T1_dag * Q^dagger')
             call zgemm ('N', 'C', N2, N2, N2, (1.d0,0.d0), T1_dag, N2, &
                         Q, N2, (0.d0,0.d0), T1_aux, N2)
!            ('T1_dag = Q^dagger * T1_aux')
             call zgemm ('C', 'N', N2, N2, N2, (1.d0,0.d0), Q, N2,      &
                         T1_aux, N2, (0.d0,0.d0), T1_dag, N2)     

!            ('T1_aux = Gr_2 * T1_dag')
             call zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), Gr_2, N2,   &
                         T1_dag, N2, (0.d0,0.d0), T1_aux, N2)

!            ('Sigma = T1 * T1_aux')
             call zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), T1, N2,     &
                         T1_aux, N2, (0.d0,0.d0), Sigma, N2)

!            Free memory.
             deallocate (Gr_2)

          else ! side_rank /= 'C'

             T1(:,1:N3) = H1(:,number+1:N2) - Ei*S1(:,number+1:N2)

!            ('T1_dag = H1^dagger - Ei*S1^dagger')
             do I = 1,N3
                do J = 1,N2
                   T1_dag(I,J) = DCONJG(H1(J,number+I))                 &
                        - Ei*DCONJG(S1(J,number+I))
                enddo
             enddo

             Gr_1 = Gr(N3+1:2*N3,N3+1:2*N3)
          
             if (side_rank /= '0') then

!               (Gr_1 = (Gr(1:N3,1:N3))^-1')
                call CHECKzgetrf (N3, Gr_1, IPIV)
                call CHECKzgetri (N3, Gr_1, IPIV)

!               (Gr_1 = (Gr_1 - h0corr)^-1')
                Gr_1 = Gr_1 - h0corr(number+1:N2,number+1:N2)
                call CHECKzgetrf (N3, Gr_1, IPIV)
                call CHECKzgetri (N3, Gr_1, IPIV)

             endif

             iF (DIMAG(smear) > 1.d-7) then ! PB: it will never happen
                                            ! since 'smear'is set to 0?

!               Allocate matrix.
                allocate (Gr_square(N3,N3))

!               ('Gr_square = Gr_1 * Gr_1')
                call zgemm ('N', 'N', N3, N3, N3, (1.d0,0.d0),          &
                            Gr_1, N3, Gr_1, N3, (0.d0,0.d0),            &
                            Gr_square, N3)

                Gr_1 = Gr_1 - smear*Gr_square

!               Free memory.
                deallocate (Gr_square)

             endif

!            ('T1_aux = Gr_1 * T1_dag')
             foo32 = T1_dag(1:N3,:)
             call zgemm ('N', 'N', N3, N2, N3, (1.d0,0.d0), Gr_1, N3,   &
                         foo32, N3, (0.d0,0.d0), aux32, N3)

!            ('Sigma = T1 * T1_aux')
             foo23 = T1(:,1:N3)
             call zgemm ('N', 'N', N2, N2, N3, (1.d0,0.d0), foo23, N2,  &
                         aux32, N3, (0.d0,0.d0), Sigma, N2)

          endif ! if (side_rank == 'C')

       EndIf ! If (side == 'L')

    ELSE

       RETURN

    ENDIF ! IF (info == 0)

!   Free memory.
    deallocate (T1_aux, H1_aux, H1_dag_aux, h0corr)
    deallocate (Gr)
    deallocate (Gr_1)
    deallocate (IPIV)
    deallocate (IPIV2)
    deallocate (foo23)
    deallocate (foo32, aux32)
    deallocate (T1, T1_dag)


  end subroutine selfenergy


!  *******************************************************************  !
!                               freeleads                               !
!  *******************************************************************  !
!  Description: free allocated matrices and vectors.                    !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer nspin                : Number of spin components             !
!  *******************************************************************  !
  subroutine freeleads

!
!   Modules
!
    use idsrdr_options,  only: nspin


!   Free memory.
    deallocate (S0_L, S1_L)
    deallocate (H0_L, H1_L)
    deallocate (S0_R, S1_R)
    deallocate (H0_R, H1_R)
    deallocate (numberL, numberR)
    deallocate (side_rankL, side_rankR)
    deallocate (QL, VHL, VSL)
    deallocate (QR, VHR, VSR)
    deallocate (Sigma_L, Sigma_R)


  end subroutine freeleads


!  *******************************************************************  !


END MODULE idsrdr_leads
