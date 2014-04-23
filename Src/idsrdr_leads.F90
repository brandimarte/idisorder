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
  use idsrdr_io,       only: 
  use idsrdr_string,   only: 
  use idsrdr_engrid,   only: 
  use fdf

  implicit none
  
  PUBLIC  :: NL, NR, numberL, numberR, EfLead, S0_L, S0_R, S1_L, S1_R,  &
             H0_L, H0_R, H1_L, H1_R, VHL, VHR, VSL, VSR, QL, QR,        &
             Sigma_L, Sigma_R, side_rankL, side_rankR,                  &
             readleads, leadsSelfEn, freeleads
  PRIVATE :: selfenergy, selfenergy2, kofewrap, kofec, kofe_svdlr2,     &
             get_options_kofe, svdm, addnoise, geigenvalues,            &
             geigenvalues2, alphaofz

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
!  integer MPI_Comm_MyWorld     : MPI communicator                      !
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
          call STRpaste (directory, 'bulklft', file)
          call zhsunits (nspin, nspin, NL, nsc, 1, 0, 0, .true.,        &
                         1, temp, H0_L, H1_L, S0_L, S1_L, file)
          call STRpaste (directory, 'bulkrgt', file)
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
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode                 : True if it is the I/O node          !
!  integer SigmaMethod            : Method for calculating self-energy  !
!  ****************************** INPUT ******************************  !
!  integer ispin                  : Spin component index                !
!  real*8 Ei                      : Energy grid point                   !
!  ***************************** OUTPUT ******************************  !
!  complex*8 Sigma_L(NL,NL)       : Left-lead self-energy               !
!  complex*8 Sigma_R(NR,NR)       : Right-lead self-energy              !
!  ************************ OUTPUT TO MODULES ************************  !
!  real*8 :: deltaEn              : Energy step size                    !
!  *******************************************************************  !
  subroutine leadsSelfEn (Ei, ispin)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_options,  only: SigmaMethod, deltaEn

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    integer :: INFO, NCHAN
    real(8) :: deltaene

    if (IOnode) then
       write (6,'(/,a,f8.5,a,i1,a,/)') '[energy(Ry),spin] = [',         &
            Ei, ',', ispin, ']'
       write (6,'(a)', advance='no')                                    &
            '      computing leads self-energies... '
    endif

    if (SigmaMethod == 0) then
       call selfenergy ('L', side_rankL(ispin), numberL(ispin), NL,     &
                        DCMPLX(Ei), H0_L(:,:,ispin), VHL(:,:,ispin),    &
                        S0_L, VSL(:,:,ispin), Sigma_L, QL(:,:,ispin),   &
                        INFO, NCHAN)

       if (INFO == 0) then
          call selfenergy ('R', side_rankR(ispin), numberR(ispin), NR,  &
                           DCMPLX(Ei), H0_R(:,:,ispin), VHR(:,:,ispin), &
                           S0_R, VSR(:,:,ispin), Sigma_R,               &
                           QR(:,:,ispin), INFO, NCHAN)
       endif

    else ! SigmaMethod == 1

       INFO = 0
       deltaene = deltaEn
       call selfenergy2 ('L', NL, DCMPLX(Ei), H0_L(:,:,ispin),          &
                         H1_L(:,:,ispin), S0_L, S1_L,                   &
                         Sigma_L, NCHAN, deltaene)

       deltaene = deltaEn
       call selfenergy2 ('R', NR, DCMPLX(Ei), H0_R(:,:,ispin),          &
                         H1_R(:,:,ispin), S0_R, S1_R,                   &
                         Sigma_R, NCHAN, deltaene)          

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
!  complex*8 H0(N2,N2)        : Hamiltonian of the lead lead            !
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


!   Input variables.
    integer, intent(in) :: number, N2
    integer, intent(out) :: INFO, nrchan
    complex(8), intent(in) :: Ei
    complex(8), dimension (N2,N2), intent(in) :: H0, H1, S0, S1, Q
    complex(8), dimension (N2,N2), intent(out) :: Sigma
    character(len=1), intent(in) :: SIDE, side_rank

!   Local variables.
    integer :: N3, nlchan, I, J
    complex(8), parameter :: smear = (0.d0,0.d0)
    complex(8), allocatable, dimension (:,:) :: T1_aux, T1_dag, T1,     &
                                                H0_aux, H1_aux,         &
                                                H1_dag_aux, h0corr,     &
                                                Gr_2, Gr_square
    complex(8), allocatable, dimension (:,:) :: Gr
    complex(8), allocatable, dimension (:,:) :: Gr_1
    complex(8), allocatable, dimension (:,:) :: foo23
    complex(8), allocatable, dimension (:,:) :: foo32, aux32
    external :: DECIMATE_LEADS, LEADS, HI_zgemm, HI_zgeInvert

!   Allocate matrices and arrays.
    allocate (T1_aux(N2,N2), H0_aux(N2,N2), H1_aux(N2,N2),              &
              H1_dag_aux(N2,N2), h0corr(N2,N2))
    allocate (Gr(2*(N2-number),2*(N2-number)))
    allocate (Gr_1(N2-number,N2-number))

!   Initialize variables.
    N3 = N2 - number
    H0_aux = H0 - Ei*S0
    Gr = (0.d0,0.d0)
    nrchan = 0
    nlchan = 0

!   ('T1_aux = Q^dagger * H0_aux')
    call HI_zgemm ('C', 'N', N2, N2, N2, (1.d0,0.d0), Q, N2,            &
                   H0_aux, N2, (0.d0,0.d0), T1_aux, N2)
!   ('H0_aux = T1_aux * Q')
    call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), T1_aux, N2,       &
                   Q, N2, (0.d0,0.d0), H0_aux, N2)

    IF (side_rank /= '0') THEN

!      Compute 'H1_dag_aux' according to the direction of decimation.
       If (side_rank == 'N') Then

          T1_aux = H1 - Ei*S1

!         ('H1_aux = Q^dagger * T1_aux')
          call HI_zgemm ('C', 'N', N2, N2, N2, (1.d0,0.d0), Q, N2,      &
                         T1_aux, N2, (0.d0,0.d0), H1_aux, N2)

!         ('T1_aux = H1^dagger - Ei*S1^dagger')
          do I = 1,N2
             do J = 1,N2
                T1_aux(I,J) = DCONJG(H1(J,I)) - Ei*DCONJG(S1(J,I))
             enddo
          enddo

!         ('H1_dag_aux = T1_aux * Q')
          call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), T1_aux, N2, &
                         Q, N2, (0.d0,0.d0), H1_dag_aux, N2)

       Else

          T1_aux = H1 - Ei*S1

!         ('H1_aux = T1_aux * Q')
          call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), T1_aux, N2, &
                         Q, N2, (0.d0,0.d0), H1_aux, N2)

!         ('T1_aux = H1^dagger - Ei*S1^dagger')
          do I = 1,N2
             do J = 1,N2
                T1_aux(I,J) = DCONJG(H1(J,I)) - Ei*DCONJG(S1(J,I))
             enddo
          enddo

!         ('H1_dag_aux = Q^dagger * T1_aux')
          call HI_zgemm ('C', 'N', N2, N2, N2, (1.d0,0.d0), Q, N2,      &
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

    IF (INFO == 0) THEN
       If (side == 'L') Then
          if (side_rank == 'N') then

!            Allocate matrix.
             allocate (Gr_2(N2,N2))

             T1(1:N3,:) = H1_aux(number+1:N2,:)
             T1_dag(:,1:N3) = H1_dag_aux(:,number+1:N2)
             Gr_1 = Gr(1:N3,1:N3)

!            ('T1_aux = Gr_1 * T1')
             foo32 = T1(1:N3,:)
             call HI_zgemm ('N', 'N', N3, N2, N3,( 1.d0,0.d0),          &
                            Gr_1, N3, foo32, N3, (0.d0,0.d0), aux32, N3)

!            ('Sigma = T1_dag * T1_aux')
             foo23 = T1_dag(:,1:N3)
             call HI_zgemm ('N', 'N', N2, N2, N3, (1.d0,0.d0),          &
                            foo23, N2, aux32, N3, (0.d0,0.d0), Sigma, N2)

!            (Gr_2 = (-h0corr - Sigma)^-1')
             Gr_2 = -h0corr - Sigma
             call HI_zgeInvert (Gr_2, N2)

             if (DIMAG(smear) > 1.d-7) then ! PB: it will never happen
                                            ! since 'smear'is set to 0?

!               Allocate matrix.
                allocate (Gr_square(N2,N2))

!               ('Gr_square = Gr_2 * Gr_2')
                call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), Gr_2, &
                               N2, Gr_2, N2, (0.d0,0.d0), Gr_square, N2)

                Gr_2 = Gr_2 - smear*Gr_square

!               Free memory.
                deallocate (Gr_square)

             endif

             T1 = H1 - Ei*S1

!            ('T1_aux = T1 * Q^dagger')
             call HI_zgemm ('N', 'C', N2, N2, N2, (1.d0,0.d0), T1, N2,  &
                            Q, N2, (0.d0,0.d0), T1_aux, N2)
!            ('T1 = Q^dagger * T1_aux')
             call HI_zgemm ('C', 'N', N2, N2, N2, (1.d0,0.d0), Q, N2,   &
                            T1_aux, N2, (0.d0,0.d0), T1, N2)     

!            ('T1_dag = H1^dagger - Ei*S1^dagger')
             do i = 1,N2
                do j = 1,N2
                   T1_dag(I,J) = DCONJG(H1(J,I)) - Ei*DCONJG(S1(J,I))
                enddo
             enddo

!            ('T1_aux = T1_dag * Q')
             call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0),          &
                            T1_dag, N2, Q, N2, (0.d0,0.d0), T1_aux, N2)
!            ('T1_dag = Q * T1_aux')
             call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), Q, N2,   &
                            T1_aux, N2, (0.d0,0.d0), T1_dag, N2)     

!            ('T1_aux = Gr_2 * T1')
             call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0),          &
                            Gr_2, N2, T1, N2, (0.d0,0.d0), T1_aux, N2)

!            ('Sigma = T1_dag * T1_aux')
             call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), T1_dag,  &
                            N2, T1_aux, N2, (0.d0,0.d0), Sigma, N2)

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
                call HI_zgeInvert (Gr_1, N3)

!               (Gr_1 = (Gr_1 - h0corr)^-1')
                Gr_1 = Gr_1 - h0corr(number+1:N2,number+1:N2)
                call HI_zgeInvert (Gr_1, N3)

             endif

             if (DIMAG(smear) > 1.d-7) then ! PB: it will never happen
                                            ! since 'smear'is set to 0?

!               Allocate matrix.
                allocate (Gr_square(N3,N3))

!               ('Gr_square = Gr_1 * Gr_1')
                call HI_zgemm ('N', 'N', N3, N3, N3, (1.d0,0.d0),       &
                               Gr_1, N3, Gr_1, N3, (0.d0,0.d0),         &
                               Gr_square, N3)

                Gr_1 = Gr_1 - smear*Gr_square

!               Free memory.
                deallocate (Gr_square)

             endif

!            ('T1_aux = Gr_1 * T1')
             foo32 = T1(1:N3,:)
             call HI_zgemm ('N', 'N', N3, N2, N3, (1.d0,0.d0),          &
                            Gr_1, N3, foo32, N3, (0.d0,0.d0), aux32, N3)

!            ('Sigma = T1_dag * T1_aux')
             foo23 = T1_dag(:,1:N3)
             call HI_zgemm ('N', 'N', N2, N2, N3, (1.d0,0.d0),          &
                            foo23, N2, aux32, N3, (0.d0,0.d0), Sigma, N2)
 
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
             call HI_zgemm ('N', 'N', N3, N2, N3,( 1.d0,0.d0),          &
                            Gr_1, N3, foo32, N3, (0.d0,0.d0), aux32, N3)

!            ('Sigma = T1 * T1_aux')
             foo23 = T1(:,1:N3)
             call HI_zgemm ('N', 'N', N2, N2, N3, (1.d0,0.d0),          &
                            foo23, N2, aux32, N3, (0.d0,0.d0), Sigma, N2)

!            (Gr_2 = (-h0corr - Sigma)^-1')
             Gr_2 = -h0corr - Sigma
             call HI_zgeInvert (Gr_2, N2)

             if (DIMAG(smear) > 1.d-7) then ! PB: it will never happen
                                            ! since 'smear'is set to 0?

!               Allocate matrix.
                allocate (Gr_square(N2,N2))

!               ('Gr_square = Gr_2 * Gr_2')
                call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), Gr_2, &
                               N2, Gr_2, N2, (0.d0,0.d0), Gr_square, N2)

                Gr_2 = Gr_2 - smear*Gr_square

!               Free memory.
                deallocate (Gr_square)

             endif
                                               
             T1 = H1 - Ei*S1

!            ('T1_aux = T1 * Q')
             call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), T1, N2,  &
                            Q, N2, (0.d0,0.d0), T1_aux, N2)
!            ('T1 = Q * T1_aux')
             call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), Q, N2,   &
                            T1_aux, N2, (0.d0,0.d0), T1, N2)     

!            ('T1_dag = H1^dagger - Ei*S1^dagger')
             do I = 1,N2
                do J = 1,N2
                   T1_dag(I,J) =  H1(J,I) - Ei*S1(J,I)
                enddo
             enddo

!            ('T1_aux = T1_dag * Q^dagger')
             call HI_zgemm ('N', 'C', N2, N2, N2, (1.d0,0.d0),          &
                            T1_dag, N2, Q, N2, (0.d0,0.d0), T1_aux, N2)
!            ('T1_dag = Q^dagger * T1_aux')
             call HI_zgemm ('C', 'N', N2, N2, N2, (1.d0,0.d0), Q, N2,   &
                            T1_aux, N2, (0.d0,0.d0), T1_dag, N2)     

!            ('T1_aux = Gr_2 * T1_dag')
             call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), Gr_2,    &
                            N2, T1_dag, N2, (0.d0,0.d0), T1_aux, N2)

!            ('Sigma = T1 * T1_aux')
             call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0), T1, N2,  &
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
                call HI_zgeInvert (Gr_1, N3)

!               (Gr_1 = (Gr_1 - h0corr)^-1')
                Gr_1 = Gr_1 - h0corr(number+1:N2,number+1:N2)
                call HI_zgeInvert (Gr_1, N3)

             endif

             iF (DIMAG(smear) > 1.d-7) then ! PB: it will never happen
                                            ! since 'smear'is set to 0?

!               Allocate matrix.
                allocate (Gr_square(N3,N3))

!               ('Gr_square = Gr_1 * Gr_1')
                call HI_zgemm ('N', 'N', N3, N3, N3, (1.d0,0.d0),       &
                               Gr_1, N3, Gr_1, N3, (0.d0,0.d0),         &
                               Gr_square, N3)

                Gr_1 = Gr_1 - smear*Gr_square

!               Free memory.
                deallocate (Gr_square)

             endif

!            ('T1_aux = Gr_1 * T1_dag')
             foo32 = T1_dag(1:N3,:)
             call HI_zgemm ('N', 'N', N3, N2, N3, (1.d0,0.d0),          &
                            Gr_1, N3, foo32, N3, (0.d0,0.d0), aux32, N3)

!            ('Sigma = T1 * T1_aux')
             foo23 = T1(:,1:N3)
             call HI_zgemm ('N', 'N', N2, N2, N3, (1.d0,0.d0),          &
                            foo23, N2, aux32, N3, (0.d0,0.d0), Sigma, N2)

          endif ! if (side_rank == 'C')

       EndIf ! If (side == 'L')

    ELSE

       RETURN

    ENDIF ! IF (INFO == 0)

!   Free memory.
    deallocate (T1_aux, H1_aux, H1_dag_aux, h0corr)
    deallocate (Gr)
    deallocate (Gr_1)
    deallocate (T1, T1_dag)
    deallocate (foo32, aux32)
    deallocate (foo23)


  end subroutine selfenergy


!  *******************************************************************  !
!                              selfenergy2                              !
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
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode             : True if it is the I/O node              !
!  integer MPI_Comm_MyWorld   : MPI communicator                        !
!  real*8 imDelta             : Small imaginary part                    !
!  ****************************** INPUT ******************************  !
!  character(len=1) SIDE      : 'L' (left) or 'R' (right) depending on  !
!                               the self-energy calculated              !
!  integer N2                 : Number of basis orbitals on the lead    !
!  complex*8 Ei               : Energy                                  !
!  complex*8 H0(N2,N2)        : Hamiltonian of the lead lead            !
!  complex*8 H1(N2,N2)        : Coupling Matrix                         !
!  complex*8 S0(N2,N2)        : On-site Overlap Matrix                  !
!  complex*8 S1(N2,N2)        : Neighbour Overlap Matrix                !
!  ************************** INPUT/OUTPUT ***************************  !
!  real*8 deltaene            : Energy "step" size                      !
!  ***************************** OUTPUT ******************************  !
!  complex*8 Sigma            : Calculated self-energy                  !
!  integer nrchan             : Number of open channels                 !
!  *******************************************************************  !
  subroutine selfenergy2 (SIDE, N2, Ei, H0, H1, S0, S1,                 &
                          Sigma, nrchan, deltaene)

!
!   Modules
!
#ifdef MPI
    use parallel,        only: IOnode, MPI_Comm_MyWorld
#else
    use parallel,        only: IOnode
#endif
    use idsrdr_options,  only: imDelta
    use fdf

#ifdef MPI
    include "mpif.h"
#endif

!   Input variables.
    integer, intent(in) :: N2
    integer, intent(out) :: nrchan
    real(8), intent(inout) :: deltaene
    complex(8), intent(in) :: Ei
    complex(8), dimension (N2,N2), intent(in) :: H0, H1, S0, S1
    complex(8), dimension (N2,N2), intent(out) :: Sigma
    character(len=1), intent(in) :: SIDE

!   Local variables.
    integer :: j, neig
    real(8) :: dsigma, maxsigma, dsigmar, dgammam
    real(8), parameter :: pi = 3.1415926535897932384626433832795028841D0
    complex(8) :: Ei0, tracegs0, tracegsm1, alpha
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8), allocatable, dimension (:) :: zeig
    complex(8), allocatable, dimension (:,:) :: k0, k1, km1, k1o, k1ob, &
                                                Gr_2, Gr_2t, Aux,       &
                                                rho0, rhom1
    logical, save :: pdosgs, skipright
    logical, save :: firstcall = .true.
    external :: HI_zgeInvert, HI_zgemm
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif


!   Initialize variables.
    Ei0 = Ei
    if (firstcall) then

       firstcall = .false.

       if (IOnode) then
          pdosgs = fdf_boolean ('Sigma.PDOS', .false.)
          skipright = fdf_boolean ('Sigma.SkipRight', .false.)
       endif

#ifdef MPI
       call MPI_Bcast (pdosgs, 1, MPI_Logical, 0,                       &
                       MPI_Comm_MyWorld, MPIerror)
       call MPI_Bcast (skipright, 1, MPI_logical, 0,                    &
                       MPI_Comm_MyWorld, MPIerror)
#endif
    endif

!   Allocate matrices and arrays.
    allocate (k0(N2,N2), k1(N2,N2), km1(N2,N2), k1o(N2,N2),             &
              k1ob(N2,N2), Gr_2(N2,N2), Gr_2t(N2,N2),                   &
              rho0(N2,N2), rhom1(N2,N2))
    allocate (zeig(2*N2))

!   Check energy step size.
    if (deltaene > 10d-4) deltaene = 1.d-4 ! shouldn't it be:
                                           ! if (deltaene > 1.d-4) ?

!   Set energy and hamiltonians.
    Ei0 = Ei0 + zi * imDelta
    k0 = H0 - Ei0*S0
    k1 = H1 - Ei0*S1
    km1 = DCONJG(TRANSPOSE(H1)) - Ei0 * DCONJG(TRANSPOSE(S1))

    if (skipright .and. SIDE == 'R' ) then

!      Free memory.
       deallocate (k0, k1, km1)
       deallocate (zeig)

       RETURN

    endif

!   Compute lead self-energy.
    call kofewrap (SIDE, 'S', k0, k1, km1, N2,                          &
                   zeig, neig, Ei0, dsigma, nrchan)

    if (pdosgs .and. SIDE == 'L') then

!      Allocate matrices and arrays.
       allocate (k1o(N2,N2), k1ob(N2,N2), Gr_2(N2,N2), Gr_2t(N2,N2),    &
                 rho0(N2,N2), rhom1(N2,N2), Aux(N2,N2))

       k1o = H0 - Ei0 * S0
       k1 = H1 - Ei0 * S1
       km1 = DCONJG(TRANSPOSE(H1)) - Ei0 * DCONJG(TRANSPOSE(S1))

       dgammam = MAXVAL(ABS(DREAL(k0-DCONJG(TRANSPOSE(k0)))))
       dgammam = 0D0
       do j = 1,n2
          dgammam = dgammam + DREAL(k0(j,j))
       enddo

       k1o = - k1o - k0
       call HI_zgeInvert (k1o, n2)
       Gr_2 = k1o

!      ('k1ob = km1 * Gr_2')
       call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0),                &
                      km1, N2, Gr_2, N2, (0.d0,0.d0), k1ob, N2)

!      ('k1o = k1ob * k1')
       call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0),                &
                      k1ob, N2, k1, N2, (0.d0,0.d0), k1o, N2)

       k1ob = k0 - k1o
       dsigma = MAXVAL(ABS(k1ob))
       maxsigma = MAXVAL(ABS(k0))
       dsigmar = dsigma / maxsigma

       Gr_2t = DCONJG(TRANSPOSE(Gr_2))
       rho0 = (0.5d0 * zi / pi) * (Gr_2 - Gr_2t)

!    **('rhom1 = (0.5*zi/pi) * (Gr_2t*k1^dagger*Gr_2 - Gr_2*km1*Gr_2)')**

!      ('Aux = - (0.5*zi/pi) * Gr_2 * km1')
       alpha = - 0.5d0 * zi / pi
       call HI_zgemm ('N', 'N', N2, N2, N2, alpha,                      &
                      Gr_2, N2, km1, N2, (0.d0,0.d0), Aux, N2)

!      ('rhom1 = Aux * Gr_2')
       call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0),                &
                      Aux, N2, Gr_2, N2, (0.d0,0.d0), rhom1, N2)

!      ('Aux = (0.5*zi/pi) * Gr_2t * k1^dagger')
       alpha = 0.5d0 * zi / pi
       call HI_zgemm ('N', 'C', N2, N2, N2, alpha,                      &
                      Gr_2t, N2, k1, N2, (0.d0,0.d0), Aux, N2)

!      ('rhom1 = Aux * Gr_2 + rhom1')
       call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0),                &
                      Aux, N2, Gr_2, N2, (1.d0,0.d0), rhom1, N2)

!      ('rho0 = rho0 * S0')
       Aux = rho0
       call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0),                &
                      Aux, N2, S0, N2, (0.d0,0.d0), rho0, N2)

!      ('rhom1 = rhom1 * S1')
       Aux = rhom1
       call HI_zgemm ('N', 'N', N2, N2, N2, (1.d0,0.d0),                &
                      Aux, N2, S1, N2, (0.d0,0.d0), rhom1, N2)

       tracegs0 = 0.d0
       do j = 1,n2
          tracegs0 = tracegs0 + rho0(j,j)
       enddo
       tracegsm1 = 0.d0
       do j=1,n2
          tracegsm1 = tracegsm1 + rhom1(j,j)
       enddo
       write (12347,'(a," ",e," ",9(e," "))') "trrss = ", DREAL(Ei0),   &
            DREAL(tracegs0), DREAL(tracegs0+tracegsm1),                 &
            DREAL(tracegsm1), DIMAG(tracegsm1), dsigma, dsigmar,        &
            maxsigma, DIMAG(Ei0), dgammam

!      Free memory.
       deallocate (k1o, k1ob, Gr_2, Gr_2t, rho0, rhom1, Aux)

    endif ! if (pdosgs .and. SIDE == 'L')

    Sigma = k0

!   Free memory.
    deallocate (k0, k1, km1)
    deallocate (zeig)


  end subroutine selfenergy2


!  *******************************************************************  !
!                               kofewrap                                !
!  *******************************************************************  !
!  Description:                                                         !
!                                                                       !
!  From SMEAGOL code - 2003.                                            !
!                                                                       !
!  Written by Alexandre Reily Rocha, Jun 2003.                          !
!  Computational Spintronics Group                                      !
!  Trinity College Dublin                                               !
!  e-mail: reilya@ift.unesp.br                                          !
!  ***************************** HISTORY *****************************  !
!  Original version:    June 2003                                       !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode                : True if it is the I/O node           !
!  integer Node                  : Actual node (MPI_Comm_rank)          !
!  integer MPI_Comm_MyWorld      : MPI communicator                     !
!  ****************************** INPUT ******************************  !
!  character(len=1) SIDE         : 'L' (left) or 'R' (right) depending  !
!                                  on the self-energy calculated        !
!  character(len=1) gf           :                                      !
!  complex*8 k0(n,n)             : 'H0 - zi*Ei*S0'                      !
!  complex*8 k1in(n,n)           : 'H1 - zi*Ei*S1'                      !
!  complex*8 km1in(n,n)          : 'H1^dagger - zi*Ei*S1^dagger'        !
!  integer n                     : # of basis orbitals on the lead      !
!  complex*8 zeig(2*n)           :                                      !
!  integer neig                  :                                      !
!  complex*8 Ei                  : Energy                               !
!  real*8 dsigma                 :                                      !
!  ************************** INPUT/OUTPUT ***************************  !
!  complex*8 k0(n,n)          : 'H0 - zi*Ei*S0' (return self-energy)    !
!  ***************************** OUTPUT ******************************  !
!  integer nrchan             : Number of open channels                 !
!  *******************************************************************  !
  subroutine kofewrap (SIDE, gf, k0, k1in, km1in, n,                    &
                       zeig, neig, Ei, dsigma, nrchan)

!
!   Modules
!
#ifdef MPI
    use parallel,        only: IOnode, Node, MPI_Comm_MyWorld
#else
    use parallel,        only: IOnode, Node
#endif
    use fdf

#ifdef MPI
    include "mpif.h"
#endif

!   Input variables.
    integer, intent(in) :: n
    integer, intent(inout) :: neig
    integer, intent(out) :: nrchan
    real(8), intent(inout) :: dsigma
    complex(8), intent(in) :: Ei
    complex(8), dimension (n,n), intent(in) :: k1in, km1in
    complex(8), dimension (n,n), intent(inout) :: k0
    complex(8), dimension (2*n), intent(out) :: zeig
    character(len=1), intent(in) :: SIDE, gf

!   Local variables.
    integer :: i
    real(8) :: svdtolzi2
    real(8), save :: svdtolzi, dsigmamax, eimag
    complex(8), allocatable, dimension (:,:) :: k0in
    logical, save :: firstcall = .true.
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

!   Initialize variables.
    if (firstcall) then

       firstcall = .false.

       if (IOnode) then
          svdtolzi = fdf_double ('Sigma.SVDTolZero', 0.d0)
          dsigmamax = fdf_double ('Sigma.DMax', 1.d-5)
          eimag = fdf_double ('Sigma.EImag', 0.d-5) ! 1.d-5?
       endif
#ifdef MPI
       call MPI_Bcast (svdtolzi, 1, MPI_Double_Precision, 0,            &
                       MPI_Comm_MyWorld, MPIerror)
       call MPI_Bcast (dsigmamax, 1, MPI_Double_Precision, 0,           &
                       MPI_Comm_MyWorld, MPIerror)
       call MPI_Bcast (eimag, 1, MPI_Double_Precision, 0,               &
                       MPI_Comm_MyWorld, MPIerror)
#endif
    endif
    svdtolzi2 = svdtolzi

!   Allocate auxiliary matrix.
    allocate (k0in(n,n))

    if (svdtolzi2 == 0.d0) then 

       call kofec (SIDE, gf, k0, k1in, km1in, n,                        &
                   zeig, neig, Ei, dsigma, nrchan)
    else

       k0in = k0
    
       call kofe_svdlr2 (SIDE, k0, k1in, km1in, n, neig,                &
                         Ei, svdtolzi2, dsigma, nrchan)

       if (dsigma > dsigmamax .and. dsigma > 1.d1*eimag) then

          write (12347,*) "svdtr1=", svdtolzi2, dsigma, 1.d1*eimag

          do i = 1,3

             svdtolzi2 = svdtolzi2 / 1.d2
             if (svdtolzi2 < 1.d-13) exit
             k0 = k0in

             call kofe_svdlr2 (SIDE, k0, k1in, km1in, n, neig,          &
                               Ei, svdtolzi2, dsigma, nrchan)

             if (dsigma < dsigmamax .or. dsigma < 1.d1*eimag) exit

             write (12347,*) "svdtr1=", i, svdtolzi2, dsigma, 1.d1*eimag

          enddo

          write(12347,*) "svdtexit=", svdtolzi2, dsigma, Node

       endif

    endif ! if (svdtolzi2 == 0.d0)

!   Free memory.
    deallocate (k0in)


  end subroutine kofewrap


!  *******************************************************************  !
!                                 kofec                                 !
!  *******************************************************************  !
!  Description:                                                         !
!                                                                       !
!  From SMEAGOL code - 2003.                                            !
!                                                                       !
!  Written by Alexandre Reily Rocha, Jun 2003.                          !
!  Computational Spintronics Group                                      !
!  Trinity College Dublin                                               !
!  e-mail: reilya@ift.unesp.br                                          !
!  ***************************** HISTORY *****************************  !
!  Original version:    June 2003                                       !
!  *********************** INPUT FROM MODULES ************************  !
!  integer Node                  : Actual node (MPI_Comm_rank)          !
!  ****************************** INPUT ******************************  !
!  character(len=1) SIDE         : 'L' (left) or 'R' (right) depending  !
!                                  on the self-energy calculated        !
!  character(len=1) gf           :                                      !
!  complex*8 k0(n,n)             : 'H0 - zi*Ei*S0'                      !
!  complex*8 k1in(n,n)           : 'H1 - zi*Ei*S1'                      !
!  complex*8 km1(n,n)            : 'H1^dagger - zi*Ei*S1^dagger'        !
!  integer n                     : # of basis orbitals on the lead      !
!  complex*8 zeig(2*n)           :                                      !
!  integer neig                  :                                      !
!  complex*8 Ei                  : Energy                               !
!  real*8 dsigma                 :                                      !
!  ************************** INPUT/OUTPUT ***************************  !
!  complex*8 k0(n,n)          : 'H0 - zi*Ei*S0' (return self-energy)    !
!  ***************************** OUTPUT ******************************  !
!  integer nrchan             : Number of open channels                 !
!  *******************************************************************  !
  subroutine kofec (SIDE, gf, k0, k1in, km1in, n,                       &
                    zeig, neig, Ei, dsigmam, nrchanm)

!
!   Modules
!
    use parallel,        only: Node

!   Input variables.
    integer, intent(in) :: n
    integer, intent(inout) :: neig
    integer, intent(out) :: nrchanm
    real(8), intent(inout) :: dsigmam
    complex(8), intent(in) :: Ei
    complex(8), dimension (n,n), intent(in) :: k1in, km1in
    complex(8), dimension (n,n), intent(inout) :: k0
    complex(8), dimension (2*n), intent(out) :: zeig
    character(len=1), intent(in) :: SIDE, gf

!   Local variables.
    integer :: nrchan, i, i1, i2, il, ir, iloc, iroc, nls, nrs, INFO
               
    integer, dimension (n) :: vml, vmr, IPIV
    real(8) :: normvr, normvr1, normvr2, svdtol, lz, alphak,            &
               dsigma, signs, ANORM, RCOND, svdtolm
    real(8), save :: svdtolmax, svdtolmin, dsigmamax,                   &
                     tolki, rnoise, skipsvd
    real(8), allocatable, dimension (:) :: rwork, rwork3
    complex(8) :: zvu, vg
    complex(8), parameter :: alpha = (1.d0,0.d0)
    complex(8), dimension (n) :: zvl, zvr, vri, vrbuf, vritd, vli, vlitd
    complex(8), dimension (2*n) :: rworkc, walpha, wbeta
    complex(8), dimension (n*n) :: worki
    complex(8), dimension (n,n) :: k1, grfunct, k1input, k0input,       &
                                   matout, km1, km1input, al, ar,       &
                                   sigma, vlg, vrg, vlgb, vrgb
    complex(8), dimension (n,2*n) :: vr
    complex(8), dimension (2*n,2*n) :: vr2
    complex(8), allocatable, dimension (:) :: zv, work, work3
    complex(8), allocatable, dimension (:,:) :: k1o, k1ob, k1t, m2s,    &
                                                m2sb, m2sbuf, m2sbbuf,  &
                                                Aux
    logical :: callsvd, usehinv1, calcsucceed
    logical, save :: usehinv, usevinv ,callsvdin
    logical, save :: firstcall = .true.
    real(8), external :: ZLANGE
    external :: ZGETRF, ZGETRI, ZGECON, ZGEMV, HI_zgeInvert, HI_zgemm

!   Initialize variables.
    nrchanm = 0
    if (SIDE == 'R') then
       k1 = km1in
       km1 = k1in
       signs = -1.d0
    else ! left lead
       k1 = k1in
       km1 = km1in
       signs = 1.d0
    endif
    k1input = k1
    km1input = km1
    k0input = k0
    svdtol = 0.d0
    svdtolm = svdtol
    dsigmam = 0.d0
    usehinv1 = .false.
    calcsucceed = .false.
    if (firstcall) then
       firstcall = .false.
       rnoise = 0.d0
       call get_options_kofe (usehinv, usevinv, tolki, svdtolmax,       &
                              svdtolmin, dsigmamax, callsvdin,          &
                              rnoise, skipsvd)
    endif
    callsvd = callsvdin

!   Allocate matrices and arrays.
    allocate (rwork(5*n))
    allocate (rwork3(16*n))
    allocate (zv(2*n))
    allocate (work(8*n))
    allocate (work3(8*n))
    allocate (k1o(n,n), k1ob(n,n), k1t(n,n), Aux(n,n))
    allocate (m2s(2*n,2*n), m2sb(2*n,2*n),                              &
              m2sbuf(2*n,2*n), m2sbbuf(2*n,2*n))


102 if (callsvd) then

       if (svdtol == 0.d0) svdtol = svdtolmin / skipsvd
       svdtol = svdtol * skipsvd
       k1 = k1input
       km1 = km1input
       k0 = k0input

       if (.not. usehinv) then
          call svdm (k1, k1o, n, svdtol, rnoise*svdtol, .false., .true.)
!!$          call addnoise (k1, n, rnoise*svdtol)
       else
          call svdm (k1, k1o, n, svdtol, rnoise*svdtol, .true., .false.)
!!$          call addnoise (k0, n, rnoise*svdtol)
!!$          call addnoise (k1, n, rnoise*svdtol)
!!$          km1 = km1input
!!$          call addnoise (km1, n, rnoise*svdtol)
!!$          k1o = k1
!!$          call ZGETRF (n, n, k1o, n, IPIV, INFO)
!!$          call ZGETRI (n, k1o, n, IPIV, worki, n*n, INFO)
       endif
!!$       k0 = 0.5d0 * (k0 + DCONJG(TRANSPOSE(k0)))

       if (rnoise /= 0.d0) then
          call addnoise (k0, n, rnoise*svdtol)
          call addnoise (km1, n, rnoise*svdtol)
       endif

    else ! callsvd == .false.

       callsvd = .true.
       svdtol = 0.d0
       svdtolm = svdtol
       k1 = k1input
       km1 = km1input

       if (usehinv) then
          k0 = k0input
!!$          k0 = 0.5d0 * (k0 + DCONJG(TRANSPOSE(k0)))
          k1o = k1
          call HI_zgeInvert (k1o, n)

       endif

       if (rnoise /= 0.d0) then
          call addnoise (k0, n, rnoise*1.d-17)
          call addnoise (k1, n, rnoise*1.d-17)
          call addnoise (km1, n, rnoise*1.d-17)
       endif

    endif ! if (callsvd) then

!!$    km1 = DCONJG(transpose(k1))

110 if (usehinv .or. usehinv1) then

       usehinv1 = .false.

       m2s = 0.d0

!      ('Aux = - k1o * k0')
       call HI_zgemm ('N', 'N', n, n, n, (-1.d0,0.d0),                  &
                      k1o, n, k0, n, (0.d0,0.d0), Aux, n)
       m2s(1:n,1:n) = Aux

!      ('Aux = - k1o * km1')
       call HI_zgemm ('N', 'N', n, n, n, (-1.d0,0.d0),                  &
                      k1o, n, km1, n, (0.d0,0.d0), Aux, n)
       m2s(1:n,n+1:2*n) = Aux

       do i = 1,n
          m2s(n+i,i) = 1.d0
       enddo
       call geigenvalues2 (m2s, zv, 2*n, vr2, INFO)

!!$       do i1 = 1,2*n
!!$          m2sbbuf(:,i1) = zv(i1) * vr2(:,i1)
!!$       enddo
!!$       m2sbuf = matmul(m2s,vr2)
!!$       m2s = m2sbuf - m2sbbuf
!!$       dsigma = MAXVAL(ABS(m2s))
!!$       write(*,*) "svdresqevp=", dsigma

    else

       m2s = 0.d0
       m2s(1:n,1:n) = k0(1:n,1:n)
       m2s(1:n,1+n:2*n) = km1(1:n,1:n)
       do i = 1,n
          m2s(n+i,i) = 1.d0
       enddo
       m2sb = 0.d0
       m2sb(1:n,1:n) = -k1
       do i = 1,n
          m2sb(n+i,n+i) = 1.d0
       enddo

       call geigenvalues (m2s, m2sb, walpha, wbeta, 2*n, vr2, INFO)

!+++ The following lines should be commented
!+++ out if H1^-1 is to be calculated.
!+++         if (INFO /= 0 .and. svdtol /= 0.d0) then
!+++            call svdm3 (Ei, k1, k1o, n, svdtol)
!+++            km1 = DCONJG(transpose(k1))
!+++            usehinv1 = .true.
!+++            goto 110
!+++         endif
!+++ end H1^-1

       do i1 = 1,2*n
          if (ABS(wbeta(i1)) == 0.d0) then
             zv(i1) = 1.d20
          elseif (ABS(walpha(i1)) == 0.d0) then
             zv(i1) = 1.d-20
          else
             zv(i1) = walpha(i1) / wbeta(i1)
          endif
       enddo

!!$       do i1 = 1,2*n
!!$          m2sbbuf(:,i1) = zv(i1) * vr2(:,i1)
!!$       enddo
!!$       m2sbuf = matmul(m2s,vr2)
!!$       m2sbbuf = matmul(m2sb,m2sbbuf)
!!$       m2s = m2sbuf - m2sbbuf
!!$       dsigma = MAXVAL(ABS(m2s))
!!$       write (12347,*) "resqevp=", dreal(Ei), dsigma, svdtol

    endif ! if (usehinv .or. usehinv1)


    if (INFO /= 0 .and. svdtol < svdtolmax) then
       write (12347,*) "warning: increasing svdtol info",               &
            DREAL(Ei), svdtol, INFO
       goto 102
    elseif (INFO /= 0) then 
       write (12347,*) "warning: last not increasing svdtol info",      &
            DREAL(Ei), svdtol, INFO
       goto 103
    endif

    do i1 = 1,2*n

       normvr1 = SQRT(DOT_PRODUCT(vr2(1:n,i1),vr2(1:n,i1)))
       normvr2 = SQRT(DOT_PRODUCT(vr2(n+1:2*n,i1),vr2(n+1:2*n,i1)))

       if(normvr1 >= normvr2) then
          vr(:,i1) = vr2(1:n,i1)
          normvr = normvr1
       else
          vr(:,i1) = vr2(n+1:2*n,i1)
          normvr = normvr2
       endif
       vr(:,i1) = vr(:,i1) / normvr

    enddo
    neig = 2*n
    zeig(1:neig) = zv(1:neig)

!   Calculate group velocities 'Vg = i ( H_{1} z - H_{-1}/z )'.
    ir = 1
    il = 1
    iroc = 0
    iloc = 0
    nrchan = 0

!!$    do i1 = 1,2*n
!!$       call alphaofz (alphak, lz, zv(i1))
!!$
!!$       if (SIDE == 'R') then
!!$          write (*,*) "nallr=", DREAL(Ei), DIMAG(Ei), il, lz, alphak
!!$       else
!!$          write (*,*) "nalll=", DREAL(Ei), DIMAG(Ei), il, lz, alphak
!!$       endif
!!$    enddo

    do i1 = 1,2*n

       call alphaofz (alphak, lz, zv(i1))

       if (ABS(lz) < tolki) then

!!$          write (*,*) "open channel:", zv(i1)
          k1t = signs * (k1 * zv(i1) - km1 / zv(i1))
          vri = vr(:,i1)

!         ('vrbuf = k1t * vri')
          call ZGEMV ('N', n, n, (1.d0,0.d0), k1t, n, vri, 1,           &
                      (0.d0,0.d0), vrbuf, 1)

          vg = 0.d0
          do i2 = 1,n
             vg = vg + (0.d0,1.d0) *  DCONJG(vri(i2)) * vrbuf(i2)
          enddo
          if (DREAL(vg) > 0.d0) then
!!$             write (*,*) "right going vg = ", DREAL(vg)
             if (ir > n) exit
             vmr(ir) = i1
             ir = ir + 1
             iroc = iroc + 1
             nrchan = nrchan + 1
          else
!!$             write (*,*) "left going vg = ", DIMAG(vg)
             if (il > n) exit
             vml(il) = i1
             il = il + 1
             iloc = iloc + 1
          endif

       else ! ABS(lz) > tolki

!!$          write (*,*) "closed channel: ", zv(i1)
          if (signs*lz < 0.d0) then
!!$             write (*,*) "right decaying state"
             if (ir > n) exit
             vmr(ir) = i1
             ir = ir + 1
          else
!!$             write (*,*) "left decaying state"
             if (il > n) exit
             vml(il) = i1
             il = il + 1
          endif

       endif ! if (ABS(lz) < tolki)

    enddo

    nrs = ir - 1
    nls = il - 1
    if ((nrs /= nls) .or. (nrs+nls /= 2*n)) then
       write (12347,*) "warning:nrs,nls=", DREAL(Ei), nrs,              &
            nls, iroc, iloc, Node, SIDE, svdtol

       if (svdtol < svdtolmax) then
          write (12347,*) "warning:increasing svdtol", DREAL(Ei),       &
               svdtol, Node
          goto 102
       else
          write (12347,*) "warning:problem in last", DREAL(Ei)
          goto 103
       endif

    endif

    do il = 1,n
       vrg(:,il) = vr(:,vmr(il))
       vlg(:,il) = vr(:,vml(il))
       zvr(il) = zv(vmr(il))
       zvl(il) = zv(vml(il))
    enddo

!!$    if (.false.) then ! ignore the following code...
!!$
!!$       do il = 1,n
!!$
!!$          if (SIDE == 'R') then
!!$             zvu = 1.d0 / zvr(il)
!!$          else
!!$             zvu = zvr(il)
!!$          endif
!!$
!!$          call alphaofz (alphak, lz, zvu)
!!$
!!$          if (SIDE == 'R') then
!!$             write (12346,*) "kbsr=", DREAL(Ei), DIMAG(Ei), il,      &
!!$                  lz, alphak, svdtol
!!$          else
!!$             write (12346,*) "kbsrfl=", DREAL(Ei), DIMAG(Ei), il,    &
!!$                  lz, alphak, svdtol
!!$          endif
!!$
!!$          if (SIDE == 'R') then
!!$             zvu = 1.d0 / zvl(il)
!!$          else
!!$             zvu = zvl(il)
!!$          endif
!!$
!!$          call alphaofz (alphak, lz, zvu)
!!$
!!$          if (SIDE == 'L') then
!!$             write (12346,*) "kbsl=", DREAL(Ei), DIMAG(Ei), il,      &
!!$                  lz, alphak, svdtol
!!$          else
!!$             write (12346,*) "kbslfr=", DREAL(Ei), DIMAG(Ei), il,    &
!!$                  lz, alphak, svdtol
!!$          endif
!!$
!!$       enddo
!!$
!!$    endif

!!$-    k1o = vrg
!!$-    call svdmtest (Ei, vrg, vrgb, n, 1.d-8, 0.d0, .true., .false.)
!!$      call addnoise (vrg,n,0d-10)
!!$      !vrg = !!k1o
!!$      dsigma = MAXVAL(ABS(k1o-vrg))
!!$      write (12347,*) "noisechange=", dsigma

    vrgb = vrg
!!$    ANORM = ZLANGE ('1', n, n, vrgb, n, rworkc)
!-    call ZGETRF (n, n, vrgb, n, IPIV, INFO)
!!$    k1o = vrgb
!!$    IPIVs = IPIV
!!$    call ZGECON ('1', n, k1o, n, ANORM, RCOND, walpha, rworkc, INFO)
!-    call ZGETRI (n, vrgb, n, IPIV, worki, n*n, INFO)
    call HI_zgeInvert (vrgb, n)

!!$    if (.false.) then ! ignore the following code...
!!$
!!$       k1ob = 0.d0
!!$       do i = 1,n
!!$          k1ob(i,i) = 1.d0
!!$       enddo
!!$       k1o = matmul(vrg,vrgb)-k1ob
!!$       dsigma = MAXVAL(ABS(k1o))
!!$
!!$       call svdmk (Ei, k1o, k1ob, n, svdtol, rnoise*svdtol,          &
!!$                  .false., .true., kappa, deltaainvn, sn)
!!$       k1o = vrgb
!!$       call svdmk (Ei, k1o, k1ob, n, svdtol, rnoise*svdtol,          &
!!$                   .false., .true., kappa, ainvn, sn)
!!$
!!$       dsigma2 = MAXVAL(ABS(vrgb))
!!$       k1o = vrg
!!$       call svdmk (Ei, k1o, k1ob, n, svdtol, rnoise*svdtol,          &
!!$                   .false., .true., kappa, s1, sn)
!!$       ANORM2 = ZLANGE ('1', n, n, vrgb, n, rworkc)
!!$
!!$       write (12347,*) "RCONDr=", DREAL(Ei), dsigma, dsigma/dsigma2,&
!!$            rcond, 1.d0/kappa, INFO, ANORM, ANORM2,                  &
!!$            s1, sn, deltaainvn, deltaainvn, ainvn
!!$
!!$    endif

!!$!    k1ob  =0.d0
!!$!    do i = 1,n
!!$!       k1ob(i,i) = 1.d0
!!$!    enddo
!!$!    call ZGERFS ('N', n, n, vrg, n, k1o, n, IPIVs, k1ob, n,         &
!!$!                 vrgb, n, FERR, BERR, work4, rwork4, INFO)
!!$    k1o = vrg
!!$    call svdm (k1o, vrgb, n, 1.d-10, svdtol, .true., .false.)

!!$-    call svdmtest (Ei, vlg, vlgb, n, 1.d-8, 0.d0, .true., .false.)
!!$-    call addnoise2 (vlg, n, 1.d-10)
    vlgb = vlg
!-    call ZGETRF (n, n, vlgb, n, IPIV, INFO)
!!$    k1o = vlgb
!!$    IPIVs = IPIV
!-    call ZGETRI (n, vlgb, n, IPIV, worki, n*n, INFO)
    call HI_zgeInvert (vlgb, n)

!!$!    k1ob = 0.d0
!!$!    do i = 1,n
!!$!       k1ob(i,i) = 1.d0
!!$!    enddo
!!$!    call ZGERFS ('N', n, n, vlg, n, k1o, n, IPIVs, k1ob, n,         &
!!$!                 vlgb, n, FERR, BERR, work4, rwork4, INFO)
!!$    k1o = vlg
!!$    call svdm (k1o, vlgb, n, 1.d-10, svdtol, .true., .false.)


    if (SIDE == 'L') then

       al = 0.d0
       do il = 1,n
          vli(:) = vlg(:,il)
          vlitd(:) = vlgb(il,:)
          do i1 = 1,n
             do i2 = 1,n
                al(i1,i2) = al(i1,i2) + vli(i1) * vlitd(i2) / zvl(il)
             enddo
          enddo
       enddo

    else ! side == 'R'

       ar = 0.d0
       do il = 1,n
          zvu = 1.d0 / zvr(il)
          vri(:) = vrg(:,il)
          vritd(:) = vrgb(il,:)
          do i1 = 1,n
             do i2 = 1,n
                ar(i1,i2) = ar(i1,i2) + vri(i1) * vritd(i2) * zvu
             enddo
          enddo
       enddo

    endif


    if (usevinv) then

       if (SIDE == 'L') then
 
          ar = 0.d0
          do il = 1,n
             vri(:) = vrg(:,il)
             vritd(:) = vrgb(il,:)
             do i1 = 1,n
                do i2 = 1,n
                   ar(i1,i2) = ar(i1,i2) + vri(i1) * vritd(i2) * zvr(il)
                enddo
             enddo
          enddo

!         ('k1ob = - al * ar')
          call HI_zgemm ('N', 'N', n, n, n, (-1.d0,0.d0),               &
                         al, n, ar, n, (0.d0,0.d0), k1ob, n)
          do i = 1,n
             k1ob(i,i) = 1.d0 + k1ob(i,i)
          enddo

       else ! side == 'R'

          al = 0.d0
          do il = 1,n
             zvu = 1.d0 / zvl(il)
             vli(:) = vlg(:,il)
             vlitd(:) = vlgb(il,:)
             do i1 = 1,n
                do i2 = 1,n
                   al(i1,i2) = al(i1,i2) + vli(i1) * vlitd(i2) / zvu
                enddo
             enddo
          enddo

!         ('k1ob = - ar * al')
          call HI_zgemm ('N', 'N', n, n, n, (-1.d0,0.d0),               &
                         ar, n, al, n, (0.d0,0.d0), k1ob, n)
          do i = 1,n
             k1ob(i,i) = 1.d0 + k1ob(i,i)
          enddo

       endif

       if(SIDE == 'L')then

          ar = 0.d0
          do il = 1,n
             vri(:) = vrg(:,il)
             vritd(:) = vrgb(il,:)
             do i1 = 1,n
                do i2 = 1,n
                   ar(i1,i2) = ar(i1,i2) + vri(i1) * vritd(i2) / zvr(il)
                enddo
             enddo
          enddo

       else ! side == 'R'

          al = 0.d0
          do il = 1,n
             zvu = 1.d0 / zvl(il)
             vli(:) = vlg(:,il)
             vlitd(:) = vlgb(il,:)
             do i1 = 1,n
                do i2 = 1,n
                   al(i1,i2) = al(i1,i2) + vli(i1) * vlitd(i2) * zvu
                enddo
             enddo
          enddo

       endif

!      ('k1o = km1 * signs * (ar - al)')
       Aux = signs * (ar - al)
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      km1, n, Aux, n, (0.d0,0.d0), k1o, n)

       call HI_zgeInvert (k1o, n)

!!$       tracegs = 0.d0
!!$       do i = 1,n
!!$          tracegs = tracegs + k1o(i,i) - DCONJG(k1o(i,i))
!!$       enddo
!!$       tracegs = -0.5d0 * tracegs / pi
!!$       if (SIDE == 'L') then
!!$          write (12347,*) "trgsl=",                                  &
!!$               DREAL(Ei), DIMAG(tracegs), DREAL(tracegs), n
!!$          write (*,*) "trgsl=",                                      &
!!$               DREAL(Ei), DIMAG(tracegs), DREAL(tracegs), n
!!$       else
!!$          write (12347,*) "trgsr=",                                  &
!!$               DREAL(Ei), DIMAG(tracegs), DREAL(tracegs), n
!!$          write (*,*) "trgsr=",                                      &
!!$          DREAL(Ei), DIMAG(tracegs), DREAL(tracegs), n
!!$       endif

!      ** ('sigma = km1 * k1ob * k1o * k1') **

!      ('grfunct = k1ob * k1o')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      k1ob, n, k1o, n, (0.d0,0.d0), grfunct, n)
       
!      ('k1o = grfunct * k1')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      grfunct, n, k1, n, (0.d0,0.d0), k1o, n)

!      ('sigma = km1 * k1o')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      km1, n, k1o, n, (0.d0,0.d0), sigma, n)

    else ! usevinv == .false.

       if (SIDE == 'L') then

!         ('sigma = km1 * al')
          call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                &
                         km1, n, al, n, (0.d0,0.d0), sigma, n)
!!$          sigma = matmul(DCONJG(TRANSPOSE(k1input)),al)

       else

!         ('sigma = km1 * ar')
          call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                &
                         km1, n, ar, n, (0.d0,0.d0), sigma, n)
!!$          sigma = matmul(DCONJG(TRANSPOSE(k1input)),ar)

       endif

    endif ! if (usevinv)

!!$    k0 = 0.5d0 * (k0input + DCONJG(TRANSPOSE(k0input)))
!!$    k1o = -k0 -sigma
    grfunct = -k0input -sigma

    ANORM = ZLANGE ('1', n, n, grfunct, n, rworkc)
    call ZGETRF (n, n, grfunct, n, IPIV, INFO)
    k1o = grfunct
    call ZGECON ('1', n, k1o, n, ANORM, RCOND, walpha, rworkc, INFO)
    call ZGETRI (n, grfunct, n, IPIV, worki, n*n, INFO)

!!$    if (.false.) then ! ignore the following code...
!!$       k1ob = 0.d0
!!$       do i = 1,n
!!$          k1ob(i,i) = 1.d0
!!$       enddo
!!$       k1o = matmul(-k0input-sigma,grfunct) - k1ob
!!$
!!$       dsigma = MAXVAL(ABS(k1o))
!!$
!!$       k1o = matmul(grfunct,k1o)
!!$       call svdmk (Ei, k1o, k1ob, n, svdtol, rnoise*svdtol,          &
!!$                   .false., .true., kappa, deltaainvn, sn)
!!$
!!$       k1o = grfunct
!!$       call svdmk (Ei, k1o, k1ob, n, svdtol, rnoise*svdtol,          &
!!$                   .false., .true., kappa, ainvn, sn)
!!$
!!$       dsigma2 = MAXVAL(ABS(grfunct))
!!$       k1o = -k0input - sigma
!!$       call svdmk (Ei, k1o, k1ob, n, svdtol, rnoise*svdtol,          &
!!$                   .false., .true., kappa, s1, sn)
!!$       ANORM2 = ZLANGE ('1', n, n, k1o, n, rworkc)
!!$
!!$       write (12347,*) "RCONDgf=", DREAL(Ei), dsigma,                &
!!$            dsigma/dsigma2, 1.d0/dsigma2,   &
!!$            rcond,1.d0/kappa, rcond*kappa, dsigma*rcond,             &
!!$            dsigma/kappa, INFO, ANORM, ANORM2,   &
!!$            s1, sn, deltaainvn/ainvn, deltaainvn, ainvn
!!$    endif

!!$    k1o = grfunct
!!$
!!$    do il = 1,10000
!!$       k1ob = k1o
!!$       call ZGETRF (n, n, k1o, n, IPIV, INFO)
!!$       call ZGETRI (n, k1o, n, IPIV, worki, n*n, INFO)
!!$       call ZGETRF (n, n, k1o, n, IPIV, INFO)
!!$       call ZGETRI (n, k1o, n, IPIV, worki, n*n, INFO)
!!$
!!$       dsigma2 = MAXVAL(ABS(k1o - k1ob))
!!$       k1ob = k1o - grfunct
!!$       dsigma = MAXVAL(ABS(k1ob))
!!$       write (*,*) "deltav=", il, dsigma, dsigma2
!!$
!!$    enddo

!!$    if (dsigma > 1.D-3) stop

!   ('k1ob = km1input * grfunct')
    call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                      &
                   km1input, n, grfunct, n, (0.d0,0.d0), k1ob, n)

!   ('k1o = k1ob * k1input')
    call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                      &
                   k1ob, n, k1input, n, (0.d0,0.d0), k1o, n)

    k1ob = sigma - k1o
    dsigma = MAXVAL(ABS(k1ob))

!!$    write (12346,*) "dsigmal=", DREAL(Ei), svdtol, dsigma

    if (dsigma < dsigmamax .or. dsigma < dsigmam                        &
         .or. dsigmam == 0.d0) then 

       calcsucceed = .true.
       nrchanm = nrchan
       if (gf == 'G') then
          matout = grfunct
       else
          matout = sigma
       endif
       if (dsigma < dsigmamax) then
          dsigmam = dsigma
          svdtolm = svdtol
       endif

    endif

    if (dsigma > dsigmamax .and. svdtol < svdtolmax) then

       write (12347,'(a," ",d," ",d," ",d)')                            &
            "rcev", DREAL(Ei), dsigma, svdtol
       if (dsigma < dsigmam .or. dsigmam == 0.d0) then
          dsigmam = dsigma
          svdtolm = svdtol
       endif
       goto 102

    elseif (svdtol /= 0.d0) then

       write (12347,'(a," ",d," ",d," ",d)')                            &
            "rcevlast", DREAL(Ei), dsigma, svdtol
       if (dsigma < dsigmam .or. dsigmam == 0.d0) then
          dsigmam = dsigma
          svdtolm = svdtol
       endif

    endif

103 if (.not. calcsucceed) then
       write (*,*) "warning: calculation of selfenergies failed"
       write (*,*) "eneside", SIDE, dreal(Ei), dsigmam
       dsigmam = 10.d0
       write (*,*) "enesideout", SIDE, dreal(Ei), dsigmam
    else
       write (12347,'(a," ",a," ",e," ",e," ",e," ",i4)')               &
            "f,sv=", SIDE, dreal(Ei), svdtolm, dsigmam, Node
       k0 = matout
    endif

!   Free memory.
    deallocate (rwork)
    deallocate (rwork3)
    deallocate (zv)
    deallocate (work)
    deallocate (work3)
    deallocate (k1o, k1ob, k1t, Aux)
    deallocate (m2s, m2sb, m2sbuf, m2sbbuf)


  end subroutine kofec


!  *******************************************************************  !
!                              kofe_svdlr2                              !
!  *******************************************************************  !
!  Description:                                                         !
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
!  character(len=1) SIDE         : 'L' (left) or 'R' (right) depending  !
!                                  on the self-energy calculated        !
!  complex*8 k0(n,n)             : 'H0 - zi*Ei*S0'                      !
!  complex*8 k1in(n,n)           : 'H1 - zi*Ei*S1'                      !
!  complex*8 km1(n,n)            : 'H1^dagger - zi*Ei*S1^dagger'        !
!  integer n                     : # of basis orbitals on the lead      !
!  complex*8 zeig(2*n)           :                                      !
!  integer neig                  :                                      !
!  complex*8 Ei                  : Energy                               !
!  real*8 svdtolz                :                                      !
!  real*8 dsigma                 :                                      !
!  ************************** INPUT/OUTPUT ***************************  !
!  complex*8 k0(n,n)          : 'H0 - zi*Ei*S0' (return self-energy)    !
!  ***************************** OUTPUT ******************************  !
!  integer nrchan             : Number of open channels                 !
!  *******************************************************************  !
  subroutine kofe_svdlr2 (SIDE, k0, k1, km1, n, neig,                   &
                          Ei, svdtol, dsigma, nrchan)


!   Input variables.
    integer, intent(in) :: n
    integer, intent(inout) :: neig
    integer, intent(out) :: nrchan
    real(8), intent(inout) :: svdtol, dsigma
    complex(8), intent(in) :: Ei
    complex(8), dimension (n,n), intent(in) :: k1, km1
    complex(8), dimension (n,n), intent(inout) :: k0
    character(len=1), intent(in) :: SIDE

!   Local variables.
    integer :: m, m2, i, INFO, nwork
    real(8) :: rmmax, svmax
    real(8), allocatable, dimension (:) :: rwork, sv, sv2
    complex(8), dimension (n,n) :: grf
    complex(8), allocatable, dimension (:) :: zv2, WORK
    complex(8), allocatable, dimension (:,:) :: k1o, mata, matb,        &
                                                matc, matd, matdi,      &
                                                v1c, v1n, vm1c, vm1n,   &
                                                v1e, vm1e, v0e,         &
                                                k0t, k1t, km1t, k0r2,   &
                                                u, vt, u2, vt2,         &
                                                Aux, AuxDag
    external :: ZGESVD, HI_zgeInvert, HI_zgemm


!   Allocate matrices and arrays.
    nwork = 8
    allocate (k1o(n,n), k0t(n,n), k1t(n,n), km1t(n,n),                  &
              u(n,n), vt(n,n), u2(n,n), vt2(n,n))
    allocate (sv(n), sv2(n))
    allocate (WORK(nwork*n))
    allocate (rwork(5*n))

!!$    km1 = DCONJG(transpose(k1))
!!$    rmmax = MAXVAL(ABS(DIMAG(k0)))
!!$    write (12347,*) "imk0=", rmmax
!!$    rmmax = MAXVAL(ABS(DIMAG(k1)))
!!$    write (12347,*) "imk1=", rmmax
!!$    rmmax = MAXVAL(ABS(DIMAG(km1)))
!!$    write (12347,*) "imkm1=", rmmax
!!$    rmmax = MAXVAL(ABS(DIMAG(k1 - DCONJG(TRANSPOSE(km1)))))
!!$    write (12347,*) "symmk1i=", rmmax
!!$    rmmax = MAXVAL(ABS(DREAL(k1 - DCONJG(TRANSPOSE(km1)))))
!!$    write (12347,*) "symmk1r=", rmmax

    k1o = k1
    if (SIDE == 'L') then
       call ZGESVD ('N', 'A', n, n, k1o, n, sv, u, n, vt, n,            &
                    WORK, nwork*n, rwork, INFO)
    else
       call ZGESVD ('A', 'N', n, n, k1o, n, sv, u, n, vt, n,            &
                    WORK, nwork*n, rwork, INFO)
    endif

    svmax = sv(1)

    m = n
    do i = 1,n
       if (sv(i)/svmax < svdtol) then
          m = i - 1
          exit
       endif
    enddo

    k1o = DCONJG(TRANSPOSE(km1))
    if (SIDE == 'L') then
       call ZGESVD ('N', 'A', n, n, k1o, n, sv2, u2, n, vt2, n,         &
                    WORK, nwork*n, rwork, INFO)
    else
       call ZGESVD ('A', 'N', n, n, k1o, n, sv2, u2, n, vt2, n,         &
                    WORK, nwork*n, rwork, INFO)
    endif

    deallocate (rwork)

    svmax = sv2(1)

    m2 = n
    do i = 1,n
       if (sv2(i)/svmax < svdtol) then
          m2 = i - 1
          exit
       endif
    enddo

    write (12347,*) "msvd2=", DREAL(Ei), m, m2, n, svdtol

    rmmax = 0.d0
    do i = 1,n
!!$       write (12347,*) "svsv2=", i, sv(i), sv2(i), sv(i)-sv2(i)
       if (abs(sv(i) - sv2(i)) > rmmax) rmmax = abs(sv(i) - sv2(i))
    enddo
    write (12347,*) "sdiffmax=", rmmax

    if (m2 < m) m = m2

    
    if (m == n) then
       write (*,*) "no need to decimate, k1 is invertable"
       allocate (v0e(n,n), v1e(n,n), vm1e(n,n))
       v0e = k0
       v1e = k1
       vm1e = km1
       goto 106
    endif

!   Allocate auxiliary matrices.
    allocate (Aux(n,n), AuxDag(n,n))

    if (SIDE == 'R') then

       allocate (v1c(m,m), v1n(m,n-m), vm1c(m,m), vm1n(n-m,m))

!      ('AuxDag = u^dagger')
       AuxDag = DCONJG(TRANSPOSE(u))

!      ('Aux = k1 * u2')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      k1, n, u2, n, (0.d0,0.d0), Aux, n)

!      ('k1t = u^dagger * Aux')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      AuxDag, n, Aux, n, (0.d0,0.d0), k1t, n)

!      ('Aux = km1 * u2')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      km1, n, u2, n, (0.d0,0.d0), Aux, n)

!      ('km1t = u^dagger * Aux')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      AuxDag, n, Aux, n, (0.d0,0.d0), km1t, n)

!      ('Aux = k0 * u2')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      k0, n, u2, n, (0.d0,0.d0), Aux, n)

!      ('k0t = u^dagger * Aux')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      AuxDag, n, Aux, n, (0.d0,0.d0), k0t, n)
 
       v1c = k1t(1:m,1:m)
       v1n = k1t(1:m,m+1:n)
       vm1c = km1t(1:m,1:m)
       vm1n = km1t(m+1:n,1:m)

    else ! SIDE == 'L'

       allocate (v1c(m,m), v1n(n-m,m), vm1c(m,m), vm1n(m,n-m))

!      ('AuxDag = vt^dagger')
       AuxDag = DCONJG(TRANSPOSE(vt))

!      ('Aux = vt2 * k1')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      vt2, n, k1, n, (0.d0,0.d0), Aux, n)

!      ('k1t = Aux * vt^dagger')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      Aux, n, AuxDag, n, (0.d0,0.d0), k1t, n)

!      ('Aux = vt2 * km1')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      vt2, n, km1, n, (0.d0,0.d0), Aux, n)

!      ('km1t = Aux * vt^dagger')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      Aux, n, AuxDag, n, (0.d0,0.d0), km1t, n)

!      ('Aux = vt2 * k0')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      vt2, n, k0, n, (0.d0,0.d0), Aux, n)

!      ('k0t = Aux * vt^dagger')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      Aux, n, AuxDag, n, (0.d0,0.d0), k0t, n)

       v1c = k1t(1:m,1:m)
       v1n = k1t(m+1:n,1:m)
       vm1c = km1t(1:m,1:m)
       vm1n = km1t(1:m,m+1:n)

    endif ! if (SIDE == 'R')

!   Free memory.
    deallocate (Aux, AuxDag)

    allocate (mata(m,m), matb(m,n-m), matc(n-m,m), matd(n-m,n-m),       &
              matdi(n-m,n-m), v1e(m,m), vm1e(m,m), v0e(m,m))
    mata = k0t(1:m,1:m)
    matb = k0t(1:m,m+1:n)
    matc = k0t(m+1:n,1:m)
    matd = k0t(m+1:n,m+1:n)
    matdi = matd

    call HI_zgeInvert (matdi, n-m)

!   Allocate auxiliary matrix.
    allocate (Aux(n-m,m))

    if (SIDE == 'R') then

!      ('Aux(n-m,m) = matdi(n-m,n-m) * matc(n-m,m)')
       call HI_zgemm ('N', 'N', n-m, m, n-m, (1.d0,0.d0),               &
                      matdi, n-m, matc, n-m, (0.d0,0.d0), Aux, n-m)

!      ('v1e(m,m) = - v1n(m,n-m) * Aux(n-m,m)')
       call HI_zgemm ('N', 'N', m, m, n-m, (-1.d0,0.d0),                &
                      v1n, m, Aux, n-m, (0.d0,0.d0), v1e, m)

       v1e = v1e + v1c

!      ('v0e(m,m) = mata(m,m) - matb(m,n-m) * Aux(n-m,m)')
       v0e = mata
       call HI_zgemm ('N', 'N', m, m, n-m, (-1.d0,0.d0),                &
                      matb, m, Aux, n-m, (1.d0,0.d0), v0e, m)

!      ('Aux(n-m,m) = matdi(n-m,n-m) * vm1n(n-m,m)')
       call HI_zgemm ('N', 'N', n-m, m, n-m, (1.d0,0.d0),               &
                      matdi, n-m, vm1n, n-m, (0.d0,0.d0), Aux, n-m)

!      ('vm1e(m,m) = - matb(m,n-m) * Aux(n-m,m)')
       call HI_zgemm ('N', 'N', m, m, n-m, (-1.d0,0.d0),                &
                      matb, m, Aux, n-m, (0.d0,0.d0), vm1e, m)

       vm1e = vm1e + vm1c

!      ('v0e(m,m) = v0e(m,m) - v1n(m,n-m) * Aux(n-m,m)')
       call HI_zgemm ('N', 'N', m, m, n-m, (-1.d0,0.d0),                &
                      v1n, m, Aux, n-m, (1.d0,0.d0), v0e, m)

    else ! SIDE == 'L'

!      ('Aux(n-m,m) = matdi(n-m,n-m) * v1n(n-m,m)')
       call HI_zgemm ('N', 'N', n-m, m, n-m, (1.d0,0.d0),               &
                      matdi, n-m, v1n, n-m, (0.d0,0.d0), Aux, n-m)

!      ('v1e(m,m) = - matb(m,n-m) * Aux(n-m,m)')
       call HI_zgemm ('N', 'N', m, m, n-m, (-1.d0,0.d0),                &
                      matb, m, Aux, n-m, (0.d0,0.d0), v1e, m)

       v1e = v1e + v1c

!      ('v0e(m,m) = mata(m,m) - vm1n(m,n-m) * Aux(n-m,m)')
       v0e = mata
       call HI_zgemm ('N', 'N', m, m, n-m, (-1.d0,0.d0),                &
                      vm1n, m, Aux, n-m, (1.d0,0.d0), v0e, m)

!      ('Aux(n-m,m) = matdi(n-m,n-m) * matc(n-m,m)')
       call HI_zgemm ('N', 'N', n-m, m, n-m, (1.d0,0.d0),               &
                      matdi, n-m, matc, n-m, (0.d0,0.d0), Aux, n-m)

!      ('vm1e(m,m) = - vm1n(m,n-m) * Aux(n-m,m)')
       call HI_zgemm ('N', 'N', m, m, n-m, (-1.d0,0.d0),                &
                      vm1n, m, Aux, n-m, (0.d0,0.d0), vm1e, m)

       vm1e = vm1e + vm1c

!      ('v0e(m,m) = v0e(m,m) - matb(m,n-m) * Aux(n-m,m)')
       call HI_zgemm ('N', 'N', m, m, n-m, (-1.d0,0.d0),                &
                      matb, m, Aux, n-m, (1.d0,0.d0), v0e, m)

    endif ! if (SIDE == 'R')

106 allocate (k0r2(m,m), zv2(2*m))

    k0r2 = v0e
!!$    write (*,*) "svmax=", svdtol, svmax, 1.d-5*svdtol*svmax
    call addnoise (v1e, m, 1.d-5*svdtol*svmax)
    call addnoise (vm1e, m, 1.d-5*svdtol*svmax)
    call addnoise (k0r2, m, 1.d-5*svdtol*svmax)
    call kofec (SIDE, 'S', k0r2, v1e, vm1e, m,                          &
                zv2, neig, Ei, dsigma, nrchan)
    deallocate (zv2)

    if (m == n) then
       k0 = k0r2

!      Free memory.
       deallocate (k0r2)
       deallocate (k1o, k0t, k1t, km1t, u, vt, u2, vt2)
       deallocate (sv, sv2)
       deallocate (WORK)
       deallocate (v0e, v1e, vm1e)

       RETURN

    endif

    if (SIDE == 'R') then

!      ('Aux(n-m,m) = matdi(n-m,n-m) * vm1n(n-m,m)')
       call HI_zgemm ('N', 'N', n-m, m, n-m, (1.d0,0.d0),               &
                      matdi, n-m, vm1n, n-m, (0.d0,0.d0), Aux, n-m)

!      ('v0e(m,m) = - v1n(m,n-m) * Aux(n-m,m)')
       call HI_zgemm ('N', 'N', m, m, n-m, (-1.d0,0.d0),                &
                      v1n, m, Aux, n-m, (0.d0,0.d0), v0e, m)

    else ! SIDE == 'L'

!      ('Aux(n-m,m) = matdi(n-m,n-m) * v1n(n-m,m)')
       call HI_zgemm ('N', 'N', n-m, m, n-m, (1.d0,0.d0),               &
                      matdi, n-m, v1n, n-m, (0.d0,0.d0), Aux, n-m)

!      ('v0e(m,m) = - vm1n(m,n-m) * Aux(n-m,m)')
       call HI_zgemm ('N', 'N', m, m, n-m, (-1.d0,0.d0),                &
                      vm1n, m, Aux, n-m, (0.d0,0.d0), v0e, m)

    endif ! if (SIDE == 'R')

    k0r2 = k0r2 + v0e

    k1o = 0.d0
    k1o(1:m,1:m) = k0r2

!   Re-allocate auxiliary matrix.
    deallocate (Aux)
    allocate (Aux(n,n))

    if (SIDE == 'R') then

!      ('Aux = u * k1o')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      u, n, k1o, n, (0.d0,0.d0), Aux, n)

!      ('k1o = Aux * u2^dagger')
       call HI_zgemm ('N', 'C', n, n, n, (1.d0,0.d0),                   &
                      Aux, n, u2, n, (0.d0,0.d0), k1o, n)

    else ! SIDE == 'L'

!      ('Aux = k1o * vt')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      k1o, n, vt, n, (0.d0,0.d0), Aux, n)

!      ('k1o = vt2^dagger * Aux')
       call HI_zgemm ('C', 'N', n, n, n, (1.d0,0.d0),                   &
                      vt2, n, Aux, n, (0.d0,0.d0), k1o, n)

    endif ! if (SIDE == 'R')

    grf = - k0 - k1o
    call HI_zgeInvert (grf, n)

    if (SIDE == 'R') then

!      ('Aux = grf * km1')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      grf, n, km1, n, (0.d0,0.d0), Aux, n)

!      ('grf = k1 * Aux')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      k1, n, Aux, n, (0.d0,0.d0), grf, n)

    else ! SIDE == 'L'

!      ('Aux = grf * k1')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      grf, n, k1, n, (0.d0,0.d0), Aux, n)

!      ('grf = km1 * Aux')
       call HI_zgemm ('N', 'N', n, n, n, (1.d0,0.d0),                   &
                      km1, n, Aux, n, (0.d0,0.d0), grf, n)

    endif

    grf = k1o - grf
    dsigma = MAXVAL(ABS(grf))
!!$    write (*,*) "dsigma_svd=", SIDE, dsigma

    k0 = k1o

!   Free memory.
    deallocate (k0r2)
    deallocate (mata, matb, matc, matd, matdi, v1e, vm1e, v0e)
    deallocate (v1c, v1n, vm1c, vm1n)
    deallocate (k1o, k0t, k1t, km1t, u, vt, u2, vt2, Aux)
    deallocate (sv, sv2)
    deallocate (WORK)


  end subroutine kofe_svdlr2


!  *******************************************************************  !
!                           get_options_kofe                            !
!  *******************************************************************  !
!  Description: Read options for 'selfenergy2' computation.             !
!                                                                       !
!  From SMEAGOL code - 2003.                                            !
!                                                                       !
!  Written by Alexandre Reily Rocha, Jun 2003.                          !
!  Computational Spintronics Group                                      !
!  Trinity College Dublin                                               !
!  e-mail: reilya@ift.unesp.br                                          !
!  ***************************** HISTORY *****************************  !
!  Original version:    June 2003                                       !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode                : True if it is the I/O node           !
!  integer MPI_Comm_MyWorld      : MPI communicator                     !
!  ***************************** OUTPUT ******************************  !
!  logical usehinv               :                                      !
!  logical usevinv               :                                      !
!  real*8 tolki                  :                                      !
!  real*8 svdtolmax              :                                      !
!  real*8 svdtolmin              :                                      !
!  real*8 dsigmamax              :                                      !
!  logical callsvd               :                                      !
!  real*8 rnoise                 :                                      !
!  real*8 skipsvd                :                                      !
!  *******************************************************************  !
  subroutine get_options_kofe (usehinv, usevinv, tolki,                 &
                               svdtolmax, svdtolmin, dsigmamax,         &
                               callsvd, rnoise, skipsvd)

!
!   Modules
!
#ifdef MPI
    use parallel,        only: IOnode, MPI_Comm_MyWorld
#else
    use parallel,        only: IOnode
#endif
    use fdf

#ifdef MPI
    include "mpif.h"
#endif

!   Input variables.
    real(8), intent(out) :: tolki, svdtolmax, svdtolmin,                &
                            dsigmamax, rnoise, skipsvd
    logical, intent(out) :: usehinv, usevinv, callsvd
      
!   Local variables.
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

    if (IOnode) then
       usehinv = fdf_boolean ('Sigma.InvertH1', .false.)
       usevinv = fdf_boolean ('Sigma.InvertV', .false.)
       callsvd = fdf_boolean ('Sigma.CSVD', .false.)
       tolki = fdf_double ('Sigma.Dkimag', 1.d-6)
       svdtolmax = fdf_double ('Sigma.DSVDMax', 5.d-9)
       svdtolmin = fdf_double ('Sigma.DSVDMin', 1.d-15)
       dsigmamax = fdf_double ('Sigma.DMax', 1.d-5)
       rnoise = fdf_double ('Sigma.RNoise', 1.d0)
       skipsvd = fdf_double ('Sigma.SkipSVD', 1.d3)
       if (skipsvd < 2.d0) skipsvd = 10.d0
    endif

#ifdef MPI
    call MPI_Bcast (usehinv, 1, MPI_Logical, 0,                         &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (usevinv, 1, MPI_Logical, 0,                         &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (callsvd, 1, MPI_Logical, 0,                         &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (tolki, 1, MPI_Double_Precision, 0,                  &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (svdtolmax, 1, MPI_Double_Precision, 0,              &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (svdtolmin, 1, MPI_Double_Precision, 0,              &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (dsigmamax, 1, MPI_Double_Precision, 0,              &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (rnoise, 1, MPI_Double_Precision, 0,                 &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (skipsvd, 1, MPI_Double_Precision, 0,                &
                    MPI_Comm_MyWorld, MPIerror)
#endif

!!$      write (*,*) "Sigma.Options =", usehinv, usevinv, tolki,        &
!!$           svdtolmax, svdtolmin, dsigmamax, callsvd


    end subroutine get_options_kofe


!  *******************************************************************  !
!                                 svdm                                  !
!  *******************************************************************  !
!  Description:                                                         !
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
!  integer n                     :                                      !
!  real*8 svdtol                 :                                      !
!  real*8 noisemax               :                                      !
!  logical fmatinv               :                                      !
!  logical faddnoise             :                                      !
!  ***************************** OUTPUT ******************************  !
!  complex*8 matin(n,n)          :                                      !
!  complex*8 invmat(n,n)         :                                      !
!  *******************************************************************  !
    subroutine svdm (matin, invmat, n, svdtol,                          &
                     noisemax, fmatinv, faddnoise)


!   Input variables.
    integer, intent(in) :: n
    real(8), intent(in) :: svdtol, noisemax
    complex(8), dimension (n,n), intent(out) :: matin, invmat
    logical, intent(in) :: fmatinv, faddnoise

!   Local variables.
    integer :: m, INFO, i, j, i1
    real(8) :: svmax, condn
    real(8), dimension (n) :: sv
    real(8), allocatable, dimension (:) :: rwork
    complex(8), parameter :: alpha = (1.d0,0.d0)
    complex(8), allocatable, dimension (:) :: work
    complex(8), dimension (n,n) :: u, vt
    external :: ZGESVD


!   Allocate auxiliary arrays.
    allocate (work(8*n))
    allocate (rwork(5*n))


    call ZGESVD ('A', 'A', n, n, matin, n, sv, u, n, vt, n,             &
                 work, 8*n, rwork, INFO)
    if (info /= 0) write(*,*) "infosvd=", info

    svmax = sv(1)
    m = n
    do i = 1,n
       condn = svmax/sv(i)
!!$       write(*,*) "sigmak1=", DREAL(Ei), i, sv(i), sv(i)/svmax, condn
       if (sv(i)/svmax < svdtol) then
          m = i - 1
          exit
       endif
    enddo

!!$    write(12347,*) "m,svmax=", m, svmax, svdtol
!!$    write(*,*) "m,svmax=", m, svmax, svdtol

    do i = 1,n
       if (sv(i)/svmax < svdtol) then
          sv(i) = svdtol * svmax 
       endif
!!$       sv(i) = sv(i) + svdtol * svmax
!!$       if (sv(i) < 1.d-15) then
!!$          sv(i) = sv(i) / 1.d0
!!$       endif
    enddo

    matin = 0.d0
    do i = 1,n
       do j = 1,n
          do i1 = 1,n
             matin(i,j) = matin(i,j) + u(i,i1) * sv(i1)  * vt(i1,j)
          enddo
       enddo
    enddo

!!$    k1o = 0.d0
!!$    do i = 1,n
!!$       do j = 1,n
!!$          do i1 = 1,n
!!$             k1o(i,j) = k1o(i,j) +                                   &
!!$                  DCONJG(vt(i1,i)) *  DCONJG(u(j,i1)) / sv(i1) 
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$    invmat = 0.d0
!!$    do i = 1,n
!!$       invmat(i,i) = 1.d0 / sv(i)
!!$    enddo
!!$    k1ob = matmul(invmat,k1o)
!!$
!!$    call geigenvalues2 (k1ob, zv, n, k1o, info)
!!$
!!$    do i = 1,n
!!$       call alphaofz (alphak, lz, zv(i))
!!$       write(*,*) "kapprox=", DREAL(Ei), lz, alphak, svdtol
!!$    enddo

    if (fmatinv) then
       invmat = 0.d0
       do i = 1,n
          do j = 1,n
             do i1 = 1,n
                invmat(i,j) = invmat(i,j) +                             &
                     DCONJG(vt(i1,i)) *  DCONJG(u(j,i1)) / sv(i1) 
             enddo
          enddo
       enddo
    endif

    if (faddnoise) then
       call addnoise (matin, n, noisemax)
       call addnoise (invmat, n, noisemax)
    endif

!   Free memory.
    deallocate (rwork)
    deallocate (work)


  end subroutine svdm


!  *******************************************************************  !
!                               addnoise                                !
!  *******************************************************************  !
!  Description:                                                         !
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
!  integer n                     :                                      !
!  real*8 delta                  :                                      !
!  ***************************** OUTPUT ******************************  !
!  complex*8 mat(n,n)            :                                      !
!  *******************************************************************  !
  subroutine addnoise (mat, n, delta)


!   Input variables.
    integer, intent(in) :: n
    real(8), intent(in) :: delta
    complex(8), dimension (n,n), intent(out) :: mat

!   Local variables.
    integer :: i, j
    real(8) :: rnum, inum
    complex(8) :: noise
    logical, save :: initseed = .true.

    if (initseed) then
       write (12347,*) "initseed"
       call random_seed
       initseed = .false.
    endif

    do i = 1,n
       do j = 1,n
          call random_number (rnum)
          call random_number (inum)
!!$          write (*,*) "rnum=", rnum
          noise = (rnum - 0.5d0 + (0.d0,1.d0) * (inum - 0.5d0)) * delta
          mat(i,j) = mat(i,j) + noise
       enddo
    enddo


  end subroutine addnoise


!  *******************************************************************  !
!                             geigenvalues                              !
!  *******************************************************************  !
!  Description: call lapack 'zggev' for computing the generalized       !
!  eigenvalues for a pair of nonsymmetric matrices.                     !
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
!  complex*8 mat(n,n)            :                                      !
!  complex*8 mat2(n,n)           :                                      !
!  integer n                     :                                      !
!  ***************************** OUTPUT ******************************  !
!  complex*8 walpha(n)           :                                      !
!  complex*8 wbeta(n)            :                                      !
!  complex*8 vrg(n,n)            :                                      !
!  integer info                  :                                      !
!  *******************************************************************  !
  subroutine geigenvalues (mat, mat2, walpha, wbeta, n, vrg, info)


!   Input variables.
    integer, intent(in) :: n
    integer, intent(out) :: info
    complex(8), dimension (n), intent(out) :: walpha, wbeta
    complex(8), dimension (n,n), intent(in) :: mat, mat2
    complex(8), dimension (n,n), intent(out) :: vrg

!   Local variables.
    real(8), allocatable, dimension (:) :: rwork
    complex(8), allocatable, dimension (:) :: work
    complex(8), allocatable, dimension (:,:) :: matbuf, mat2buf, vlg
    external :: zggev

!   Allocate auxiliary array and matrices.
    allocate (rwork(8*n))
    allocate (work(4*n))
    allocate (matbuf(n,n), mat2buf(n,n), vlg(n,n))

    matbuf = mat
    mat2buf = mat2
    call zggev ('N', 'V', n, matbuf ,n, mat2buf, n, walpha, wbeta,      &
                vlg, n, vrg, n, work, 4*n, rwork, info)
    if (info /= 0) write (*,*) "info_c2=", info

!   Free memory.
    deallocate (rwork)
    deallocate (work)
    deallocate (matbuf, mat2buf, vlg)


  end subroutine geigenvalues


!  *******************************************************************  !
!                             geigenvalues2                             !
!  *******************************************************************  !
!  Description: call lapack 'zgeevx' for computing the eigenvalues and  !
!  left and right eigenvectors of a general matrix, with preliminary    !
!  matrix balancing, and computes reciprocal condition numbers for the  !
!  eigenvalues and right eigenvectors.                                  !
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
!  complex*8 mat(n,n)            :                                      !
!  integer n                     :                                      !
!  ***************************** OUTPUT ******************************  !
!  complex*8 zv(n)               :                                      !
!  complex*8 vrg(n,n)            :                                      !
!  integer info                  :                                      !
!  *******************************************************************  !
  subroutine geigenvalues2 (mat, zv, n, vrg, info)


!   Input variables.
    integer, intent(in) :: n
    integer, intent(out) :: info
    complex(8), dimension (n), intent(out) :: zv
    complex(8), dimension (n,n), intent(in) :: mat
    complex(8), dimension (n,n), intent(out) :: vrg

!   Local variables.
    integer ilo, ihi
    real(8) :: ABNRM
    real(8), allocatable, dimension (:) :: scaleev, rconde, rcondv, rwork
    complex(8), allocatable, dimension (:) :: work
    complex(8), allocatable, dimension (:,:) :: matbuf, vlg
    external :: zgeevx

!   Allocate auxiliary array and matrices.
    allocate (scaleev(n), rconde(n), RCONDV(n))
    allocate (rwork(2*n))
    allocate (work(4*n))
    allocate (matbuf(n,n), vlg(n,n))
      
    matbuf = mat
    call zgeevx ('P', 'N', 'V', 'N', N, matbuf, N, zv,                  &
                 vlg, N, vrg, N, ilo, ihi, scaleev, abnrm,              &
                 rconde, rcondv, work, 4*n, rwork, info)

!   Free memory.
    deallocate (scaleev, rconde, rcondv)
    deallocate (rwork)
    deallocate (work)
    deallocate (matbuf, vlg)


  end subroutine geigenvalues2


!  *******************************************************************  !
!                               alphaofz                                !
!  *******************************************************************  !
!  Description:                                                         !
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
!  real*8 zev                    :                                      !
!  ***************************** OUTPUT ******************************  !
!  complex*8 alpha               :                                      !
!  complex*8 kappa               :                                      !
!  *******************************************************************  !
  subroutine alphaofz (alpha, kappa, zev)


!   Input variables.
    real(8), intent(out) :: alpha, kappa
    complex(8), intent(in) :: zev

!   Local variables.
    real(8) :: lz, zr, zi, alpha0
    real(8), parameter :: pi = 3.1415926535897932384626433832795028841D0


!   Initialize variables.
    lz = CDABS(zev)
    zr = DREAL(zev)
    zi = DIMAG(zev)

    if (zr >= 0.d0 .and. zi >= 0.d0) then
       alpha0 = 0.d0
    elseif (zr >= 0.d0 .and. zi < 0.d0) then
       zr = -DIMAG(zev)
       zi = DREAL(zev)
       alpha0 = - pi / 2.d0
    elseif (zr < 0.d0 .and. zi >= 0.d0) then
       zr = DIMAG(zev)
       zi = -DREAL(zev)
       alpha0 = pi / 2.d0
    elseif (zr < 0.d0 .and. zi < 0.d0) then
       zr = -DREAL(zev)
       zi = -DIMAG(zev)
       alpha0 = -pi
    endif

    if (zr == 0.d0) then
       alpha = pi / 2.d0 + alpha0
    else
       alpha = ATAN(zi / zr) + alpha0
!!$       alpha = ATAN((zi / lz) / (zr / lz)) + alpha0
    endif
    kappa = DLOG(lz)


  end subroutine alphaofz


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
