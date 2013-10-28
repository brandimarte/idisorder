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

  implicit none
  
  PUBLIC ! default is public

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
    use parallel,        only: IOnode
    use idsrdr_options,  only: nspin, temp, directory

    include "mpif.h"

!   Input variables.
    integer, intent(in) :: nsc(2)

!   Local variables.
    integer :: iu, nspinu, maxnh
    real(8) :: EfLeadR
    character (len=30) :: slabeli
    character(len=72), external :: paste
    external :: io_assign, io_close, zhsunits, ranksvd
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

    if (IOnode) then

       write (6,'(/,30("*"),a,31("*"))')                              &
            ' Leads input data '

       call io_assign (iu)
       open (iu, file=paste(directory,'bulklft.DAT'), status='old')
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
       call io_close (iu)
       call io_assign (iu)
       open (iu, file=paste(directory,'bulkrgt.DAT'), status='old')
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
       write (6,6)                                                      &
            'readleads: Lead Fermi energy                       ' //    &
            '      =', EfLead, ' Ry'
       call io_close (iu)

       write (6,'(2a)') 'readleads: ', repeat('*', 68)

    endif

#ifdef MPI
    call MPI_Bcast (NL, 1, MPI_Integer, 0, MPI_Comm_world, MPIerror)
    call MPI_Bcast (NR, 1, MPI_Integer, 0, MPI_Comm_world, MPIerror)
    call MPI_Bcast (EfLead, 1 ,MPI_Double_Precision, 0,                 &
                    MPI_Comm_world, MPIerror)
#endif

!   Allocate hamiltonians and overlap matrices.
    allocate (S0_L(NL,NL), S1_L(NL,NL))
    allocate (H0_L(NL,NL,nspin), H1_L(NL,NL,nspin))
    allocate (S0_R(NR,NR), S1_R(NR,NR))
    allocate (H0_R(NR,NR,nspin), H1_R(NR,NR,nspin))

!   Reads and calculates the Hamiltonians
!   and overlaps matrices of the leads.
    if (IOnode) then
       call zhsunits (nspin, nspin, NL, nsc, 1, 0, 0, .true., 1, temp,  &
                      H0_L, H1_L, S0_L, S1_L, paste(directory,'bulklft'))
       call zhsunits (nspin, nspin, NR, nsc, 1, 0, 0, .true., 1, temp,  &
                      H0_R, H1_R, S0_R, S1_R, paste(directory,'bulkrgt'))
    endif

#ifdef MPI
    call MPI_Bcast (H0_L, NL*NL*nspin, MPI_Double_Complex, 0,           &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (S0_L, NL*NL, MPI_Double_Complex, 0,                 &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (H1_L, NL*NL*nspin, MPI_Double_Complex, 0,           &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (S1_L, NL*NL, MPI_Double_Complex, 0,                 &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (H0_R, NR*NR*nspin, MPI_Double_Complex, 0,           &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (S0_R, NR*NR, MPI_Double_Complex, 0,                 &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (H1_R, NR*NR*nspin, MPI_Double_Complex, 0,           &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (S1_R, NR*NR, MPI_Double_Complex, 0,                 &
                    MPI_Comm_world, MPIerror)
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
!  integer INFO                   : Is '0' if self-energy calculated    !
!                                   with success, 'NOT 0' otherwise     !
!  integer NCHAN                  : Number of open channels             !
!  *******************************************************************  !
  subroutine leadsSelfEn (Ei, ispin, INFO, NCHAN)


!   Input variables.
    integer, intent(in) :: ispin
    integer, intent(out) :: INFO
    integer, intent(out) :: NCHAN
    real(8), intent(in) :: Ei

!   Local variables.
    external :: SELFENERGY

    call SELFENERGY ('L', side_rankL(ispin), numberL(ispin), NL,        &
                     DCMPLX(Ei), H0_L(:,:,ispin), VHL(:,:,ispin),       &
                     S0_L, VSL(:,:,ispin), Sigma_L, QL(:,:,ispin),      &
                     INFO, NCHAN)

    if (info == 0) then
       call SELFENERGY ('R', side_rankR(ispin), numberR(ispin), NR,     &
                        DCMPLX(Ei), H0_R(:,:,ispin), VHR(:,:,ispin),    &
                        S0_R, VSR(:,:,ispin), Sigma_R, QR(:,:,ispin),   &
                        INFO, NCHAN)
    endif


  end subroutine leadsSelfEn


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
