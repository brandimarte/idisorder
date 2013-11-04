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
!                        MODULE idsrdr_ephcoupl                         !
!  *******************************************************************  !
!  Description: reads the electron-phonon coupling matrix, the          !
!  vibrational frequencies (which is given in eV) and the dynamic       !
!  orbitals indexes, for each deffect unit.                             !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *******************************************************************  !

MODULE idsrdr_ephcoupl

!
!   Modules
!
  use parallel,        only: 
  use idsrdr_options,  only: 

  implicit none

  PRIVATE ! default is private
  PUBLIC :: eph, neph, nModes, norbDyn, idxF, idxL, freq, Meph,         &
            EPHread, EPHfree

  logical :: eph ! Inelastic calculation?

  integer :: neph ! # of units with e-ph interaction
  integer, dimension(:), allocatable :: nModes ! # of vibrational modes
  integer, dimension(:), allocatable :: norbDyn ! # of orbitals from
                                                ! dynamic atoms
  integer, dimension(:), allocatable :: idxF ! First dynamic atom orbital
  integer, dimension(:), allocatable :: idxL ! Last dynamic atom orbital

  TYPE ephFreq
     real(8), pointer :: F(:) ! pointer to frequency array
  END TYPE ephFreq

  TYPE(ephFreq), dimension(:), allocatable :: freq ! Vibrational mode's
                                                   ! frequencies (energy)

  TYPE ephCplng
     complex(8), pointer :: M(:,:,:,:) ! pointer to eph coupling matrices
  END TYPE ephCplng

  TYPE(ephCplng), dimension(:), allocatable :: Meph ! Electron-phonon
                                                    ! coupling matrix


CONTAINS


!  *******************************************************************  !
!                                EPHread                                !
!  *******************************************************************  !
!  Description: reads the electron-phonon coupling matrix, the          !
!  vibrational frequencies (which is given in eV) and the dynamic       !
!  orbitals indexes.                                                    !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2013                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode               : True if it is the I/O node            !
!  integer nspin                : Number of spin components             !
!  character(60) directory      : Working directory                     !
!  ****************************** INPUT ******************************  !
!  integer ntypeunits                    : Number of unit types         !
!  character(30) fileunits(ntypeunits+2) : Units files                  !
!  integer unitdimensions(ntypeunits+2)  : Units number of orbitals     !
!  integer ephIndic(ntypeunits+2)        : E-ph interaction indicator   !
!  ***************************** OUTPUT ******************************  !
!  logical eph                 : Inelastic calculation?
!  integer neph                : Number of units with e-ph interaction  !
!  integer nModes(neph)        : Number of vibrational modes            !
!  integer norbDyn(neph)       : Number of orbitals from dynamic atoms  !
!  integer idxF(neph)          : First dynamic atom orbital             !
!  integer idxL(neph)          : Last dynamic atom orbital              !
!  TYPE(ephFreq) freq(neph)%F(nModes) : [real*8] Vibrational mode's     !
!                                       frequencies                     !
!  TYPE(ephCplng) Meph(neph)%M(norbDyn,norbDyn,nModes*nspin) :          !
!                                   [complex*8] E-ph coupling matrices  !
!  *******************************************************************  !
  subroutine EPHread (ntypeunits, fileunits, unitdimensions, ephIndic)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_options,  only: nspin, directory

#ifdef MPI
    include "mpif.h"
#endif

!   Input variables.
    integer, intent(in) :: ntypeunits
    integer, dimension (ntypeunits+2), intent(in) :: unitdimensions
    integer, dimension (ntypeunits+2), intent(in) :: ephIndic
    character(len=30), dimension (ntypeunits+2), intent(in) :: fileunits

!   Local variables.
    integer :: iu, idx, nu, i, l, s, sDyn, nDyn
    complex(8), dimension(:), allocatable :: aux
    character(len=100), external :: paste
    logical :: found
    external :: io_assign, io_close
#ifdef MPI
    integer :: MPIerror
#endif

!   Number of units with e-ph interaction.
    if (IOnode) then
       neph = SUM (ephIndic)
       write (6,'(/)',advance='no')
    endif

#ifdef MPI
    call MPI_Bcast (neph, 1, MPI_Integer, 0, MPI_Comm_World, MPIerror)
#endif

    if (neph == 0) then
       if (IOnode) write (6,'(/,a)')                                    &
            'EPHread: Calculation without e-ph interaction!'
       return
    endif

!   Allocate arrays.
    allocate (nModes(neph))
    allocate (norbDyn(neph))
    allocate (idxF(neph))
    allocate (idxL(neph))
    allocate (freq(neph))
    allocate (Meph(neph))

!   Read e-ph interaction data of the required units.
    idx = 1
    do nu = 1,ntypeunits+2

       if (ephIndic(nu) == 1) then ! consider eph interaction

!         Check if electron-phonon coupling file exists.
          inquire (file=paste(directory, paste(fileunits(nu),'.Meph')), &
               exist=found)
          if (.not.found) go to 123

          if (IOnode) then

             write(6,'(a,i3,a)',advance='no') "EPHread: Reading"   //   &
                  " electron-phonon coupling matrix ", idx,  "..."

!            Opens the file.
             call io_assign (iu)
             open (iu,                                                  &
                  file=paste(directory, paste(fileunits(nu),'.Meph')),  &
                  form='formatted', status='old')

!            Reads: # of spin (sDyn), # of dynamic atoms (nDyn), # of
!            dynamic orbitals (norbDyn), first orbital index (idxF) and
!            last orbital index (idxL).
             read (iu,*) sDyn, nDyn, norbDyn(idx), idxF(idx), idxL(idx)

!            Verifies spin number.
             if (sDyn /= nspin) then
                write (6,'(/,a,/)')                                     &
                     "EPHread: ERROR: Electron-phonon coupling" //      &
                     " calculated with different spin number."
#ifdef MPI
                call MPI_Abort (MPI_Comm_World, 1, MPIerror)
                stop
#else
                stop
#endif
             endif

!            Verifies orbital indexes.
             if (idxL(idx) > unitdimensions(nu)) then
                write (6,'(/,a,/)')                                     &
                     "EPHread: ERROR: Last orbital index from"     //   &
                     " electron-phonon coupling matrix is greater" //   &
                     " than the total number of orbitals."
#ifdef MPI
                call MPI_Abort (MPI_Comm_World, 1, MPIerror)
                stop
#else
                stop
#endif
             endif

!            Total number of modes.
             nModes(idx) = 3 * nDyn

!            Allocates memory.
             allocate (freq(idx)%F(nModes(idx)))
             allocate (Meph(idx)%M(norbDyn(idx),norbDyn(idx),           &
                                   nspin,nModes(idx)))
             allocate (aux(norbDyn(idx)))

!            Reads the mode's frequencies (energies).
             read (iu,*) freq(idx)%F(1:nModes(idx))
             do l = 1,3*nDyn ! don't consider null modes
                if (freq(idx)%F(l) <= 0.D0) nModes(idx) = nModes(idx) - 1
             enddo

!            Reads the electron-phonon coupling matrix.
             do l = 1,nModes(idx)
                do s = 1,nspin
                   do i = 1,norbDyn(idx)
                      read (iu,*) aux
                      Meph(idx)%M(i,1:norbDyn(idx),s,l) = aux
                   enddo
                enddo
             enddo

!            Energies in Ry (from CODATA - 2012).
             freq(idx)%F = freq(idx)%F / 13.60569253D0
             Meph(idx)%M = Meph(idx)%M / 13.60569253D0

!            Closes the file.
             call io_close (iu)

!            Free memory.
             deallocate (aux)

             write(6,'(a)') " ok!"

          endif ! if (IOnode)

          idx = idx + 1

       endif ! if (ephIndic == 1)

    enddo ! do nu = 1,ntypeunits+2

!   Broadcast read variables.
#ifdef MPI
    call MPI_Bcast (nModes, neph, MPI_Integer, 0,                       &
                    MPI_Comm_World, MPIerror)
    call MPI_Bcast (norbDyn, neph, MPI_Integer, 0,                      &
                    MPI_Comm_World, MPIerror)
    call MPI_Bcast (idxF, neph, MPI_Integer, 0,                         &
                    MPI_Comm_World, MPIerror)
    call MPI_Bcast (idxL, neph, MPI_Integer, 0,                         &
                    MPI_Comm_World, MPIerror)
    if (.not. IOnode) then
       do idx = 1,neph
          allocate (freq(idx)%F(nModes(idx)))
          allocate (Meph(idx)%M(norbDyn(idx),norbDyn(idx),              &
                                nspin,nModes(idx)))
       enddo
    endif
    do idx = 1,neph
       call MPI_Bcast (freq(idx)%F, nModes(idx), MPI_Double_Precision,  &
                       0, MPI_Comm_World, MPIerror)
       call MPI_Bcast (Meph(idx)%M, norbDyn(idx)*norbDyn(idx)           &
                       *nspin*nModes(idx), MPI_Double_Complex, 0,       &
                       MPI_Comm_World, MPIerror)
    enddo
#endif

    RETURN

123 if (IOnode) then
       write (6,'(/,a,/)') "EPHread: ERROR: Couldn't open " //          &
            "electron-phonon coupling file."
#ifdef MPI
       call MPI_Abort (MPI_Comm_World, 1, MPIerror)
       stop
#else
       stop
#endif
    else
       RETURN
    endif


  end subroutine EPHread


!  *******************************************************************  !
!                                EPHfree                                !
!  *******************************************************************  !
!  Description: free allocated vectors.                                 !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *******************************************************************  !
  subroutine EPHfree

!   Local variables.
    integer :: I

!   Free memory.
    if (eph) then
       deallocate (nModes)
       deallocate (norbDyn)
       deallocate (idxF)
       deallocate (idxL)
!      First deallocates pointed arrays and matrices.
       do I = 1,neph
          deallocate (freq(I)%F)
          deallocate (Meph(I)%M)
       enddo
       deallocate (freq)
       deallocate (Meph)
    endif


  end subroutine EPHfree


!  *******************************************************************  !


END MODULE idsrdr_ephcoupl

