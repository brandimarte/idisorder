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
!                          MODULE idsrdr_units                          !
!  *******************************************************************  !
!  Description: read and build disorder units (layers).                 !
!                                                                       !
!  Use FDF (Flexible Data Format) package of J.M.Soler and A.Garcia.    !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *******************************************************************  !

MODULE idsrdr_units

!
!   Modules
!
  use parallel,        only: 
  use idsrdr_options,  only: 
  use idsrdr_leads,    only: 
  use idsrdr_init,     only: 
  use idsrdr_ephcoupl, only: 
  use fdf

  implicit none
  
  PUBLIC  :: N1, unit_type, unitdimensions, unitshift, theta,           &
             S1unit, H1unit, Sunits, Hunits, makeunits, freeunits,      &
             nunitseph, eph_type, ephIndic
  PRIVATE :: readunits, TBreadunits, buildunits

  integer :: N1 ! Leads total number of orbitals (left + right)
  integer :: nunitseph ! Number of units with eph
  integer, allocatable, dimension (:) :: unit_type ! Units types
  integer, allocatable, dimension (:) :: unitdimensions ! Units number
                                                        ! of orbitals
  integer, allocatable, dimension (:) :: ephIndic ! E-ph interaction in
                                                  ! unit type (0 or 1)
  integer, allocatable, dimension (:) :: eph_type ! Units types with eph

  real(8), allocatable, dimension (:) :: theta ! 
  real(8), allocatable, dimension (:) :: unitshift ! Units shift
  real(8), allocatable, dimension (:,:) :: S1unit ! Unit coupling overlap
  real(8), allocatable, dimension (:,:,:) :: H1unit ! Unit coupling
                                                    ! Hamiltonian

  TYPE unitS
     real(8), pointer :: S(:,:) ! pointer to units overlaps
  END TYPE unitS

  TYPE(unitS), allocatable, dimension (:) :: Sunits ! Units overlap

  TYPE unitH
     real(8), pointer :: H(:,:,:) ! pointer to units hamiltonians
  END TYPE unitH

  TYPE(unitH), allocatable, dimension (:) :: Hunits ! Units hamiltonian

CONTAINS


!  *******************************************************************  !
!                               makeunits                               !
!  *******************************************************************  !
!  Description: main subroutine that computes a random distribution of  !
!  the deffects and, read and build the units (size, weight, shift,     !
!  overlap and hamiltonian).                                            !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode               : True if it is the I/O node            !
!  integer NDeffects            : Number of deffects blocks             !
!  integer ntypeunits           : Number of unit types                  !
!  real*8 avgdist               : Average deffect distance              !
!  integer nspin                : Number of spin components             !
!  real*8 temp                  : Electronic temperature                !
!  character(60) directory      : Working directory                     !
!  integer NL                   : Number of left lead orbitals          !
!  integer NR                   : Number of right lead orbitals         !
!  integer nsc(2)               : Number of unit cells along parallel   !
!                                 directions                            !
!  ***************************** OUTPUT ******************************  !
!  integer N1                           : Leads total number of         !
!                                         orbitals (left + right)       !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  real*8 theta(NDeffects+1)            :                               !
!  real*8 unitshift(ntypeunits+2)       : Units shift                   !
!  *******************************************************************  !
  subroutine makeunits

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_options,  only: NDeffects, ntypeunits, avgdist,          &
                               nspin, temp, directory
    use idsrdr_leads,    only: NL, NR
    use idsrdr_init,     only: nsc

!   Local variables.
    real(8), allocatable, dimension (:) :: dist ! Deffects distribution
    real(8), allocatable, dimension (:) :: unitlength ! Units size (in z)
    real(8), allocatable, dimension (:) :: unitweight ! Units weight
    character (len=30), dimension (:), allocatable :: fileunits ! Units
                                                                ! files
    external :: randon_d

    if (IOnode) write (6,'(/,25("*"),a,26("*"),/)')                     &
         ' Reading and Building units '

!   Allocate arrays.
    allocate (dist(NDeffects+1))
    allocate (theta(NDeffects+1))
    allocate (fileunits(ntypeunits+2))
    allocate (unitlength(ntypeunits+2))
    allocate (unitshift(ntypeunits+2))
    allocate (unitweight(ntypeunits+2))
    allocate (unitdimensions(ntypeunits+2))

!   Initialize arrays.
    dist = 0.d0
    theta = 0.d0
    fileunits = ""
    unitlength = 0.d0
    unitshift = 0.d0 
    unitweight = 0.d0
    unitdimensions = 0
    unitdimensions(ntypeunits+1) = NL
    unitdimensions(ntypeunits+2) = NR

    call random_d (Ndeffects+1, avgdist, dist, theta)

    if (IOnode) write (6,'(a,f15.10,/)')                                &
         'makeunits: temperature before = ', temp
    call readunits (nspin, ntypeunits, nsc, temp, unitdimensions,       &
                    unitlength, unitshift, unitweight, fileunits,       &
                    directory) 
    if (IOnode) write (6,'(a,f15.10,/)')                                &
         'makeunits: temperature after = ', temp

!   Leads total number of orbitals (left + right).
    N1 = unitdimensions(ntypeunits+1) + unitdimensions(ntypeunits+2)
   
    call buildunits (NDeffects, ntypeunits, unitlength, unitweight, dist)

!   Free memory.
    deallocate (dist)
    deallocate (fileunits)
    deallocate (unitlength)
    deallocate (unitweight)

    if (IOnode) write (6,'(/,2a)') 'makeunits: ', repeat('*', 68)


  end subroutine makeunits


!  *******************************************************************  !
!                               readunits                               !
!  *******************************************************************  !
!  Description: Reads the units from files and calculates their         !
!  overlap and hamiltonian matrices.                                    !
!                                                                       !
!  Use FDF (Flexible Data Format) package of J.M.Soler and A.Garcia.    !
!                                                                       !
!  Written by Alexandre Reily Rocha, --- 2007.                          !
!  Instituto de Fisica Teorica                                          !
!  Universidade Estadual de Sao Paulo                                   !
!  e-mail: reilya@ift.unesp.br                                          !
!  ***************************** HISTORY *****************************  !
!  Original version:    --- 2007                                        !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode               : True if it is the I/O node            !
!  logical tightbinding         : Tight-binding calculation?            !
!  ****************************** INPUT ******************************  !
!  integer nspin                : Number of spin components             !
!  integer ntypeunits           : Number of unit types                  !
!  integer nsc(2)               : Number of unit cells along parallel   !
!                                 directions                            !
!  real*8 temp                  : Electronic temperature                !
!  character(60) directory      : Working directory                     !
!  ************************** INPUT/OUTPUT ***************************  !
!  integer unitdimensions(ntypeunits+2)  : Units number of orbitals     !
!  real*8 unitlength(ntypeunits+2)       : Units size (in z direction)  !
!  real*8 unitshift(ntypeunits+2)        : Units shift                  !
!  real*8 unitweight(ntypeunits+2)       : Units weight                 !
!  character(30) fileunits(ntypeunits+2) : Units files                  !
!  ***************************** OUTPUT ******************************  !
!  integer ephIndic(ntypeunits+2)        : E-ph interaction in          !
!                                          unit type (0 or 1)           !
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
!  *******************************************************************  !
  subroutine readunits (nspin, ntypeunits, nsc, temp,                   &
                        unitdimensions, unitlength, unitshift,          &
                        unitweight, fileunits, directory)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_ephcoupl, only: EPHread
    use idsrdr_options,  only: tightbinding
    use fdf

#ifdef MPI
    include "mpif.h"
#endif

!   Input variables.
    integer, intent(in) :: nspin, ntypeunits
    integer, dimension (2), intent(in) :: nsc
    integer, dimension (ntypeunits+2), intent(out) :: unitdimensions
    real(8), intent(in) :: temp
    real(8), dimension (ntypeunits+2), intent(out) :: unitlength,       &
                                                      unitshift,        &
                                                      unitweight
    character(len=30), dimension (ntypeunits+2), intent(out) :: fileunits
    character(len=60), intent(in) :: directory

!   Local variables.
    integer :: iu, I, nspinu, no, nuo, maxnh
    integer, dimension (2) :: nscu
    real(8) :: tempu, efu
    real(8), allocatable, dimension (:,:,:) :: H0aux, H1aux
    real(8), allocatable, dimension (:,:) :: S0aux, S1aux
    logical :: therearefiles
    character(len=30) :: slabel
    character(len=100), external :: paste
    external :: io_assign, io_close, hsunits
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

!   Allocate electron-phonon interaction indicator array.
    allocate (ephIndic(ntypeunits+2))

!   Initialize arrays.  
    unitlength = 0.d0
    unitshift = 0.d0
    ephIndic = 0

!   Read 'UnitsFiles' block (length, shift, weight and file name).
    if (IOnode) then
       therearefiles = fdf_block ('UnitsFiles', iu)
       If (therearefiles) Then

          do I = 1,ntypeunits+2
             read (iu,*) unitlength(I), unitshift(I),                   &
                  unitweight(I), ephIndic(I), fileunits(I) 
             write (6,'(a,i3,a,a)') 'readunits: Unit ', I, ' - ',       &
                  fileunits(I)
          enddo
          write (6,'(/)',advance='no')
       Else
          write (6,'(a)') 'readunits: File names not present'
#ifdef MPI
          call MPI_Abort (MPI_Comm_world, 1, MPIerror)
#else
          stop
#endif
       EndIf
    endif

#ifdef MPI
    call MPI_Bcast (fileunits, 30*(ntypeunits+2), MPI_Character, 0,     &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (unitshift, ntypeunits+2, MPI_Double_Precision, 0,   &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (ephIndic, ntypeunits+2, MPI_Integer, 0,             &
                    MPI_Comm_world, MPIerror)
#endif

!   Read unit's dimensions (number of orbitals).
    do I = 1,ntypeunits+2
       if (IOnode) then
          call io_assign (iu)
          open (iu, file=paste(directory, paste(fileunits(I),'.DAT')),  &
                status='old')
          write (6,'(a,a)') 'readunits: Reading file = ',               &
               paste(fileunits(I),'.DAT')
          read (iu,*) slabel, nuo, nspinu, maxnh, efu,                  &
                      tempu, nscu(1), nscu(2), no
          unitdimensions(I) = nuo
          call io_close (iu)
       endif
    enddo
#ifdef MPI
    call MPI_Bcast (unitdimensions, ntypeunits+2, MPI_Integer, 0,       &
                    MPI_Comm_world, MPIerror)
#endif

!   Read electron-phonon interaction data.
    call EPHread (ntypeunits, fileunits, unitdimensions, ephIndic)

!   Allocate and initialize units hamiltonian and overlap matrices.
    allocate (S1unit(unitdimensions(ntypeunits),                        &
                     unitdimensions(ntypeunits)))
    allocate (H1unit(unitdimensions(ntypeunits),                        &
                     unitdimensions(ntypeunits),nspin))
    S1unit = 0.d0
    H1unit = 0.d0
    allocate (Sunits(ntypeunits+2))
    allocate (Hunits(ntypeunits+2))
    do I = 1,ntypeunits+2
       allocate (Sunits(I)%S(unitdimensions(I),unitdimensions(I)))
       allocate (Hunits(I)%H(unitdimensions(I),unitdimensions(I),nspin))
       Sunits(I)%S = 0.d0
       Hunits(I)%H = 0.d0
    enddo

!   Read units on node 0.
    if (IOnode) then

       If (tightbinding) Then

!         Creates overlap and hamiltonians for tight-binding calculation.
          call TBreadunits (nspin, ntypeunits, unitdimensions, ephIndic)

       Else ! standart calculation (not tight-binding)

!         Allocate auxiliary matrices.
          allocate (H0aux(maxval(unitdimensions),                       &
                          maxval(unitdimensions),nspin))
          allocate (H1aux(maxval(unitdimensions),                       &
                          maxval(unitdimensions),nspin))
          allocate (S0aux(maxval(unitdimensions),maxval(unitdimensions)))
          allocate (S1aux(maxval(unitdimensions),maxval(unitdimensions)))

!         Read hamiltonian and overlap matrices for each unit.
          do I = 1,ntypeunits+2

             write (6,'(a,i3,a,i5,a,a,a,a)') 'readunits: Unit ', I,     &
                  ' - ', unitdimensions(I), '  ', trim(fileunits(I)),   &
                  '  ', trim(paste(directory,fileunits(I)))

!            Initialize auxiliary matrices.
             H0aux = 0.d0
             H1aux = 0.d0
             S0aux = 0.d0
             S1aux = 0.d0

!            Read hamiltonian and overlap matrices.
             call hsunits (nspin, nspin, unitdimensions(I), nsc, 1,     &
                      0, 0, .true., 1, temp,                            &
                      H0aux(1:unitdimensions(I),1:unitdimensions(I),:), &
                      H1aux(1:unitdimensions(I),1:unitdimensions(I),:), &
                      S0aux(1:unitdimensions(I),1:unitdimensions(I)),   &
                      S1aux(1:unitdimensions(I),1:unitdimensions(I)),   &
                      paste(directory,fileunits(I)))

!            Assign hamiltonian and overlap matrices from auxiliaries.
             Hunits(I)%H = H0aux(1:unitdimensions(I),                   &
                                 1:unitdimensions(I),:)
             Sunits(I)%S = S0aux(1:unitdimensions(I),1:unitdimensions(I))

             If (I == ntypeunits) Then
                H1unit = H1aux(1:unitdimensions(I),1:unitdimensions(I),:)
                S1unit = S1aux(1:unitdimensions(I),1:unitdimensions(I))
             EndIf

          enddo

!         Free memory.
          deallocate (H0aux, H1aux, S0aux, S1aux)

       EndIF ! If (tightbinding)

    endif ! if (IOnode)

!   Broadcast read unit's hamiltonians and overlap matrices.
#ifdef MPI
    do I = 1,ntypeunits+2
       call MPI_Bcast (Hunits(I)%H,                                     &
                       unitdimensions(I)*unitdimensions(I)*nspin,       &
                       MPI_Double_Precision, 0, MPI_Comm_world, MPIerror)
       call MPI_Bcast (Sunits(I)%S,                                     &
                       unitdimensions(I)*unitdimensions(I),             &
                       MPI_Double_Precision, 0, MPI_Comm_world, MPIerror)

       If (I == ntypeunits) Then
          call MPI_Bcast (H1unit, unitdimensions(I)*unitdimensions(I)   &
                          *nspin, MPI_Double_Precision, 0,              &
                          MPI_Comm_world, MPIerror)
          call MPI_Bcast (S1unit, unitdimensions(I)*unitdimensions(I),  &
                          MPI_Double_Precision, 0,                      &
                          MPI_Comm_world, MPIerror)
       EndIf

    enddo
#endif


  end subroutine readunits


!  *******************************************************************  !
!                              TBreadunits                              !
!  *******************************************************************  !
!  Description: Creates the unit's overlap and hamiltonian matrices     !
!  for tight-binding calculation.                                       !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2013.                              !
!  Instituto de Fisica Teorica                                          !
!  Universidade Estadual de Sao Paulo                                   !
!  e-mail: reilya@ift.unesp.br                                          !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode               : True if it is the I/O node            !
!  integer ephIdx(ntypeunits+2)        : Unit index (those with e-ph)   !
!  integer idxF(neph)                  : First dynamic atom orbital     !
!  integer idxL(neph)                  : Last dynamic atom orbital      !
!  real*8 TBenerg0              : Tight-binding site energy             !
!  real*8 TBenerg1              : Tight-binding deffect site energy     !
!  real*8 TBcoupl0              : Tight-binding site couplings          !
!  real*8 TBcoupl1              : Tight-binding deffect coupling        !
!  ****************************** INPUT ******************************  !
!  integer nspin                : Number of spin components             !
!  integer ntypeunits           : Number of unit types                  !
!  integer unitdimensions(ntypeunits+2)  : Units number of orbitals     !
!  integer ephIndic(ntypeunits+2)        : E-ph interaction indicator   !
!  ***************************** OUTPUT ******************************  !
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
!  *******************************************************************  !
  subroutine TBreadunits (nspin, ntypeunits, unitdimensions, ephIndic)

!
!   Modules
!
    use idsrdr_ephcoupl, only: ephIdx, idxF, idxL
    use idsrdr_options,  only: TBenerg0, TBenerg1, TBcoupl0, TBcoupl1
    use fdf

#ifdef MPI
    include "mpif.h"
#endif

!   Input variables.
    integer, intent(in) :: nspin, ntypeunits
    integer, dimension (ntypeunits+2), intent(in) :: unitdimensions
    integer, dimension (ntypeunits+2), intent(in) :: ephIndic

!   Local variables.
    integer :: I, w, r, s

!   Build hamiltonian and overlap matrices for each unit.
    w = 1
    I = ntypeunits+1
    Sunits(I)%S = 0.d0
    Hunits(I)%H = 0.d0
    do r = 1,unitdimensions(I)-1
       Hunits(I)%H(r,r,1) = TBenerg0
       Hunits(I)%H(r,r+1,1) = TBcoupl0
       Hunits(I)%H(r+1,r,1) = TBcoupl0
       Sunits(I)%S(r,r) = 1.d0
    enddo
    Hunits(I)%H(r,r,1) = TBenerg0
    Sunits(I)%S(r,r) = 1.d0
    If (ephIndic(I) == 1) Then
       r = idxF(w) - 1
       Hunits(I)%H(r,r+1,1) = TBcoupl1
       Hunits(I)%H(r+1,r,1) = TBcoupl1
       do r = idxF(w),idxL(w)
          Hunits(I)%H(r,r,1) = TBenerg1
          Hunits(I)%H(r,r+1,1) = TBcoupl1
          Hunits(I)%H(r+1,r,1) = TBcoupl1
       enddo
       w = w + 1
    EndIf

    do I = 1,ntypeunits
       Sunits(I)%S = 0.d0
       Hunits(I)%H = 0.d0
       do r = 1,unitdimensions(I)-1
          Hunits(I)%H(r,r,1) = TBenerg0
          Hunits(I)%H(r,r+1,1) = TBcoupl0
          Hunits(I)%H(r+1,r,1) = TBcoupl0
          Sunits(I)%S(r,r) = 1.d0
       enddo
       Hunits(I)%H(r,r,1) = TBenerg0
       Sunits(I)%S(r,r) = 1.d0
       If (ephIndic(I) == 1) Then
          r = idxF(w) - 1
          Hunits(I)%H(r,r+1,1) = TBcoupl1
          Hunits(I)%H(r+1,r,1) = TBcoupl1
          do r = idxF(w),idxL(w)
             Hunits(I)%H(r,r,1) = TBenerg1
             Hunits(I)%H(r,r+1,1) = TBcoupl1
             Hunits(I)%H(r+1,r,1) = TBcoupl1
          enddo
          w = w + 1
       EndIf
    enddo
    H1unit = 0.d0
    S1unit = 0.d0
    H1unit(unitdimensions(I),1,1) = TBcoupl0
    
    I = ntypeunits+2
    Sunits(I)%S = 0.d0
    Hunits(I)%H = 0.d0
    do r = 1,unitdimensions(I)-1
       Hunits(I)%H(r,r,1) = TBenerg0
       Hunits(I)%H(r,r+1,1) = TBcoupl0
       Hunits(I)%H(r+1,r,1) = TBcoupl0
       Sunits(I)%S(r,r) = 1.d0
    enddo
    Hunits(I)%H(r,r,1) = TBenerg0
    Sunits(I)%S(r,r) = 1.d0
    If (ephIndic(I) == 1) Then
       r = idxF(w) - 1
       Hunits(I)%H(r,r+1,1) = TBcoupl1
       Hunits(I)%H(r+1,r,1) = TBcoupl1
       do r = idxF(w),idxL(w)
          Hunits(I)%H(r,r,1) = TBenerg1
          Hunits(I)%H(r,r+1,1) = TBcoupl1
          Hunits(I)%H(r+1,r,1) = TBcoupl1
       enddo
    EndIf

    do s = 2,nspin ! over other spin components

       w = 1
       I = ntypeunits+1
       do r = 1,unitdimensions(I)-1
          Hunits(I)%H(r,r,s) = TBenerg0
          Hunits(I)%H(r,r+1,s) = TBcoupl0
          Hunits(I)%H(r+1,r,s) = TBcoupl0
       enddo
       Hunits(I)%H(r,r,s) = TBenerg0
       If (ephIndic(I) == 1) Then
          r = idxF(w) - 1
          Hunits(I)%H(r,r+1,s) = TBcoupl1
          Hunits(I)%H(r+1,r,s) = TBcoupl1
          do r = idxF(w),idxL(w)
             Hunits(I)%H(r,r,s) = TBenerg1
             Hunits(I)%H(r,r+1,s) = TBcoupl1
             Hunits(I)%H(r+1,r,s) = TBcoupl1
          enddo
          w = w + 1
       EndIf

       do I = 1,ntypeunits
          do r = 1,unitdimensions(I)-1
             Hunits(I)%H(r,r,s) = TBenerg0
             Hunits(I)%H(r,r+1,s) = TBcoupl0
             Hunits(I)%H(r+1,r,s) = TBcoupl0
          enddo
          Hunits(I)%H(r,r,s) = TBenerg0
          If (ephIndic(I) == 1) Then
             r = idxF(w) - 1
             Hunits(I)%H(r,r+1,s) = TBcoupl1
             Hunits(I)%H(r+1,r,s) = TBcoupl1
             do r = idxF(w),idxL(w)
                Hunits(I)%H(r,r,s) = TBenerg1
                Hunits(I)%H(r,r+1,s) = TBcoupl1
                Hunits(I)%H(r+1,r,s) = TBcoupl1
             enddo
             w = w + 1
          EndIf
       enddo
       H1unit(unitdimensions(I),1,s) = TBcoupl0
    
       I = ntypeunits+2
       do r = 1,unitdimensions(I)-1
          Hunits(I)%H(r,r,s) = TBenerg0
          Hunits(I)%H(r,r+1,s) = TBcoupl0
          Hunits(I)%H(r+1,r,s) = TBcoupl0
       enddo
       Hunits(I)%H(r,r,s) = TBenerg0
       If (ephIndic(I) == 1) Then
          r = idxF(w) - 1
          Hunits(I)%H(r,r+1,s) = TBcoupl1
          Hunits(I)%H(r+1,r,s) = TBcoupl1
          do r = idxF(w),idxL(w)
             Hunits(I)%H(r,r,s) = TBenerg1
             Hunits(I)%H(r,r+1,s) = TBcoupl1
             Hunits(I)%H(r+1,r,s) = TBcoupl1
          enddo
       EndIf

    enddo ! do s = 2,nspin

!   Write matrices on screen.
    do s = 1,nspin

       do I = 1,ntypeunits+2
          write (6,'(a,i5)') 'Hunit', I
          do r = 1,unitdimensions(I)
             do w = 1,unitdimensions(I)-1
                write (6,'(e17.8e3)',advance='no') Hunits(I)%H(r,w,s)
             enddo
             write (6,'(e17.8e3)') Hunits(I)%H(r,w,s)
          enddo
       enddo

       do I = 1,ntypeunits+2
          write (6,'(a,i5)') 'Sunit', I
          do r = 1,unitdimensions(I)
             do w = 1,unitdimensions(I)-1
                write (6,'(e17.8e3)',advance='no') Sunits(I)%S(r,w)
             enddo
             write (6,'(e17.8e3)') Sunits(I)%S(r,w)
          enddo
       enddo
       
       I = ntypeunits
       write (6,'(a)') 'H1unit'
       do r = 1,unitdimensions(I)
          do w = 1,unitdimensions(I)-1
             write (6,'(e17.8e3)',advance='no') H1unit(r,w,s)
          enddo
          write (6,'(e17.8e3)') H1unit(r,w,s)
       enddo
       write (6,'(a)') 'S1unit'
       do r = 1,unitdimensions(I)
          do w = 1,unitdimensions(I)-1
             write (6,'(e17.8e3)',advance='no') S1unit(r,w)
          enddo
          write (6,'(e17.8e3)') S1unit(r,w)
       enddo

    enddo
    write (6,'(a)') ""


  end subroutine TBreadunits


!  *******************************************************************  !
!                              buildunits                               !
!  *******************************************************************  !
!  Description: calculate the segments lengths and set the 'unit_type'  !
!  index vector.                                                        !
!                                                                       !
!  Use FDF (Flexible Data Format) package of J.M.Soler and A.Garcia.    !
!                                                                       !
!  Written by Alexandre Reily Rocha, --- 2007.                          !
!  Instituto de Fisica Teorica                                          !
!  Universidade Estadual de Sao Paulo                                   !
!  e-mail: reilya@ift.unesp.br                                          !
!  ***************************** HISTORY *****************************  !
!  Original version:    --- 2007                                        !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode               : True if it is the I/O node            !
!  logical readunitstf          : Read 'UnitIndex' block?               !
!  ******************** INPUT/OUTPUT FROM MODULES ********************  !
!  integer nunits               : Total number of units                 !
!  ****************************** INPUT ******************************  !
!  integer NDeffects                    : Number of deffects blocks     !
!  integer ntypeunits                   : Number of unit types          !
!  real*8 unitlength(ntypeunits+2)      : Units size (in z direction)   !
!  real*8 unitweight(ntypeunits+2)      : Units weight                  !
!  ************************** INPUT/OUTPUT ***************************  !
!  real*8 dist(NDeffects+1)             : Deffects random distribution  !
!  ***************************** OUTPUT ******************************  !
!  integer unit_type(nunits+2)          : Units types                   !
!  integer nunitseph                    : Number of units with eph      !
!  integer eph_type(nunitseph)          : Units types with eph          !
!  *******************************************************************  !
  subroutine buildunits (NDeffects, ntypeunits, unitlength,             &
                         unitweight, dist)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_options,  only: readunitstf, nunits
    use fdf

#ifdef MPI
    include "mpif.h"
#endif

!   Input variables.
    integer, intent(in) :: NDeffects, ntypeunits
    real(8), dimension (ntypeunits+2), intent(in) :: unitlength,        &
                                                     unitweight
    real(8), dimension (NDeffects+1), intent(inout) :: dist

!   Local variables.
    integer :: I, J, nunits_aux, iu, index, ueph
    integer, allocatable, dimension (:) :: randomdeffects
    real(8) :: length
    real(8), allocatable, dimension(:) :: total_length ! Segment length
    integer, external :: irandomizedeffects, irandomize_index
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

!   Allocate arrays.
    allocate (randomdeffects(NDeffects))

    IF (IOnode) THEN
       J = irandomizedeffects (ntypeunits-1,                            &
            unitweight(1:ntypeunits-1), .true.)

       If (.not. readunitstf) Then
          if (nunits == 0) then
             do I = 1,NDeffects+1
                J = irandomizedeffects (ntypeunits-1,                   &
                     unitweight(1:ntypeunits-1), .false.)
                if (I <= NDeffects) randomdeffects(I) = J
!               We can set to a random number if
!               if we have more than one deffect.
                length = 0.d0
                do while (length < dist(I)-unitlength(J))
                   length = length + unitlength(ntypeunits)
                   nunits = nunits + 1
                enddo
                if (I <= NDeffects) nunits = nunits + 1
             enddo
          elseif (nunits /= 0) then
             do I = 1,NDeffects
                J = irandomizedeffects (ntypeunits-1,                   &
                     unitweight(1:ntypeunits-1), .false.)
                randomdeffects(I) = J
             enddo
          endif
       EndIf
    ENDIF

#ifdef MPI
    call MPI_Bcast (nunits, 1, MPI_Integer, 0, MPI_Comm_world, MPIerror)
#endif

!   Allocate arrays.
    allocate (unit_type(nunits+2))
    allocate (total_length(NDeffects+1))

!   Build units types and calcutates the 'total_length'.
    IF (IOnode) THEN
       If (readunitstf) Then
          if (fdf_block ('UnitIndex', iu)) then
             do I = 1,nunits+2
                read (iu,*) unit_type(I)
             enddo
          else
#ifdef MPI
             call MPI_Abort (MPI_Comm_world, 1, MPIerror)
#else
             stop "buildunits: ERROR: No UnitIndex block defined!"
#endif
          endif
       ElseIf (nunits == 0) Then
          unit_type(1) = -1
          unit_type(2) = -2
          nunits_aux = 0

          do I = 1,NDeffects+1
             if (I <= NDeffects) J = randomdeffects(I)
             length = 0.0d0
             do while (length < dist(I)-unitlength(J))
                length = length + unitlength(ntypeunits)
                nunits_aux = nunits_aux + 1
                unit_type(nunits_aux+1) = ntypeunits
             enddo
             if (I <= NDeffects) then
                nunits_aux = nunits_aux + 1
                unit_type(nunits_aux+1) = J
                total_length(I) = length + unitlength(J)
             else
                total_length(I) = length
             endif
          enddo

       ElseIf (nunits /= 0) Then
          unit_type = ntypeunits
          unit_type(1) = -1
          unit_type(nunits+2) = -2

          do I = 1,NDeffects
             index = 1
             do while (unit_type(index) /= ntypeunits)
                index = irandomize_index (nunits, .false.)
             enddo
             write (6,'(/,a,i3,i3)') 'buildunits: index ', index,       &
                  randomdeffects(I)
             unit_type(index) = randomdeffects(I)
          enddo

          total_length = 0.0d0
          index = 2
          do I = 1,NDeffects+1
             write (6,'(a,i3,i3)') 'buildunits: index teste', index,    &
                  unit_type(index)
             do while (unit_type(index) == ntypeunits)
                total_length(I) = total_length(I)                       &
                     + unitlength(ntypeunits)
                index = index + 1
             enddo
             write (6,'(a,i3,i3)') '                       ', index,    &
                  unit_type(index)
             if (unit_type(index) /= -2) then
                total_length(I) = total_length(I)                       &
                     + unitlength(unit_type(index))
             endif
             index = index + 1
          enddo
       EndIf ! If (readunitstf)
    ENDIF ! IF (IOnode)

!   Read segments distribution and lengths.
    IF (readunitstf) THEN
       If (IOnode) Then
          if (fdf_block ('SegmentLengths', iu)) then
             do I = 1,NDeffects+1
                read (iu,*) dist(I), total_length(I)
             enddo
          else
#ifdef MPI
             call MPI_Abort (MPI_Comm_world, 1, MPIerror)
#else
             stop "buildunits: ERROR: Segment lengths not defined!"
#endif
          endif
       EndIf

#ifdef MPI
       call MPI_Bcast (dist, NDeffects+1, MPI_Double_Precision, 0,      &
                       MPI_Comm_world, MPIerror)
#endif
    ENDIF

#ifdef MPI
    call MPI_Bcast (unit_type, nunits+2, MPI_Integer, 0,                &
                    MPI_Comm_world, MPIerror)
#endif

!   Calculate the number of units with eph 'nunitseph'.
    nunitseph = 0
    if (ephIndic(ntypeunits+1)) eph_type(1) = unit_type(1)
    do I = 2,nunits+1
       if (ephIndic(unit_type(I))) nunitseph = nunitseph + 1
    enddo
    if (ephIndic(ntypeunits+2)) nunitseph = nunitseph + 1

!   Asign the array 'eph_type' indicating the units type with eph.
    allocate (eph_type(nunitseph))
    ueph = 1
    if (ephIndic(ntypeunits+1)) then
       eph_type(ueph) = ntypeunits + 1
       ueph = ueph + 1
    endif
    do I = 2,nunits+1
       if (ephIndic(unit_type(I))) then
          eph_type(ueph) = unit_type(I)
          ueph = ueph + 1
       endif
    enddo
    if (ephIndic(ntypeunits+2)) then
       eph_type(ueph) = ntypeunits+2
    endif

    if (IOnode) then
       write (6,'(a,i6)') 'buildunits: Number of segments', NDeffects+1
       do I = 1,NDeffects+1
          write (6,'(a,d16.8,a,d16.8)')                                 &
               'buildunits: Desired/Real segment length', dist(I),      &
               ' /', total_length(I)
       enddo
       write (6,'(/,a,i3)') 'buildunits: Total number of units = ',     &
            nunits
       write (6,'(a,i3)') 'buildunits: Total number of units with ' //  &
            'eph = ', nunitseph
       ueph = 1
       write (6,'(a,i3)', advance="no")                              &
            'buildunits: Unit type (eph) = ', unit_type(1)
       if (ephIndic(ntypeunits+1)) then
          write (6,'(a,i3,a)') ' (', eph_type(ueph), '  )'
          ueph = ueph + 1
       else
          write (6,'(a)') ' (  -  )'
       endif
       do I = 2,nunits+1
          write (6,'(a,i3)', advance="no")                              &
               'buildunits: Unit type (eph) = ', unit_type(I)
          if (ephIndic(unit_type(I))) then
             write (6,'(a,i3,a)') ' (', eph_type(ueph), '  )'
             ueph = ueph + 1
          else
             write (6,'(a)') ' (  -  )'
          endif
       enddo
       write (6,'(a,i3)', advance="no")                              &
            'buildunits: Unit type (eph) = ', unit_type(nunits+2)
       if (ephIndic(ntypeunits+2)) then
          write (6,'(a,i3,a)') ' (', eph_type(ueph), '  )'
          ueph = ueph + 1
       else
          write (6,'(a)') ' (  -  )'
       endif
    endif

!   Free memory.
    deallocate (randomdeffects)
    deallocate (total_length)


  end subroutine buildunits


!  *******************************************************************  !
!                               freeunits                               !
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
  subroutine freeunits

!
!   Modules
!
    use idsrdr_options,  only: ntypeunits

!   Local variables.
    integer :: I

!   Free memory.
    deallocate (theta)
    deallocate (unitshift)
    deallocate (S1unit)
    deallocate (H1unit)
!   First deallocates pointed matrices.
    do I = 1,ntypeunits
       deallocate (Sunits(I)%S)
       deallocate (Hunits(I)%H)
    enddo
    deallocate (Sunits)
    deallocate (Hunits)
    deallocate (unit_type)
    deallocate (unitdimensions)
    deallocate (ephIndic)
    deallocate (eph_type)


  end subroutine freeunits


!  *******************************************************************  !


END MODULE idsrdr_units
