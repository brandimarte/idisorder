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
!                         MODULE idsrdr_options                         !
!  *******************************************************************  !
!  Description: read input option.                                      !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !

MODULE idsrdr_options

!
!   Modules
!
  use parallel,        only: 
  use fdf

  implicit none
  
  PUBLIC ! default is public

  logical :: calcdos             ! Calculate total DOS?

  integer :: NDeffects           ! Number of deffects blocks
  integer :: NTenerg             ! Number of transmission energy points
  integer :: nspin               ! Number of spin components
  integer :: ntypeunits          ! Number of unit types
  integer :: symmetry            ! 
  integer :: norbitals           ! Number of orbitals
  integer :: numberrings         ! 
  integer, parameter :: label_length = 60 ! Length of system label

  integer, allocatable, dimension (:) :: atoms_per_ring ! 

  real*8 :: avgdist              ! Average deffect distance
  real*8 :: TEnergI              ! Initial transmission energy
  real*8 :: TEnergF              ! Final transmission energy
  real*8 :: temp                 ! Electronic temperature

  character(len=60) :: directory ! Working directory
  character(len=label_length), save :: slabel ! System Label
                                              ! (to name output files)


CONTAINS


!  *******************************************************************  !
!                                readopt                                !
!  *******************************************************************  !
!  Description: subroutine to read input variables.                     !
!                                                                       !
!  Use FDF (Flexible Data Format) package of J.M.Soler and A.Garcia.    !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ***************************** OUTPUT ******************************  !
!  logical calcdos              : Calculate total DOS?                  !
!  integer NDeffects            : Number of deffects blocks             !
!  integer NTenerg              : Number of transmission energy points  !
!  integer nspin                : Number of spin components             !
!  integer ntypeunits           : Number of unit types                  !
!  integer symmetry             :                                       !
!  integer norbitals            : Number of orbitals                    !
!  integer numberrings          :                                       !
!  integer atoms_per_ring(numberrings) :                                !
!  integer label_length         : Length of system label
!  real*8 avgdist               : Average deffect distance              !
!  real*8 TEnergI               : Initial transmission energy           !
!  real*8 TEnergF               : Final transmission energy             !
!  real*8 temp                  : Electronic temperature                !
!  character(60) directory      : Working directory                     !
!  character(label_length) slabel : System Label (for output files)     !
!  *******************************************************************  !
  subroutine readopt

!
!   Modules
!
    use parallel,        only: IOnode

    include "mpif.h"

!   Local variables.
    integer :: iu, j
    logical :: spinpol, isblock
    character :: slabel_default*60
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

    if (IOnode) then
       write(6,'(/,a,16("*"),a,17("*"))')                               &
            'I-DISORDER: ', ' I-DISORDER Simulation parameters '
       write(6,'(a)') 'I-DISORDER:'
       write(6,'(a)') 'I-DISORDER:     The following are some of ' //   &
            'the I-DISORDER parameters of the'
       write(6,'(a)') 'I-DISORDER:     simulation. A complete ' //      &
            'list of the parameters used,'
       write(6,'(a,a)') 'I-DISORDER:     including default values,' //  &
            ' can be found in file fdf.log.'
       write(6,'(a)')  'I-DISORDER:'

!      Defile System Label (short name to label files).
       slabel = ""
       slabel_default = 'i-disorder'
       slabel = fdf_string ('SystemLabel', slabel_default)
!!$       label_length = size (slabel)
       write(6,2) 'redata: System label                         ' //    &
            '         =  ', slabel

!      Number of deffects blocks.
       NDeffects = fdf_integer ('NumberDeffects', 1)
       write(6,4)                                                       &
            'redata: Number of deffects blocks                     =',  &
            NDeffects

!      Average deffect distance.
       avgdist = fdf_physical ('AvgDeffectDist', 50.0d0, 'Ang')
       write(6,6)                                                       &
            'redata: Average deffect distance                      =',  &
            avgdist, ' Ang'

!      Number of spin components.
       spinpol = fdf_boolean ('SpinPolarized', .false.)
       if (spinpol) then
          nspin = 2
       else
          nspin = 1
       endif
       write(6,4)                                                       &
            'redata: Number of spin components                     =',  &
            nspin

!      Number of unit types.
       ntypeunits = fdf_integer ('NumberUnitTypes', 2)
       write(6,4)                                                       &
            'redata: Number of unit types                          =',  &
            ntypeunits

!      
       symmetry = fdf_integer ('Symmetry', 1)
       write(6,4)                                                       &
            'redata: Simmetry                                      =',  &
            symmetry

!      Number of orbitals.
       norbitals = fdf_integer ('NumberOfOrbitals', 1)
       write(6,4)                                                       &
            'redata: Number of orbitals                            =',  &
            norbitals

!      Electronic temperature.
       temp = fdf_physical ('ElectronicTemperature', 300.0d0, 'Ry')
       write(6,6)                                                       &
            'redata: Electronic temperature                        =',  &
            temp, ' Ry'

!      Number of transmission energy points.
       NTenerg = fdf_integer ('NumberTransmPoints', 100)
       write(6,4)                                                       &
            'redata: Number of transmission energy points          =',  &
            NTenerg

!      Initial transmission energy.
       TEnergI = fdf_physical ('TransmInitial', -1.d0, 'Ry')
       write(6,6)                                                       &
            'redata: Initial transmission energy                   =',  &
            TEnergI, ' Ry'

!      Final transmission energy.
       TEnergF = fdf_physical ('TransmFinal', 1.d0, 'Ry')
       write(6,6)                                                       &
            'redata: Final transmission energy                     =',  &
            TEnergF, ' Ry'

!      Calculate total DOS?
       calcdos = fdf_boolean ('CalculateDos', .false.)
       write(6,1)                                                       &
            'redata: Calculate total DOS?                          =',  &
            calcdos

!      Number of rings (?).
       numberrings = fdf_integer('NumberOfrings',1)
       write(6,4)                                                       &
            'redata: Number of rings                               =',  &
            numberrings

!      Atoms per ring (?).
       allocate (atoms_per_ring(numberrings))
       atoms_per_ring = 0
       isblock = fdf_block ('AtomsPerRing', iu)
       if (isblock) then
          do j = 1,numberrings
             read (iu,*) atoms_per_ring(j)
          enddo
          write(6,'(a,i7,a,i3)') 'redata: Atoms per ring        ' //    &
               '                        =', 1, ' - ', atoms_per_ring(1)
          do j = 2,numberrings
             write(6,'(a,i8,a,i3)') '                              ' // &
                  '                        ', j, ' - ', atoms_per_ring(j)
          enddo
       else
          write(6,2) 'redata: Atoms per ring                       ' // &
               '         = You must specify the number per ring.'
       endif

!      Working directory.
       directory = fdf_string ('Directory',' ')
       write(6,2) 'redata: Working directory                    ' //    &
            '         =', directory	

       write(6,'(2a)') 'redata: ', repeat('*', 71)

    endif

#ifdef MPI
    call MPI_Bcast (slabel, label_length, MPI_Character, 0,             &
         MPI_Comm_World, MPIerror)
    call MPI_Bcast (NDeffects, 1, MPI_Integer, 0, MPI_Comm_World,       &
         MPIerror)
    call MPI_Bcast (avgdist, 1, MPI_Double_Precision, 0,                &
         MPI_Comm_World, MPIerror)
    call MPI_Bcast (nspin, 1, MPI_Integer, 0, MPI_Comm_World, MPIerror)
    call MPI_Bcast (ntypeunits, 1, MPI_Integer, 0, MPI_Comm_World,      &
         MPIerror)
    call MPI_Bcast (symmetry, 1, MPI_Integer, 0, MPI_Comm_World,        &
         MPIerror)
    call MPI_Bcast (norbitals, 1, MPI_Integer, 0, MPI_Comm_World,       &
         MPIerror)     
    call MPI_Bcast (temp, 1, MPI_Double_Precision, 0, MPI_Comm_World,   &
         MPIerror)
    call MPI_Bcast (NTenerg, 1, MPI_Integer, 0, MPI_Comm_World, MPIerror)
    call MPI_Bcast (TEnergI, 1, MPI_Double_Precision, 0,                &
         MPI_Comm_World, MPIerror)
    call MPI_Bcast (TEnergF, 1, MPI_Double_Precision, 0,                &
         MPI_Comm_World, MPIerror)
    call MPI_Bcast (directory, 50, MPI_Character, 0,                    &
         MPI_Comm_World, MPIerror)
    call MPI_Bcast (calcdos, 1, MPI_Logical, 0, MPI_Comm_World, MPIerror)
    call MPI_Bcast (numberrings, 1, MPI_Integer, 0, MPI_Comm_World,     &
         MPIerror)
    if (.not. IOnode) allocate (atoms_per_ring(numberrings))
    call MPI_Bcast (atoms_per_ring, numberrings, MPI_Integer, 0,        &
         MPI_Comm_World, MPIerror)
#endif

1   format(a,6x,l1)
2   format(a,a)
4   format(a,i7)
5   format(a,i7,a)
6   format(a,f12.4,a)
7   format(a,f14.6,a)
8   format(a,f20.12)
9   format(a,f13.5)


  end subroutine readopt


!  Free memory subroutine to be called at the end.

!  *******************************************************************  !


END MODULE idsrdr_options
