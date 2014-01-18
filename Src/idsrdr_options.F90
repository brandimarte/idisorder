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

  logical :: calcdos      ! Calculate total DOS?
  logical :: readunitstf  ! Read 'UnitIndex' block?
  logical :: tightbinding ! Tight-binding calculation?
  logical :: writeondisk  ! Write Green's functions on disk?

  integer :: NDefects     ! Number of defects blocks
  integer :: NTenerg      ! Number of transmission energy points
  integer :: nspin        ! Number of spin components
  integer :: ntypeunits   ! Number of unit types
  integer :: nunits       ! Total number of units
  integer :: NIVP         ! Number of bias potential points
  integer :: symmetry     ! 
  integer :: norbitals    ! Number of orbitals
  integer :: numberrings  ! 
  integer :: nAsymmPts    ! Number of energy grid points
                          ! for asymmetric term integral

  integer, parameter :: label_length = 60 ! Length of system label

  integer, allocatable, dimension (:) :: atoms_per_ring ! 

  real(8) :: avgdist      ! Average defect distance
  real(8) :: TEnergI      ! Initial transmission energy
  real(8) :: TEnergF      ! Final transmission energy
  real(8) :: temp         ! Electronic temperature
  real(8) :: VInitial     ! Initial value of the bias potential
  real(8) :: VFinal       ! Final value of the bias potential
  real(8) :: dV           ! Bias potential step
  real(8) :: TBeFermi     ! Tight-binding Fermi energy
  real(8) :: TBenerg      ! Tight-binding site energy
  real(8) :: TBenergS     ! Tight-binding simple defect site energy
  real(8) :: TBenergSDB   ! Tight-binding simple defect site energy (DB)
  real(8) :: TBenergDB    ! Tight-binding dangling bond site energy
  real(8) :: TBcoupl      ! Tight-binding site couplings
  real(8) :: TBcouplS     ! Tight-binding simple defect coupling
  real(8) :: TBcouplSDB   ! Tight-binding simple defect coupling (DB)
  real(8) :: TBcouplDB    ! Tight-binding dangling bond coupling

  character(len=60) :: directory ! Working directory
  character(len=label_length), save :: slabel ! System Label
                                              ! (to name output files)
  character(len=14) :: integraltype ! Integration method


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
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer Nodes               : Total number of nodes (MPI_Comm_size)  !
!  ***************************** OUTPUT ******************************  !
!  logical calcdos           : Calculate total DOS?                     !
!  logical readunitstf       : Read 'UnitIndex' block?                  !
!  logical tightbinding      : Tight-binding calculation?               !
!  logical writeondisk       : Write Green's functions on disk?         !
!  integer NDefects          : Number of defects blocks                 !
!  integer NTenerg           : Number of transmission energy points     !
!  integer nspin             : Number of spin components                !
!  integer ntypeunits        : Number of unit types                     !
!  integer nunits            : Total number of units                    !
!  integer NIVP              : Number of bias potential points          !
!  integer symmetry          :                                          !
!  integer norbitals         : Number of orbitals                       !
!  integer numberrings       :                                          !
!  integer nAsymmPts         : Number of energy grid points             !
!                              for asymmetric term integral             !
!  integer atoms_per_ring(numberrings) :                                !
!  integer label_length      : Length of system label                   !
!  real*8 avgdist            : Average defect distance                  !
!  real*8 TEnergI            : Initial transmission energy              !
!  real*8 TEnergF            : Final transmission energy                !
!  real*8 temp               : Electronic temperature                   !
!  real*8 VInitial           : Initial value of the bias potential      !
!  real*8 VFinal             : Final value of the bias potential        !
!  real*8 dV                 : Bias potential step                      !
!  real*8 TBeFermi           : Tight-binding Fermi energy               !
!  real*8 TBenerg            : Tight-binding site energy                !
!  real*8 TBenergS           : Tight-binding simple defect site energy  !
!  real*8 TBenergSDB         : Tight-binding simple defect site energy  !
!                              (to DB)                                  !
!  real*8 TBenergDB          : Tight-binding dangling bond site energy  !
!  real*8 TBcoupl            : Tight-binding site couplings             !
!  real*8 TBcouplS           : Tight-binding simple defect coupling     !
!  real*8 TBcouplSDB         : Tight-binding simple defect coupling     !
!                              (to DB)                                  !
!  real*8 TBcouplDB          : Tight-binding dangling bond coupling     !
!  character(60) directory   : Working directory                        !
!  character(label_length) slabel : System Label (for output files)     !
!  character(14) integraltype : Integration method                      !
!  *******************************************************************  !
  subroutine readopt

!
!   Modules
!
    use parallel,        only: IOnode, Nodes

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    integer :: iu, j
    logical :: spinpol, isblock
    character :: slabel_default*60
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

    if (IOnode) then
       write (6,'(/,28("*"),a,28("*"))')                              &
            ' Simulation parameters '

!      Defile System Label (short name to label files).
       slabel = ""
       slabel_default = 'i-disorder'
       slabel = fdf_string ('SystemLabel', slabel_default)
       write (6,2) 'readopt: System label                         ' //  &
            '         =  ', slabel

!      Number of defects blocks.
       NDefects = fdf_integer ('NumberDefects', 1)
       write (6,4)                                                      &
            'readopt: Number of defects blocks                      =', &
            NDefects

!      Average defect distance.
       avgdist = fdf_physical ('AvgDefectDist', 50.0d0, 'Ang')
       write (6,6)                                                      &
            'readopt: Average defect distance                       =', &
            avgdist, ' Ang'

!      Number of spin components.
       spinpol = fdf_boolean ('SpinPolarized', .false.)
       if (spinpol) then
          nspin = 2
       else
          nspin = 1
       endif
       write (6,4)                                                      &
            'readopt: Number of spin components                     =', &
            nspin

!      Number of unit types.
       ntypeunits = fdf_integer ('NumberUnitTypes', 2)
       write (6,4)                                                      &
            'readopt: Number of unit types                          =', &
            ntypeunits

!      Read 'UnitIndex' block?
       readunitstf = fdf_boolean ('ReadUnits', .false.)
       write (6,1)                                                      &
            'readopt: Read UnitIndex and SegmentLengths blocks?     =', &
            readunitstf

!      Total number of units.
       nunits = 0
       nunits = fdf_integer ('NumberUnits', 0)
       if (readunitstf) then
          if (nunits == 0) then
#ifdef MPI
             call MPI_Abort (MPI_Comm_world, 1, MPIerror)
#else
             stop 'readopt: ERROR: Number of units is zero!'
#endif
          endif
          write (6,4)                                                   &
               'readopt: Total number of units                    ' //  &
               '     =', nunits
       endif

!      
       symmetry = fdf_integer ('Symmetry', 1)
       write (6,4)                                                      &
            'readopt: Simmetry                                      =', &
            symmetry

!      Number of orbitals.
       norbitals = fdf_integer ('NumberOfOrbitals', 1)
       write (6,4)                                                      &
            'readopt: Number of orbitals                            =', &
            norbitals

!      Electronic temperature.
       temp = fdf_physical ('ElectronicTemperature', 300.0d0, 'Ry')
       write (6,6)                                                      &
            'readopt: Electronic temperature                        =', &
            temp, ' Ry'

!      Number of bias potential points.
       NIVP = fdf_integer ('NIVPoints', 10)
       write(6,4)                                                       &
            'readopt: Number of bias potential points               =', &
            NIVP

!      Initial value of the bias potential.
       VInitial = fdf_physical ('VInitial', 0.0d0, 'Ry')
       write(6,6)                                                       &
            'readopt: Initial value of the bias potential           =', &
            VInitial, ' Ry'

!      Final value of the bias potential.
       VFinal = fdf_physical ('VFinal', 0.1d0, 'Ry')
       write(6,6)                                                       &
            'readopt: Final value of the bias potential             =', &
            VFinal, ' Ry'

!      Number of transmission energy points.
       NTenerg = fdf_integer ('NumberTransmPoints', 100)
       write (6,4)                                                      &
            'readopt: Number of transmission energy points          =', &
            NTenerg

       if (NTenerg > 0) then

!         Initial transmission energy.
          TEnergI = fdf_physical ('TransmInitial', -1.d0, 'Ry')
          write (6,6)                                                   &
               'readopt: Initial transmission energy              ' //  &
               '     =', TEnergI, ' Ry'

!         Final transmission energy.
          TEnergF = fdf_physical ('TransmFinal', 1.d0, 'Ry')
          write (6,6)                                                   &
               'readopt: Final transmission energy                ' //  &
               '     =', TEnergF, ' Ry'
       endif

!      Calculate total DOS?
       calcdos = fdf_boolean ('CalculateDos', .false.)
       write (6,1)                                                      &
            'readopt: Calculate total DOS?                          =', &
            calcdos

!      Number of rings (?).
       numberrings = fdf_integer('NumberOfrings', 1)
       write (6,4)                                                      &
            'readopt: Number of rings                               =', &
            numberrings

!      Atoms per ring (?).
       allocate (atoms_per_ring(numberrings))
       atoms_per_ring = 0
       isblock = fdf_block ('AtomsPerRing', iu)
       if (isblock) then
          do j = 1,numberrings
             read (iu,*) atoms_per_ring(j)
          enddo
          write (6,'(a,i7,a,i3)') 'readopt: Atoms per ring        ' //  &
               '                        =', 1, ' - ', atoms_per_ring(1)
          do j = 2,numberrings
             write (6,'(a,i8,a,i3)') '                            ' //  &
                  '                           ', j, ' - ',              &
                  atoms_per_ring(j)
          enddo
       else
          write (6,2) 'readopt: Atoms per ring                     ' // &
               '           = You must specify the number per ring.'
       endif

!      Working directory.
       directory = fdf_string ('Directory', ' ')
       write (6,2) 'readopt: Working directory                    ' //  &
            '         =  ', directory	

!      Integration method.
       integraltype = ""
       integraltype = fdf_string ('TypeOfIntegral', 'Sympson')
       write (6,2) 'readopt: Integration method                   ' //  &
            '         =  ', integraltype

!      Number of energy grid points for asymmetric term integral.
       nAsymmPts = fdf_integer('AsymmGridPts', 1000)
       if (nAsymmPts == 0) nAsymmPts = 1
       do while (MOD(nAsymmPts,Nodes) /= 0)
          nAsymmPts = nAsymmPts + 1
       enddo
       write (6,4)                                                      &
            'readopt: Number of points at asymmetric term integral  =', &
            nAsymmPts

!      Tight-binding calculation?
       tightbinding = fdf_boolean ('TightBinding', .false.)
       write (6,1)                                                      &
            'readopt: Tight-binding calculation?                    =', &
            tightbinding

       if (tightbinding) then

!         Tight-binding Fermi energy.
          TBeFermi = fdf_physical ('TB.eFermi', 0.d0, 'Ry')
          write (6,6)                                                   &
               'readopt: Tight-binding Fermi energy               ' //  &
               '     =', TBeFermi, ' Ry'

!         Tight-binding site energy.
          TBenerg = fdf_physical ('TB.energy', -0.5d0, 'Ry')
          write (6,6)                                                   &
               'readopt: Tight-binding site energy                ' //  &
               '     =', TBenerg, ' Ry'

!         Tight-binding site couplings.
          TBcoupl = fdf_physical ('TB.coupling', 0.7d0, 'Ry')
          write (6,6)                                                   &
               'readopt: Tight-binding site couplings             ' //  &
               '     =', TBcoupl, ' Ry'

!         Tight-binding simple defect site energy.
          TBenergS = fdf_physical ('TB.energy.S', -0.1d0, 'Ry')
          write (6,6)                                                   &
               'readopt: Tight-binding simple defect site energy  ' //  &
               '     =', TBenergS, ' Ry'

!         Tight-binding simple defect coupling.
          TBcouplS = fdf_physical ('TB.coupling.S', 0.175d0, 'Ry')
          write (6,6)                                                   &
               'readopt: Tight-binding simple defect coupling     ' //  &
               '     =', TBcouplS, ' Ry'

          if (ntypeunits > 2) then

!            Tight-binding simple defect site energy (to DB).
             TBenergSDB = fdf_physical ('TB.energy.SDB', -0.1d0, 'Ry')
             write (6,6)                                                &
                  'readopt: Tight-binding simple defect site '      //  &
                  'energy (DB)  =', TBenergSDB, ' Ry'

!            Tight-binding dangling bond site energy.
             TBenergDB = fdf_physical ('TB.energy.DB', -0.1d0, 'Ry')
             write (6,6)                                                &
                  'readopt: Tight-binding dangling bond site '      //  &
                  'energy       =', TBenergDB, ' Ry'

!            Tight-binding simple defect coupling (to DB).
             TBcouplSDB = fdf_physical ('TB.coupling.SDB', 0.175d0, 'Ry')
             write (6,6)                                                &
                  'readopt: Tight-binding simple defect coupling '  //  &
                  '(DB)     =', TBcouplSDB, ' Ry'

!            Tight-binding dangling bond coupling.
             TBcouplDB = fdf_physical ('TB.coupling.DB', 0.175d0, 'Ry')
             write (6,6)                                                &
                  'readopt: Tight-binding dangling bond coupling  ' //  &
                  '        =', TBcouplDB, ' Ry'

          endif

       endif

!      Write Green's functions on disk?
       writeondisk = fdf_boolean ('WriteOnDisk', .true.)
       write (6,1)                                                      &
            'readopt: Write Green s functions on disk?              =', &
            writeondisk

       write (6,'(2a)') 'readopt: ', repeat('*', 70)

    endif ! if (IOnode)

#ifdef MPI
    call MPI_Bcast (slabel, label_length, MPI_Character, 0,             &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (NDefects, 1, MPI_Integer, 0,                        &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (avgdist, 1, MPI_Double_Precision, 0,                &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (nspin, 1, MPI_Integer, 0, MPI_Comm_world, MPIerror)
    call MPI_Bcast (ntypeunits, 1, MPI_Integer, 0,                      &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (symmetry, 1, MPI_Integer, 0,                        &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (norbitals, 1, MPI_Integer, 0,                       &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (temp, 1, MPI_Double_Precision, 0,                   &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (NIVP, 1, MPI_Integer, 0, MPI_Comm_world, MPIerror)
    call MPI_Bcast (VInitial, 1, MPI_Double_Precision, 0,               &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (VFinal, 1, MPI_Double_Precision, 0,                 &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (NTenerg, 1, MPI_Integer, 0, MPI_Comm_world, MPIerror)
    call MPI_Bcast (TEnergI, 1, MPI_Double_Precision, 0,                &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (TEnergF, 1, MPI_Double_Precision, 0,                &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (calcdos, 1, MPI_Logical, 0, MPI_Comm_world, MPIerror)
    call MPI_Bcast (numberrings, 1, MPI_Integer, 0,                     &
                    MPI_Comm_world, MPIerror)
    if (.not. IOnode) allocate (atoms_per_ring(numberrings))
    call MPI_Bcast (atoms_per_ring, numberrings, MPI_Integer, 0,        &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (directory, 60, MPI_Character, 0,                    &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (integraltype, 14, MPI_Character, 0,                 &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (nAsymmPts, 1, MPI_Integer, 0,                       &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (readunitstf, 1, MPI_Logical, 0,                     &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (tightbinding, 1, MPI_Logical, 0,                    &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (TBeFermi, 1, MPI_Double_Precision, 0,               &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (TBenerg, 1, MPI_Double_Precision, 0,                &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (TBcoupl, 1, MPI_Double_Precision, 0,                &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (TBenergS, 1, MPI_Double_Precision, 0,               &
                    MPI_Comm_world, MPIerror)
    call MPI_Bcast (TBcouplS, 1, MPI_Double_Precision, 0,               &
                    MPI_Comm_world, MPIerror)
    if (ntypeunits > 2) then
       call MPI_Bcast (TBenergSDB, 1, MPI_Double_Precision, 0,          &
                       MPI_Comm_world, MPIerror)
       call MPI_Bcast (TBenergDB, 1, MPI_Double_Precision, 0,           &
                       MPI_Comm_world, MPIerror)
       call MPI_Bcast (TBcouplSDB, 1, MPI_Double_Precision, 0,          &
                       MPI_Comm_world, MPIerror)
       call MPI_Bcast (TBcouplDB, 1, MPI_Double_Precision, 0,           &
                       MPI_Comm_world, MPIerror)
    endif
    call MPI_Bcast (writeondisk, 1, MPI_Logical, 0,                     &
                    MPI_Comm_world, MPIerror)
!   It is not necessary to broadcast 'nunits' here.
#endif

!   Bias potential step.
    if(NIVP == 0) then
       dV = 0.d0
    else
       dV = 1.d0 * (VFinal - VInitial) / (1.d0 * NIVP)
    endif

1   format(a,6x,l1)
2   format(a,a)
4   format(a,i7)
6   format(a,f12.4,a)


  end subroutine readopt


!  *******************************************************************  !
!                                freeopt                                !
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
  subroutine freeopt


!   Free memory.
    deallocate (atoms_per_ring)


  end subroutine freeopt


!  *******************************************************************  !


END MODULE idsrdr_options
