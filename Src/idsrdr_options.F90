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

  logical :: readunitstf  ! Read 'UnitIndex' block?
  logical :: tightbinding ! Tight-binding calculation?
  logical :: writeondisk  ! Write Green's functions on disk?
  logical :: phonEqui     ! Consider the phonons in equilibrium?

  integer :: NDefects     ! Number of defects blocks
  integer :: NTenerg      ! Number of transmission energy points
  integer :: nspin        ! Number of spin components
  integer :: ntypeunits   ! Number of unit types
  integer :: nunits       ! Total number of units
  integer :: NIVP         ! Number of bias potential points
  integer :: norbitals    ! Number of orbitals
  integer :: nAsymmPts    ! Number of energy grid points
                          ! for asymmetric term integral
  integer :: ProcsPerGPU  ! Number of processes running per GPU
  integer :: SigmaMethod  ! Method for calculating the lead self-energy
  integer, parameter :: label_length = 60 ! Length of system label

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
  real(8) :: phonDamp     ! Phenomenological damping parameter
                          ! (~ inverse of phonon's lifetime)
  real(8) :: imDelta      ! Small imaginary part
  real(8) :: deltaEn      ! energy step size (defined at 'idsrdr_engrid')


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
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer Nodes               : Total number of nodes (MPI_Comm_size)  !
!  integer MPI_Comm_MyWorld    : MPI communicator                       !
!  ***************************** OUTPUT ******************************  !
!  logical readunitstf       : Read 'UnitIndex' block?                  !
!  logical tightbinding      : Tight-binding calculation?               !
!  logical writeondisk       : Write Green's functions on disk?         !
!  logical phonEqui          : Consider the phonons in equilibrium?     !
!  integer NDefects          : Number of defects blocks                 !
!  integer NTenerg           : Number of transmission energy points     !
!  integer nspin             : Number of spin components                !
!  integer ntypeunits        : Number of unit types                     !
!  integer nunits            : Total number of units                    !
!  integer NIVP              : Number of bias potential points          !
!  integer nAsymmPts         : Number of energy grid points             !
!                              for asymmetric term integral             !
!  integer ProcsPerGPU       : Number of processes running per GPU      !
!  integer SigmaMethod       : Method for calculating self-energy       !
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
!  real*8 phonDamp           : Phenomenological damping parameter       !
!                              (~ inverse of phonon's lifetime)         !
!  real*8 imDelta            : Small imaginary part                     !
!  character(60) directory   : Working directory                        !
!  character(label_length) slabel : System Label (for output files)     !
!  *******************************************************************  !
  subroutine readopt

!
!   Modules
!
#ifdef MPI
    use parallel,        only: IOnode, Nodes, MPI_Comm_MyWorld
#else
    use parallel,        only: IOnode, Nodes
#endif

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    logical :: spinpol
    character :: slabel_default*60
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

    if (IOnode) then
       write (6,'(/,28("*"),a,28("*"))')                                &
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
       avgdist = 10.0d0
!!$       avgdist = fdf_physical ('AvgDefectDist', 50.0d0, 'Ang')
!!$       write (6,6)                                                      &
!!$            'readopt: Average defect distance                       =', &
!!$            avgdist, ' Ang'

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
       nunits = fdf_integer ('NumberUnits', 1)
       if (readunitstf) then
          if (nunits == 0) then
#ifdef MPI
             call MPI_Abort (MPI_Comm_MyWorld, 1, MPIerror)
#else
             stop 'readopt: ERROR: Number of units is zero!'
#endif
          endif
          write (6,4)                                                   &
               'readopt: Total number of units                    ' //  &
               '     =', nunits
       endif

!      Electronic temperature.
       temp = fdf_physical ('ElectronicTemperature', 300.0d0, 'K')
       write (6,6)                                                      &
            'readopt: Electronic temperature                        =', &
            temp, ' Ry'

!      Number of bias potential points.
       NIVP = fdf_integer ('NIVPoints', 10)
       write(6,4)                                                       &
            'readopt: Number of bias potential points               =', &
            NIVP

!      Initial value of the bias potential.
       VInitial = fdf_physical ('VInitial', 0.0d0, 'eV')
       write(6,6)                                                       &
            'readopt: Initial value of the bias potential           =', &
            VInitial, ' Ry'

!      Final value of the bias potential.
       VFinal = fdf_physical ('VFinal', 0.1d0, 'eV')
       write(6,6)                                                       &
            'readopt: Final value of the bias potential             =', &
            VFinal, ' Ry'

!      Number of transmission energy points.
       NTenerg = fdf_integer ('NumberTransmPoints', 100)
       if (NTenerg == 0) NTenerg = 1
       write (6,4)                                                      &
            'readopt: Number of transmission energy points          =', &
            NTenerg

       if (NTenerg > 0) then

!         Initial transmission energy.
          TEnergI = fdf_physical ('TransmInitial', -1.d0, 'eV')
          write (6,6)                                                   &
               'readopt: Initial transmission energy              ' //  &
               '     =', TEnergI, ' Ry'

!         Final transmission energy.
          TEnergF = fdf_physical ('TransmFinal', 1.d0, 'eV')
          write (6,6)                                                   &
               'readopt: Final transmission energy                ' //  &
               '     =', TEnergF, ' Ry'
       endif

!      Method for calculating the lead self-energy.
       SigmaMethod = fdf_integer ('SigmaMethod', 0)
       write(6,4)                                                       &
            'readopt: Lead self-energy method                       =', &
            SigmaMethod

       if (SigmaMethod == 1) then

!         Small imaginary part.
          imDelta  = fdf_double ('DeltaImag', 1.d-4)
          write(6,9)                                                    &
               'readopt: Small imaginary part                     ' //  &
               '     =', imDelta

       endif

!      Number of energy grid points for asymmetric term integral.
       nAsymmPts = fdf_integer('AsymmGridPts', 1000)
       write (6,4)                                                      &
            'readopt: Number of points at asymmetric term integral  =', &
            nAsymmPts

!      Consider the phonons in equilibrium?
       phonEqui = fdf_boolean ('PhononEquilibrium', .true.)
       write (6,1)                                                      &
            'readopt: Consider the phonons in equilibrium?     '    //  &
            '     =', phonEqui

       if (.not. phonEqui) then
!         Phenomenological damping parameter
!         (related to the inverse of phonon's lifetime).
          phonDamp = fdf_physical ('PhononDamping', 0.05d0, 'eV')
          write (6,6)                                                   &
               'readopt: Phenomenological phonon damping parameter' //  &
               '     =', phonDamp, ' Ry'
       endif

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

!      Number of processes running in each GPU (default = 1).
       ProcsPerGPU = fdf_integer ('GPU.ProcsPerGPU', 1)
       write (6,4)                                                &
            'readopt: Number of processes running per GPU   ' //  &
            '        =', ProcsPerGPU

!      Working directory.
       directory = fdf_string ('Directory', ' ')
       write (6,2) 'readopt: Working directory                    ' //  &
            '         =  ', directory	

       write (6,'(2a)') 'readopt: ', repeat('*', 70)

    endif ! if (IOnode)

#ifdef MPI
    call MPI_Bcast (slabel, label_length, MPI_Character, 0,             &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (NDefects, 1, MPI_Integer, 0,                        &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (avgdist, 1, MPI_Double_Precision, 0,                &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (nspin, 1, MPI_Integer, 0, MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (ntypeunits, 1, MPI_Integer, 0,                      &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (temp, 1, MPI_Double_Precision, 0,                   &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (NIVP, 1, MPI_Integer, 0, MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (VInitial, 1, MPI_Double_Precision, 0,               &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (VFinal, 1, MPI_Double_Precision, 0,                 &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (NTenerg, 1, MPI_Integer, 0,                         &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (TEnergI, 1, MPI_Double_Precision, 0,                &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (TEnergF, 1, MPI_Double_Precision, 0,                &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (SigmaMethod, 1, MPI_Integer, 0,                     &
                    MPI_Comm_MyWorld, MPIerror)
    if (SigmaMethod == 1) then
       call MPI_Bcast (imDelta, 1, MPI_Double_Precision, 0,             &
                       MPI_Comm_MyWorld, MPIerror)
    endif
    call MPI_Bcast (nAsymmPts, 1, MPI_Integer, 0,                       &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (readunitstf, 1, MPI_Logical, 0,                     &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (tightbinding, 1, MPI_Logical, 0,                    &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (TBeFermi, 1, MPI_Double_Precision, 0,               &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (TBenerg, 1, MPI_Double_Precision, 0,                &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (TBcoupl, 1, MPI_Double_Precision, 0,                &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (TBenergS, 1, MPI_Double_Precision, 0,               &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (TBcouplS, 1, MPI_Double_Precision, 0,               &
                    MPI_Comm_MyWorld, MPIerror)
    if (ntypeunits > 2) then
       call MPI_Bcast (TBenergSDB, 1, MPI_Double_Precision, 0,          &
                       MPI_Comm_MyWorld, MPIerror)
       call MPI_Bcast (TBenergDB, 1, MPI_Double_Precision, 0,           &
                       MPI_Comm_MyWorld, MPIerror)
       call MPI_Bcast (TBcouplSDB, 1, MPI_Double_Precision, 0,          &
                       MPI_Comm_MyWorld, MPIerror)
       call MPI_Bcast (TBcouplDB, 1, MPI_Double_Precision, 0,           &
                       MPI_Comm_MyWorld, MPIerror)
    endif
    call MPI_Bcast (phonEqui, 1, MPI_Logical, 0,                        &
                    MPI_Comm_MyWorld, MPIerror)
    if (.not. phonEqui) then
       call MPI_Bcast (phonDamp, 1, MPI_Double_Precision, 0,            &
                       MPI_Comm_MyWorld, MPIerror)
    endif
    call MPI_Bcast (writeondisk, 1, MPI_Logical, 0,                     &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (ProcsPerGPU, 1, MPI_Integer, 0,                     &
                    MPI_Comm_MyWorld, MPIerror)
    call MPI_Bcast (directory, 60, MPI_Character, 0,                    &
                    MPI_Comm_MyWorld, MPIerror)
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
6   format(a,f14.8,a)
9   format(a,f14.8)


  end subroutine readopt


!  *******************************************************************  !


END MODULE idsrdr_options

