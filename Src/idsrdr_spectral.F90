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
!                        MODULE idsrdr_spectral                         !
!  *******************************************************************  !
!  Description: computes the spectral function and the density of       !
!  states (DOS) of dynamic units (units where electron-phonon           !
!  interaction is considered).                                          !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *******************************************************************  !

MODULE idsrdr_spectral

!
!   Modules
!
  use parallel,        only: 
  use idsrdr_options,  only: 
  use idsrdr_engrid,   only: 
  use idsrdr_ephcoupl, only: 
  use idsrdr_units,    only: 
  use idsrdr_green,    only: 
  use idsrdr_iostream, only: 

  implicit none

  real(8), dimension(:,:,:), allocatable :: spctrl ! spectral function
  real(8), dimension(:,:,:), allocatable :: dos ! density of states

  PUBLIC  :: spectralinit, spectral, writespectral, freespectral
  PRIVATE :: calcspectral, calcspectraldisk


CONTAINS


!  *******************************************************************  !
!                             spectralinit                              !
!  *******************************************************************  !
!  Description: allocate spectral and density of states arrays.         !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer nspin               : Number of spin components              !
!  integer NTenerg_div         : Number of energy grid points per node  !
!  integer nunitseph           : Number of units with eph               !
!  *******************************************************************  !
  subroutine spectralinit

!
!   Modules
!
    use idsrdr_options,  only: nspin
    use idsrdr_engrid,   only: NTenerg_div
    use idsrdr_units,    only: nunitseph

!   Allocate and initializes spectral and density of states arrays.
    allocate (spctrl(NTenerg_div,nspin,nunitseph+1))
    allocate (dos(NTenerg_div,nspin,nunitseph+1))
    spctrl = 0.d0
    dos = 0.d0


  end subroutine spectralinit


!  *******************************************************************  !
!                               spectral                                !
!  *******************************************************************  !
!  Description: call appropriate subroutine for computing the spectral  !
!  function and the density of states (DOS).                            !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode                       : True if it is the I/O node    !
!  logical writeondisk                  : Write GFs on disk?            !
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 ienergy                       : Energy grid index             !
!  *******************************************************************  !
  subroutine spectral (ienergy, ispin)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_options,  only: writeondisk

!   Input variables.
    integer, intent(in) :: ienergy, ispin

    if (IOnode) write (6,'(a)', advance='no')                           &
            '      computing spectral function... '

    if (writeondisk) then
       call calcspectraldisk (ienergy, ispin)
    else
       call calcspectral (ienergy, ispin)
    endif

    if (IOnode) write(6,'(a)') " ok!"


  end subroutine spectral


!  *******************************************************************  !
!                             calcspectral                              !
!  *******************************************************************  !
!  Description: computes the spectral function and the density of       !
!  states (DOS).                                                        !
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
!  integer ephIdx(ntypeunits+2)         : Unit index (those with e-ph)  !
!  integer norbDyn(neph)                : Number of orbitals from       !
!                                         dynamic atoms                 !
!  integer idxF(neph)                   : First dynamic atom orbital    !
!  integer idxL(neph)                   : Last dynamic atom orbital     !
!  integer unit_type(nunits+2)          : Units types                   !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  TYPE(unitS) Sunits(ntypeunits+2)%S(unitdimensions,unitdimensions) :  !
!                                         [real*8] Units overlap        !
!  integer nunitseph                    : Number of units with eph      !
!  integer eph_type(nunitseph)          : Units types with eph          !
!  TYPE(green) Gr_nn(neph)%G(unitdimensions,unitdimensions) :           !
!                                                  [complex] G^r_{n,n}  !
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 ienergy                       : Energy grid index             !
!  *******************************************************************  !
  subroutine calcspectral (ienergy, ispin)

!
!   Modules
!
    use idsrdr_options,  only: ntypeunits, nunits
    use idsrdr_ephcoupl, only: ephIdx, norbDyn, idxF, idxL
    use idsrdr_units,    only: unit_type, unitdimensions, Sunits,       &
                               nunitseph, eph_type
    use idsrdr_green,    only: Gr_nn

!   Input variables.
    integer, intent(in) :: ienergy, ispin

!   Local variables.
    integer :: I, K, idx, ueph
    complex(8), dimension(:,:), allocatable :: Aux1, Aux2, cpS
    external :: zsymm

!   Initialize variables.
    idx = ephIdx(ntypeunits+1)
    ueph = 1 ! eph units indexing

    if (idx /= 0) then

!      Allocate auxiliary matrix.
       allocate (Aux1(norbDyn(idx),norbDyn(idx)))
       allocate (Aux2(norbDyn(idx),norbDyn(idx)))
       allocate (cpS(norbDyn(idx),norbDyn(idx)))

!      Copy dynamic orbitals part of overlap matrix.
       cpS = Sunits(unit_type(ntypeunits+1))                            &
             %S(idxF(idx):idxL(idx),idxF(idx):idxL(idx))

!      ('Aux1 = Gr_nn * Saux')
       call zsymm ('R', 'L', norbDyn(idx), norbDyn(idx),                &
                   (1.D0,0.D0), cpS, norbDyn(idx), Gr_nn(ueph)%G,       &
                   norbDyn(idx), (0.D0,0.D0), Aux1, norbDyn(idx))

!      ('Aux2 = Saux * Aux1')
       call zsymm ('L', 'L', norbDyn(idx), norbDyn(idx),                &
                   (1.D0,0.D0), cpS, norbDyn(idx), Aux1,                &
                   norbDyn(idx), (0.D0,0.D0), Aux2, norbDyn(idx))

       spctrl(ienergy,ispin,ueph) = 0.D0
       dos(ienergy,ispin,ueph) = 0.D0
       do K = 1,norbDyn(idx)
          spctrl(ienergy,ispin,ueph) = spctrl(ienergy,ispin,ueph)       &
                                       - DIMAG(Gr_nn(ueph)%G(K,K))
          dos(ienergy,ispin,ueph) = dos(ienergy,ispin,ueph)             &
                                    - DIMAG(Aux2(K,K))
       enddo
       spctrl(ienergy,ispin,nunitseph+1) = spctrl(ienergy,ispin,ueph)
       dos(ienergy,ispin,nunitseph+1) = dos(ienergy,ispin,ueph)

!      Free memory.
       deallocate (cpS)
       deallocate (Aux1)
       deallocate (Aux2)

       ueph = ueph + 1 ! eph units indexing

    endif

    do I = 2,nunits+1 ! over units from left to right

       idx = ephIdx(unit_type(I))

       if (idx /= 0) then

!         Allocate auxiliary matrix.
          allocate (cpS(norbDyn(idx),norbDyn(idx)))
          allocate (Aux1(norbDyn(idx),norbDyn(idx)))
          allocate (Aux2(norbDyn(idx),norbDyn(idx)))

!         Copy dynamic orbitals part of overlap matrix.
          cpS = Sunits(unit_type(I))%S(idxF(idx):idxL(idx),             &
                                       idxF(idx):idxL(idx))

!         ('Aux1 = GrMM * Saux')
          call zsymm ('R', 'L', norbDyn(idx), norbDyn(idx),             &
                      (1.D0,0.D0), cpS, norbDyn(idx), Gr_nn(ueph)%G,    &
                      norbDyn(idx), (0.D0,0.D0), Aux1, norbDyn(idx))

!         ('Aux2 = Saux * Aux1')
          call zsymm ('L', 'L', norbDyn(idx), norbDyn(idx),             &
                      (1.D0,0.D0), cpS, norbDyn(idx), Aux1,             &
                      norbDyn(idx), (0.D0,0.D0), Aux2, norbDyn(idx))

          spctrl(ienergy,ispin,ueph) = 0.D0
          dos(ienergy,ispin,ueph) = 0.D0
          do K = 1,norbDyn(idx)
             spctrl(ienergy,ispin,ueph) = spctrl(ienergy,ispin,ueph)    &
                                          - DIMAG(Gr_nn(ueph)%G(K,K))
             dos(ienergy,ispin,ueph) = dos(ienergy,ispin,ueph)          &
                                       - DIMAG(Aux2(K,K))
          enddo
          spctrl(ienergy,ispin,nunitseph+1) =                           &
                                      spctrl(ienergy,ispin,nunitseph+1) &
                                      + spctrl(ienergy,ispin,ueph)
          dos(ienergy,ispin,nunitseph+1) =                              &
                                         dos(ienergy,ispin,nunitseph+1) &
                                         + dos(ienergy,ispin,ueph)

!         Free memory.
          deallocate (cpS)
          deallocate (Aux1)
          deallocate (Aux2)

          ueph = ueph + 1 ! eph units indexing

       endif

    enddo

    idx = ephIdx(ntypeunits+2)
    if (idx /= 0) then

!      Allocate auxiliary matrix.
       allocate (cpS(norbDyn(idx),norbDyn(idx)))
       allocate (Aux1(norbDyn(idx),norbDyn(idx)))
       allocate (Aux2(norbDyn(idx),norbDyn(idx)))

!      Copy dynamic orbitals part of overlap matrix.
       cpS = Sunits(ntypeunits+2)%S(idxF(idx):idxL(idx),                &
                                    idxF(idx):idxL(idx))

!      ('Aux1 = GrMM * Saux')
       call zsymm ('R', 'L', norbDyn(idx), norbDyn(idx),                &
                   (1.D0,0.D0), cpS, norbDyn(idx), Gr_nn(ueph)%G,       &
                   norbDyn(idx), (0.D0,0.D0), Aux1, norbDyn(idx))

!      ('Aux2 = Saux * Aux1')
       call zsymm ('L', 'L', norbDyn(idx), norbDyn(idx),                &
                   (1.D0,0.D0), cpS, norbDyn(idx), Aux1,                &
                   norbDyn(idx), (0.D0,0.D0), Aux2, norbDyn(idx))

       spctrl(ienergy,ispin,ueph) = 0.D0
       dos(ienergy,ispin,ueph) = 0.D0
       do K = 1,norbDyn(idx)
          spctrl(ienergy,ispin,ueph) = spctrl(ienergy,ispin,ueph)       &
                                       - DIMAG(Gr_nn(ueph)%G(K,K))
          dos(ienergy,ispin,ueph) = dos(ienergy,ispin,ueph)             &
                                    - DIMAG(Aux2(K,K))
       enddo
       spctrl(ienergy,ispin,nunitseph+1) =                              &
                                      spctrl(ienergy,ispin,nunitseph+1) &
                                      + spctrl(ienergy,ispin,ueph)
       dos(ienergy,ispin,nunitseph+1) = dos(ienergy,ispin,nunitseph+1)  &
                                        + dos(ienergy,ispin,ueph)

!      Free memory.
       deallocate (cpS)
       deallocate (Aux1)
       deallocate (Aux2)

    endif


  end subroutine calcspectral


!  *******************************************************************  !
!                           calcspectraldisk                            !
!  *******************************************************************  !
!  Description: computes the spectral function and the density of       !
!  states (DOS) (Green's functions are read from disk).                 !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer ntypeunits                   : Number of unit types          !
!  integer nunits                       : Total number of units         !
!  integer ephIdx(ntypeunits+2)         : Unit index (those with e-ph)  !
!  integer norbDyn(neph)                : Number of orbitals from       !
!                                         dynamic atoms                 !
!  integer idxF(neph)                   : First dynamic atom orbital    !
!  integer idxL(neph)                   : Last dynamic atom orbital     !
!  integer unit_type(nunits+2)          : Units types                   !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  TYPE(unitS) Sunits(ntypeunits+2)%S(unitdimensions,unitdimensions) :  !
!                                         [real*8] Units overlap        !
!  integer nunitseph                    : Number of units with eph      !
!  integer eph_type(nunitseph)          : Units types with eph          !
!  TYPE(green) Gr_nn(neph)%G(unitdimensions,unitdimensions) :           !
!                                                  [complex] G^r_{n,n}  !
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 ienergy                       : Energy grid index             !
!  *******************************************************************  !
  subroutine calcspectraldisk (ienergy, ispin)

!
!   Modules
!
    use idsrdr_options,  only: ntypeunits, nunits
    use idsrdr_ephcoupl, only: ephIdx, norbDyn, idxF, idxL
    use idsrdr_units,    only: unit_type, unitdimensions, Sunits,       &
                               nunitseph, eph_type
    use idsrdr_green,    only: Gr_nn_disk, greenload
    use idsrdr_iostream, only: openstream, closestream

!   Input variables.
    integer, intent(in) :: ienergy, ispin

!   Local variables.
    integer :: I, K, idx, ueph
    complex(8), dimension(:,:), allocatable :: Aux1, Aux2, cpS, cpGr
    external :: zsymm

!   Initialize variables.
    idx = ephIdx(ntypeunits+1)
    ueph = 1 ! eph units indexing

!   Open Green's functions files.
    call openstream (Gr_nn_disk%fname, Gr_nn_disk%lun)

    if (idx /= 0) then

!      Allocate auxiliary matrix.
       allocate (Aux1(norbDyn(idx),norbDyn(idx)))
       allocate (Aux2(norbDyn(idx),norbDyn(idx)))
       allocate (cpS(norbDyn(idx),norbDyn(idx)))
       allocate (cpGr(norbDyn(idx),norbDyn(idx)))

!      Copy dynamic orbitals part of overlap matrix.
       cpS = Sunits(unit_type(ntypeunits+1))                            &
             %S(idxF(idx):idxL(idx),idxF(idx):idxL(idx))

!      Copy 'Gr_nn' from file to auxiliary matrix.
       call greenload (Gr_nn_disk%lun, norbDyn(idx), norbDyn(idx),      &
                       1, norbDyn(idx), 1, norbDyn(idx),                &
                       cpGr, norbDyn(idx), norbDyn(idx))

!      ('Aux1 = Gr_nn * Saux')
       call zsymm ('R', 'L', norbDyn(idx), norbDyn(idx),                &
                   (1.D0,0.D0), cpS, norbDyn(idx), cpGr,                &
                   norbDyn(idx), (0.D0,0.D0), Aux1, norbDyn(idx))

!      ('Aux2 = Saux * Aux1')
       call zsymm ('L', 'L', norbDyn(idx), norbDyn(idx),                &
                   (1.D0,0.D0), cpS, norbDyn(idx), Aux1,                &
                   norbDyn(idx), (0.D0,0.D0), Aux2, norbDyn(idx))

       spctrl(ienergy,ispin,ueph) = 0.D0
       dos(ienergy,ispin,ueph) = 0.D0
       do K = 1,norbDyn(idx)
          spctrl(ienergy,ispin,ueph) = spctrl(ienergy,ispin,ueph)       &
                                       - DIMAG(cpGr(K,K))
          dos(ienergy,ispin,ueph) = dos(ienergy,ispin,ueph)             &
                                    - DIMAG(Aux2(K,K))
       enddo
       spctrl(ienergy,ispin,nunitseph+1) = spctrl(ienergy,ispin,ueph)
       dos(ienergy,ispin,nunitseph+1) = dos(ienergy,ispin,ueph)

!      Free memory.
       deallocate (Aux1)
       deallocate (Aux2)
       deallocate (cpS)
       deallocate (cpGr)

       ueph = ueph + 1 ! eph units indexing

    endif

    do I = 2,nunits+1 ! over units from left to right

       idx = ephIdx(unit_type(I))

       if (idx /= 0) then

!         Allocate auxiliary matrix.
          allocate (Aux1(norbDyn(idx),norbDyn(idx)))
          allocate (Aux2(norbDyn(idx),norbDyn(idx)))
          allocate (cpS(norbDyn(idx),norbDyn(idx)))
          allocate (cpGr(norbDyn(idx),norbDyn(idx)))

!         Copy dynamic orbitals part of overlap matrix.
          cpS = Sunits(unit_type(I))%S(idxF(idx):idxL(idx),             &
                                       idxF(idx):idxL(idx))

!         Copy 'Gr_nn' from file to auxiliary matrix.
          call greenload (Gr_nn_disk%lun, norbDyn(idx), norbDyn(idx),   &
                          1, norbDyn(idx), 1, norbDyn(idx),             &
                          cpGr, norbDyn(idx), norbDyn(idx))

!         ('Aux1 = GrMM * Saux')
          call zsymm ('R', 'L', norbDyn(idx), norbDyn(idx),             &
                      (1.D0,0.D0), cpS, norbDyn(idx), cpGr,             &
                      norbDyn(idx), (0.D0,0.D0), Aux1, norbDyn(idx))

!         ('Aux2 = Saux * Aux1')
          call zsymm ('L', 'L', norbDyn(idx), norbDyn(idx),             &
                      (1.D0,0.D0), cpS, norbDyn(idx), Aux1,             &
                      norbDyn(idx), (0.D0,0.D0), Aux2, norbDyn(idx))

          spctrl(ienergy,ispin,ueph) = 0.D0
          dos(ienergy,ispin,ueph) = 0.D0
          do K = 1,norbDyn(idx)
             spctrl(ienergy,ispin,ueph) = spctrl(ienergy,ispin,ueph)    &
                                          - DIMAG(cpGr(K,K))
             dos(ienergy,ispin,ueph) = dos(ienergy,ispin,ueph)          &
                                       - DIMAG(Aux2(K,K))
          enddo
          spctrl(ienergy,ispin,nunitseph+1) =                           &
                                      spctrl(ienergy,ispin,nunitseph+1) &
                                      + spctrl(ienergy,ispin,ueph)
          dos(ienergy,ispin,nunitseph+1) =                              &
                                         dos(ienergy,ispin,nunitseph+1) &
                                         + dos(ienergy,ispin,ueph)

!         Free memory.
          deallocate (Aux1)
          deallocate (Aux2)
          deallocate (cpS)
          deallocate (cpGr)

          ueph = ueph + 1 ! eph units indexing

       endif

    enddo

    idx = ephIdx(ntypeunits+2)
    if (idx /= 0) then

!      Allocate auxiliary matrix.
       allocate (Aux1(norbDyn(idx),norbDyn(idx)))
       allocate (Aux2(norbDyn(idx),norbDyn(idx)))
       allocate (cpS(norbDyn(idx),norbDyn(idx)))
       allocate (cpGr(norbDyn(idx),norbDyn(idx)))

!      Copy dynamic orbitals part of overlap matrix.
       cpS = Sunits(ntypeunits+2)%S(idxF(idx):idxL(idx),                &
                                    idxF(idx):idxL(idx))

!      Copy 'Gr_nn' from file to auxiliary matrix.
       call greenload (Gr_nn_disk%lun, norbDyn(idx), norbDyn(idx),      &
                       1, norbDyn(idx), 1, norbDyn(idx),                &
                       cpGr, norbDyn(idx), norbDyn(idx))

!      ('Aux1 = GrMM * Saux')
       call zsymm ('R', 'L', norbDyn(idx), norbDyn(idx),                &
                   (1.D0,0.D0), cpS, norbDyn(idx), cpGr,                &
                   norbDyn(idx), (0.D0,0.D0), Aux1, norbDyn(idx))

!      ('Aux2 = Saux * Aux1')
       call zsymm ('L', 'L', norbDyn(idx), norbDyn(idx),                &
                   (1.D0,0.D0), cpS, norbDyn(idx), Aux1,                &
                   norbDyn(idx), (0.D0,0.D0), Aux2, norbDyn(idx))

       spctrl(ienergy,ispin,ueph) = 0.D0
       dos(ienergy,ispin,ueph) = 0.D0
       do K = 1,norbDyn(idx)
          spctrl(ienergy,ispin,ueph) = spctrl(ienergy,ispin,ueph)       &
                                       - DIMAG(cpGr(K,K))
          dos(ienergy,ispin,ueph) = dos(ienergy,ispin,ueph)             &
                                    - DIMAG(Aux2(K,K))
       enddo
       spctrl(ienergy,ispin,nunitseph+1) =                              &
                                      spctrl(ienergy,ispin,nunitseph+1) &
                                      + spctrl(ienergy,ispin,ueph)
       dos(ienergy,ispin,nunitseph+1) = dos(ienergy,ispin,nunitseph+1)  &
                                        + dos(ienergy,ispin,ueph)

!      Free memory.
       deallocate (Aux1)
       deallocate (Aux2)
       deallocate (cpS)
       deallocate (cpGr)

    endif

!   Close Green's functions files.
    call closestream (Gr_nn_disk%fname, Gr_nn_disk%lun)


  end subroutine calcspectraldisk


!  *******************************************************************  !
!                             writespectral                             !
!  *******************************************************************  !
!  Description: writes calculated spectral function and density of      !
!  states to output files ('slabel.SPCTR' and 'slabel.DOS').            !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer Node                : Actual node (MPI_Comm_rank)            !
!  integer Nodes               : Total number of nodes (MPI_Comm_size)  !
!  integer nspin               : Number of spin components              !
!  integer label_length        : Length of system label                 !
!  character(label_length) slabel : System Label (for output files)     !
!  character(60) directory     : Working directory                      !
!  integer NTenerg_div         : Number of energy grid points per node  !
!  real*8 Ei(NTenerg_div)      : Energy grid points                     !
!  integer nunitseph           : Number of units with eph               !
!  *******************************************************************  !
  subroutine writespectral

!
!   Modules
!
    use parallel,        only: IOnode, Node, Nodes
    use idsrdr_options,  only: nspin, label_length, slabel, directory
    use idsrdr_engrid,   only: NTenerg_div, Ei
    use idsrdr_units,    only: nunitseph

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    integer :: J, n, e, s, iuSpc, iuDos
    real(8), parameter :: pi = 3.1415926535897932384626433832795028841D0
    real(8), dimension(:), allocatable :: buffEn
    real(8), dimension(:,:), allocatable :: buffSpc, buffDos
    character(len=10) :: suffix
    character(len=label_length+70) :: fnSpc, fnDos
    character(len=label_length+70), external :: paste
    character(len=10), external :: pasbias2
    external :: io_assign, io_close
#ifdef MPI
    integer :: MPIerror
    integer, dimension(MPI_Status_Size) :: MPIstatus
#endif

    do J = 1,nunitseph+1 ! over units with e-ph intereaction

       if (IOnode) then

!         Set file's names and open it.
          write (suffix,'(i3)') J
          suffix = pasbias2 (suffix, '.SPCTR')
          suffix = paste ('_', suffix)
          fnSpc = paste (slabel, suffix)
          fnSpc = paste (directory, fnSpc)
          write (suffix,'(i3)') J
          suffix = pasbias2 (suffix, '.DOS')
          suffix = paste ('_', suffix)
          fnDos = paste (slabel, suffix)
          fnDos = paste (directory, fnDos)
          call io_assign (iuSpc)
          open (iuSpc, file=fnSpc, form='formatted', status='unknown')
          call io_assign (iuDos)
          open (iuDos, file=fnDos, form='formatted', status='unknown')

          write (6, '(/,a,i3,a,a)', advance='no')                       &
               'Writing spectral function ', J, ' to: ', trim(fnSpc)
          write (6, '(/,a,i3,a,a)') 'Writing DOS ', J, ' to: ',         &
               trim(fnDos)

!         Allocate buffers arrays.
          allocate (buffEn(NTenerg_div))
          allocate (buffSpc(NTenerg_div,nspin))
          allocate (buffDos(NTenerg_div,nspin))

       endif

!      Write to the output files (energies in eV (from CODATA - 2012)).
#ifdef MPI
       do n = 0,Nodes-1
          if (Node == n .and. Node == 0) then
#endif
             do e = 1,NTenerg_div
                write (iuSpc, '(/,e17.8e3)', advance='no')              &
!!$                     Ei(e) * 13.60569253D0 !eV
                     Ei(e) !Ry
                write (iuDos, '(/,e17.8e3)', advance='no')              &
!!$                     Ei(e) * 13.60569253D0 !eV
                     Ei(e) !Ry
                do s = 1,nspin
                   write (iuSpc, '(e17.8e3)', advance='no')             &
!!$                        spctrl(e,s,J) / (13.60569253D0 * pi) !eV
                        spctrl(e,s,J) !Ry
                   write (iuDos, '(e17.8e3)', advance='no')             &
!!$                        dos(e,s,J) / (13.60569253D0 * pi) !eV
                        dos(e,s,J) !Ry
                enddo
             enddo
#ifdef MPI
          elseif (Node == n) then
             call MPI_Send (Ei, NTenerg_div, MPI_Double_Precision,      &
                            0, 1, MPI_Comm_world, MPIerror)
             call MPI_Send (spctrl(1,1,J), NTenerg_div*nspin,           &
                            MPI_Double_Precision, 0, 2,                 &
                            MPI_Comm_world, MPIerror)
             call MPI_Send (dos(1,1,J), NTenerg_div*nspin,              &
                            MPI_Double_Precision, 0, 3,                 &
                            MPI_Comm_world, MPIerror)
          elseif (Node == 0) then
             call MPI_Recv (buffEn, NTenerg_div, MPI_Double_Precision,  &
                            n, 1, MPI_Comm_world, MPIstatus, MPIerror)
             call MPI_Recv (buffSpc, NTenerg_div*nspin,                 &
                            MPI_Double_Precision, n, 2,                 &
                            MPI_Comm_world, MPIstatus, MPIerror)
             call MPI_Recv (buffDos, NTenerg_div*nspin,                 &
                            MPI_Double_Precision, n, 3,                 &
                            MPI_Comm_world, MPIstatus, MPIerror)
          endif
          if (n /= 0) then
             call MPI_Barrier (MPI_Comm_world, MPIerror)
             if (Node == 0) then
                do e = 1,NTenerg_div
                   write (iuSpc, '(/,e17.8e3)', advance='no')           &
!!$                        buffEn(e) * 13.60569253D0 !eV
                        buffEn(e) !Ry
                   write (iuDos, '(/,e17.8e3)', advance='no')           &
!!$                        buffEn(e) * 13.60569253D0 !eV
                        buffEn(e) !Ry
                   do s = 1,nspin
                      write (iuSpc, '(e17.8e3)', advance='no')          &
!!$                           buffSpc(e,s) / (13.60569253D0 * pi) !eV
                           buffSpc(e,s) !Ry
                      write (iuDos, '(e17.8e3)', advance='no')          &
!!$                           buffDos(e,s) / (13.60569253D0 * pi) !eV
                           buffDos(e,s) !Ry
                   enddo
                enddo
             endif
          endif
       enddo
#endif

!      Close files and free buffers memory.
       if (Node.eq.0) then

          call io_close (iuSpc)
          call io_close (iuDos)

          deallocate (buffEn)
          deallocate (buffSpc)
          deallocate (buffDos)

       endif

    enddo ! do J = 1,nunitseph


  end subroutine writespectral


!  *******************************************************************  !
!                             freespectral                              !
!  *******************************************************************  !
!  Description: free allocated vectors.                                 !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *******************************************************************  !
  subroutine freespectral

!   Free memory.
    deallocate (spctrl)
    deallocate (dos)


  end subroutine freespectral


!  *******************************************************************  !


END MODULE idsrdr_spectral

