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

  PUBLIC  :: spectralinit, spectral, freespectral, spctrl, dos
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
!  integer ienergy                      : Energy grid index             !
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

