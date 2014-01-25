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
!                          MODULE idsrdr_green                          !
!  *******************************************************************  !
!  Description: calculate the required non-equilibrium green's          !
!  functions with recursive technique.                                  !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *******************************************************************  !

MODULE idsrdr_green

!
!   Modules
!
  use parallel,        only: 
  use idsrdr_options,  only: 
  use idsrdr_units,    only: 
  use idsrdr_leads,    only: 
  use idsrdr_ephcoupl, only: 
  use idsrdr_check,    only: 
  use idsrdr_io,       only: 
  use idsrdr_string,   only: 

  implicit none
  
  PUBLIC  :: greeninit, greenfunctions, freegreen,                      &
             Gr_nn, Gr_1n, Gr_Mn, Gr_1M,                                &
             Gr_nn_disk, Gr_1n_disk, Gr_Mn_disk, greenload
  PRIVATE :: LRsweep , RLsweep, GFfull, GFtest,                         &
             GL_mm, GL_1m, GR_pp, GR_Mp,                                &
             LRsweepDisk , RLsweepDisk, GFfullDisk, GFtestDisk,         &
             GL_mm_disk, GL_1m_disk, GR_pp_disk, GR_Mp_disk,            &
             npos, pos_pp, pos_Mp, greenFilesSet, greenloadR

! Type for storing Green's functions on memory.
  TYPE green
     complex(8), pointer :: G(:,:) ! pointer to Green's function matrix
  END TYPE green

  TYPE(green), allocatable, dimension (:) :: GL_mm ! G^L_{n-1,n-1}
  TYPE(green), allocatable, dimension (:) :: GL_1m ! G^L_{1,n-1}
  TYPE(green), allocatable, dimension (:) :: GR_pp ! G^R_{n+1,n+1}
  TYPE(green), allocatable, dimension (:) :: GR_Mp ! G^R_{M,n+1}
  TYPE(green), allocatable, dimension (:) :: Gr_nn ! G^r_{n,n}
  TYPE(green), allocatable, dimension (:) :: Gr_1n ! G^r_{1,n}
  TYPE(green), allocatable, dimension (:) :: Gr_Mn ! G^r_{M,n}

! Type for storing Green's functions on disk.
  TYPE greenFile
     integer :: lun ! logical unit number
     character(len=80) :: fname ! file name
  END TYPE greenfile

  TYPE(greenfile) :: GL_mm_disk ! G^L_{n-1,n-1}
  TYPE(greenfile) :: GL_1m_disk ! G^L_{1,n-1}
  TYPE(greenfile) :: GR_pp_disk ! G^R_{n+1,n+1}
  TYPE(greenfile) :: GR_Mp_disk ! G^R_{M,n+1}
  TYPE(greenfile) :: Gr_nn_disk ! G^r_{n,n}
  TYPE(greenfile) :: Gr_1n_disk ! G^r_{1,n}
  TYPE(greenfile) :: Gr_Mn_disk ! G^r_{M,n}

! Index array with the positions of written
! matrices (for reading 'GR_pp' and 'GR_Mp').
  integer :: npos
  integer, allocatable, dimension (:) :: pos_pp
  integer, allocatable, dimension (:) :: pos_Mp

  complex(8), allocatable, dimension (:,:) :: Gr_1M ! G^r_{1,M}

! The following are allocated only when the last unit doesn't have e-ph.
  complex(8), allocatable, dimension (:,:) :: GL_nn ! G^L_{n,n}
  complex(8), allocatable, dimension (:,:) :: GL_1N ! G^L_{1,M}


CONTAINS


!  *******************************************************************  !
!                               greeninit                               !
!  *******************************************************************  !
!  Description: allocate Green's functions matrices.                    !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer ntypeunits                   : Number of unit types          !
!  integer nunits                       : Total number of units         !
!  logical writeondisk                  : Write GFs on disk?            !
!  integer unit_type(nunits+2)          : Units types                   !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  integer nunitseph                   : Number of units with eph       !
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  integer ephIdx(ntypeunits+2)        : Unit index (those with e-ph)   !
!  integer norbDyn(neph)               : Number of orbitals from        !
!                                        dynamic atoms                  !
!  *******************************************************************  !
  subroutine greeninit

!
!   Modules
!
    use idsrdr_options,  only: ntypeunits, nunits, writeondisk
    use idsrdr_units,    only: unit_type, unitdimensions, nunitseph
    use idsrdr_leads,    only: NL, NR
    use idsrdr_ephcoupl, only: ephIdx, norbDyn

!   Local variables.
    integer :: I, dim, dimbfr, dimaft, idx, ueph

    IF (writeondisk) THEN

!      Set Green's functions files ("LUN" and file name).
       call greenFilesSet (141, 'GL_mm', GL_mm_disk)
       call greenFilesSet (142, 'GL_1m', GL_1m_disk)
       call greenFilesSet (143, 'GR_pp', GR_pp_disk)
       call greenFilesSet (144, 'GR_Mp', GR_Mp_disk)
       call greenFilesSet (145, 'Gr_nn', Gr_nn_disk)
       call greenFilesSet (146, 'Gr_1n', Gr_1n_disk)
       call greenFilesSet (147, 'Gr_Mn', Gr_Mn_disk)

!      Allocate index array (for reading 'GR_pp' and 'GR_Mp').
       if (ephIdx(ntypeunits+2) /= 0) then
          npos = nunitseph - 1
       else
          npos = nunitseph
       endif
       allocate (pos_pp(npos))
       allocate (pos_Mp(npos))

! DUVIDA: GL_1m tem que ter dimensão NL mesmo?
! E GR_Mp tem que ter dimensão NL mesmo?

!      Allocate Green's functions matrices.
       if (ephIdx(ntypeunits+2) == 0) then
          dimbfr = unitdimensions(unit_type(nunits+1))
          allocate (GL_nn(dimbfr,dimbfr))
          allocate (GL_1N(NL,dimbfr))
       endif
       allocate (Gr_1M(NL,NR))

    ELSE ! write on memory...

!      Allocate Green's functions pointer array.
       allocate (GL_mm(nunitseph))
       allocate (GL_1m(nunitseph))
       allocate (GR_pp(nunitseph))
       allocate (GR_Mp(nunitseph))
       allocate (Gr_nn(nunitseph))
       allocate (Gr_1n(nunitseph))
       allocate (Gr_Mn(nunitseph))

! DUVIDA: GL_1m tem que ter dimensão NL mesmo?
! E GR_Mp tem que ter dimensão NL mesmo?

!      Allocate Green's functions matrices.
       dim = unitdimensions(ntypeunits+1) ! current unit dimension
       dimbfr = NL ! last unit dimension
       dimaft = unitdimensions(unit_type(2)) ! next unit dimension
       idx = ephIdx(ntypeunits+1) ! e-ph unit type
       ueph = 1
       if (idx /= 0) then
          allocate (GL_mm(ueph)%G(dimbfr,dimbfr))
          allocate (GL_1m(ueph)%G(NL,dimbfr))
          allocate (GR_pp(ueph)%G(dimaft,dimaft))
          allocate (GR_Mp(ueph)%G(NR,dimaft))
          allocate (Gr_nn(ueph)%G(norbDyn(idx),norbDyn(idx)))
          allocate (Gr_1n(ueph)%G(NL,norbDyn(idx)))
          allocate (Gr_Mn(ueph)%G(NR,norbDyn(idx)))
          ueph = ueph + 1
       endif

       do I = 2,nunits
          dimbfr = dim ! last unit dimension
          dim = unitdimensions(unit_type(I)) ! current unit dimension
          dimaft = unitdimensions(unit_type(I+1)) ! next unit dimension
          idx = ephIdx(unit_type(I)) ! e-ph unit type
          if (idx /= 0) then
             allocate (GL_mm(ueph)%G(dimbfr,dimbfr))
             allocate (GL_1m(ueph)%G(NL,dimbfr))
             allocate (GR_pp(ueph)%G(dimaft,dimaft))
             allocate (GR_Mp(ueph)%G(NR,dimaft))
             allocate (Gr_nn(ueph)%G(norbDyn(idx),norbDyn(idx)))
             allocate (Gr_1n(ueph)%G(NL,norbDyn(idx)))
             allocate (Gr_Mn(ueph)%G(NR,norbDyn(idx)))
             ueph = ueph + 1
          endif
       enddo

       dimbfr = dim ! last unit dimension
       dim = unitdimensions(unit_type(I)) ! current unit dimension
       dimaft = unitdimensions(ntypeunits+2) ! next unit dimension
       idx = ephIdx(unit_type(I)) ! e-ph unit type
       if (idx /= 0) then
          allocate (GL_mm(ueph)%G(dimbfr,dimbfr))
          allocate (GL_1m(ueph)%G(NL,dimbfr))
          allocate (GR_pp(ueph)%G(dimaft,dimaft))
          allocate (GR_Mp(ueph)%G(NR,dimaft))
          allocate (Gr_nn(ueph)%G(norbDyn(idx),norbDyn(idx)))
          allocate (Gr_1n(ueph)%G(NL,norbDyn(idx)))
          allocate (Gr_Mn(ueph)%G(NR,norbDyn(idx)))
          ueph = ueph + 1
       endif

       dimbfr = dim ! last unit dimension
       dim = unitdimensions(ntypeunits+2) ! current unit dimension
       dimaft = NR ! next unit dimension
       idx = ephIdx(ntypeunits+2) ! e-ph unit type
       if (idx /= 0) then
          allocate (GL_mm(ueph)%G(dimbfr,dimbfr))
          allocate (GL_1m(ueph)%G(NL,dimbfr))
          allocate (GR_pp(ueph)%G(dimaft,dimaft))
          allocate (GR_Mp(ueph)%G(NR,dimaft))
          allocate (Gr_nn(ueph)%G(norbDyn(idx),norbDyn(idx)))
          allocate (Gr_1n(ueph)%G(NL,norbDyn(idx)))
          allocate (Gr_Mn(ueph)%G(NR,norbDyn(idx)))
       else
          allocate (GL_nn(dimbfr,dimbfr))
          allocate (GL_1N(NL,dimbfr))
       endif
       allocate (Gr_1M(NL,NR))

    ENDIF ! IF (writeondisk)


  end subroutine greeninit


!  *******************************************************************  !
!                            greenfunctions                             !
!  *******************************************************************  !
!  Description: main subroutine that computes the required Green's      !
!  functions.                                                           !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical writeondisk                  : Write GFs on disk?            !
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 Ei                            : Energy grid point             !
!  *******************************************************************  !
  subroutine greenfunctions (Ei, ispin)

!
!   Modules
!
    use idsrdr_options,  only: writeondisk

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

    IF (writeondisk) THEN

!      Perform the left-to-right sweep.
       call LRsweepDisk (Ei, ispin)

!      Perform the right-to-left sweep.
       call RLsweepDisk (Ei, ispin)

!      Compute the required full Green's functions.
       call GFfullDisk (Ei, ispin)

#ifdef DEBUG
!      Compute the entire Green's function of scattering region.
       call GFtestDisk (Ei, ispin)
#endif

    ELSE ! write on memory...

!      Perform the left-to-right sweep.
       call LRsweep (Ei, ispin)

!      Perform the right-to-left sweep.
       call RLsweep (Ei, ispin)

!      Compute the required full Green's functions.
       call GFfull (Ei, ispin)

#ifdef DEBUG
!      Compute the entire Green's function of scattering region.
       call GFtest (Ei, ispin)
#endif

    ENDIF ! IF (writeondisk)


  end subroutine greenfunctions


!  *******************************************************************  !
!                                LRsweep                                !
!  *******************************************************************  !
!  Description: left-to-right sweep for computing 'G^L_{n-1,n-1}' and   !
!  'G^L_{1,n-1}'.                                                       !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode                       : True if it is the I/O node    !
!  integer ntypeunits                   : Number of unit types          !
!  integer nunits                       : Total number of units         !
!  integer unit_type(nunits+2)          : Units types                   !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  real*8 unitshift(ntypeunits+2)       : Units shift                   !
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
!  integer ephIndic(ntypeunits+2)        : E-ph interaction in          !
!                                          unit type (0 or 1)           !
!  integer NL                           : Number of left lead orbitals  !
!  complex(8) Sigma_L(NL,NL)            : Left-lead self-energy         !
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 Ei                            : Energy grid point             !
!  ***************************** OUTPUT ******************************  !
!  TYPE(green) GL_mm(neph)%G(unitdimensions,unitdimensions) :           !
!                                              [complex] G^L_{n-1,n-1}  !
!  TYPE(green) GL_1m(neph)%G(NL,unitdimensions) :                       !
!                                              [complex] G^L_{1,n-1}    !
!  *******************************************************************  !
  subroutine LRsweep (Ei, ispin)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_options,  only: ntypeunits, nunits
    use idsrdr_units,    only: unit_type, unitdimensions, unitshift,    &
                               S1unit, H1unit, Sunits, Hunits, ephIndic
    use idsrdr_leads,    only: NL, Sigma_L
    use idsrdr_check,    only: CHECKzsytrf, CHECKzsytri

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    integer :: I, n, utype, dim, dimbfr, ueph
    integer, allocatable, dimension (:) :: ipiv
    complex(8), allocatable, dimension (:,:) :: V ! pristine coupling
    complex(8), allocatable, dimension (:,:) :: Gbfr ! GF before
    complex(8), allocatable, dimension (:,:) :: Gaft ! GF after
    complex(8), allocatable, dimension (:,:) :: Gbfr_1m ! GF_1m before
    complex(8), allocatable, dimension (:,:) :: Gaft_1m ! GF_1m after
    complex(8), allocatable, dimension (:,:) :: aux1, aux2
    complex(8), allocatable, dimension (:,:) :: foo1, foo2, foo3
    external :: zsymm, zgemm

    if (IOnode) write (6,'(a)', advance='no')                           &
            '      computing left-to-right sweep... '

!   Initialize variables.
    utype = ntypeunits + 1 ! current unit type (first unit)
    dim = unitdimensions(utype) ! current type dimension
    n = unitdimensions(ntypeunits) ! pristine dimension
    if (ephIndic(utype)) then ! eph units indexing
       ueph = 2
    else
       ueph = 1
    endif

!   Allocate matrices.
    allocate (V(n,n))
    allocate (Gbfr(dim,dim))
    allocate (Gbfr_1m(dim,dim))
    allocate (ipiv(dim))
    allocate (aux1(n,n))
    allocate (aux2(n,n))
    allocate (foo1(NL,n))
    allocate (foo2(NL,n))

!   Coupling matrix (pristine).
    V = (Ei-unitshift(ntypeunits))*S1unit - H1unit(:,:,ispin)

!   Compute the first step.
    Gbfr = (Ei-unitshift(utype))*Sunits(utype)%S                        &
           - Hunits(utype)%H(:,:,ispin)
    Gbfr(1:NL,1:NL) = Gbfr(1:NL,1:NL) - Sigma_L
    call CHECKzsytrf (dim, 'L', Gbfr, ipiv)
    call CHECKzsytri (dim, 'L', Gbfr, ipiv)
    Gbfr_1m = Gbfr

!   Free memory.
    deallocate (ipiv)

!   Loop over the blocks from left to right.
    do I = 2,nunits+1

!      Assign auxiliary variables.
       dimbfr = dim ! old dimension
       utype = unit_type(I) ! current unit type
       dim = unitdimensions(utype) ! current type dimension

!      Store matrix if required.
       if (ephIndic(utype)) then
          GL_mm(ueph)%G = Gbfr
          GL_1m(ueph)%G = Gbfr_1m
          ueph = ueph + 1
       endif

!      ('aux2 = Gbfr*V')
       aux1 = Gbfr(dimbfr-n+1:dimbfr,dimbfr-n+1:dimbfr)
       call zsymm ('L', 'L', n, n, (1.d0,0.d0), aux1, n,                &
                   V, n, (0.d0,0.d0), aux2, n)

!      ('aux1 = V^dagger*aux2')
       call zgemm ('C', 'N', n, n, n, (1.d0,0.d0), V, n,                &
                   aux2, n, (0.d0,0.d0), aux1, n)

!      (Re-)Allocate matrices and array.
       allocate (Gaft(dim,dim))
       allocate (ipiv(dim))

!      ('Gaft = (E*S - H - V^T*Gbfr*V)^-1')
       Gaft = (Ei-unitshift(utype))*Sunits(utype)%S                     &
              - Hunits(utype)%H(:,:,ispin)
       Gaft(1:n,1:n) = Gaft(1:n,1:n) - aux1
       call CHECKzsytrf (dim, 'L', Gaft, ipiv)
       call CHECKzsytri (dim, 'L', Gaft, ipiv)

!      (Re-)Allocate matrices.
       allocate (Gaft_1m(NL,dim))
       allocate (foo3(n,dim))

!      ('foo2 = Gbfr_1m*V')
       foo1 = Gbfr_1m(1:NL,dimbfr-n+1:dimbfr)
       call zgemm ('N', 'N', NL, n, n, (1.d0,0.d0), foo1, NL,           &
                   V, n, (0.d0,0.d0), foo2, NL)

!      ('Gaft_1m = - foo2 * Gaft')
       foo3 = Gaft(1:n,:)
       call zgemm ('N', 'N', NL, dim, n, (-1.d0,0.d0), foo2, NL,        &
                   foo3, n, (0.d0,0.d0), Gaft_1m, NL)

!      Reallocate matrix.
       deallocate (Gbfr)
       deallocate (Gbfr_1m)
       allocate (Gbfr(dim,dim))
       allocate (Gbfr_1m(NL,dim))

       Gbfr = Gaft
       Gbfr_1m = Gaft_1m

!      (Re-)Free memory.
       deallocate (Gaft)
       deallocate (Gaft_1m)
       deallocate (foo3)
       deallocate (ipiv)

    enddo

!   Store matrix if required.
    if (ephIndic(ntypeunits+2)) then
       GL_mm(ueph)%G = Gbfr
       GL_1m(ueph)%G = Gbfr_1m
    else
       GL_nn = Gbfr
       GL_1N = Gbfr_1m
    endif

!   Free memory.
    deallocate (V)
    deallocate (Gbfr)
    deallocate (Gbfr_1m)
    deallocate (aux1)
    deallocate (aux2)
    deallocate (foo1)
    deallocate (foo2)

    if (IOnode) write(6,'(a)') " ok!"


  end subroutine LRsweep


!  *******************************************************************  !
!                                RLsweep                                !
!  *******************************************************************  !
!  Description: right-to-left sweep for computing 'G^R_{n+1,n+1}' and   !
!  'G^R_{M,n+1}'.                                                       !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode                       : True if it is the I/O node    !
!  integer ntypeunits                   : Number of unit types          !
!  integer nunits                       : Total number of units         !
!  integer unit_type(nunits+2)          : Units types                   !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  real*8 unitshift(ntypeunits+2)       : Units shift                   !
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
!  integer ephIndic(ntypeunits+2)        : E-ph interaction in          !
!                                          unit type (0 or 1)           !
!  integer nunitseph                   : Number of units with eph       !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Sigma_R(NR,NR)           : Right-lead self-energy         !
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 Ei                            : Energy grid point             !
!  ***************************** OUTPUT ******************************  !
!  TYPE(green) GR_pp(neph)%G(unitdimensions,unitdimensions) :           !
!                                              [complex] G^R_{n+1,n+1}  !
!  TYPE(green) GR_Mp(neph)%G(NR,unitdimensions) :                       !
!                                              [complex] G^R_{M,n+1}    !
!  *******************************************************************  !
  subroutine RLsweep (Ei, ispin)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_options,  only: ntypeunits, nunits
    use idsrdr_units,    only: unit_type, unitdimensions, unitshift,    &
                               S1unit, H1unit, Sunits, Hunits,          &
                               ephIndic, nunitseph
    use idsrdr_leads,    only: NR, Sigma_R
    use idsrdr_check,    only: CHECKzsytrf, CHECKzsytri

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    integer :: I, n, utype, dim, ueph
    integer, allocatable, dimension (:) :: ipiv
    complex(8), allocatable, dimension (:,:) :: V ! pristine coupling
    complex(8), allocatable, dimension (:,:) :: Gbfr ! GF before
    complex(8), allocatable, dimension (:,:) :: Gaft ! GF after
    complex(8), allocatable, dimension (:,:) :: Gbfr_Mp ! GF_Mp before
    complex(8), allocatable, dimension (:,:) :: Gaft_Mp ! GF_Mp after
    complex(8), allocatable, dimension (:,:) :: aux1, aux2
    complex(8), allocatable, dimension (:,:) :: foo1, foo2, foo3
    external :: zsymm, zgemm

    if (IOnode) write (6,'(a)', advance='no')                           &
            '      computing right-to-left sweep... '

!   Initialize variables.
    utype = ntypeunits + 2 ! current unit type (last unit)
    dim = unitdimensions(utype) ! current type dimension
    n = unitdimensions(ntypeunits) ! pristine dimension
    if (ephIndic(utype)) then ! eph units indexing
       ueph = nunitseph - 1
    else
       ueph = nunitseph
    endif

!   Allocate matrices.
    allocate (V(n,n))
    allocate (Gbfr(dim,dim))
    allocate (Gbfr_Mp(dim,dim))
    allocate (ipiv(dim))
    allocate (aux1(n,n))
    allocate (aux2(n,n))
    allocate (foo1(NR,n))
    allocate (foo2(NR,n))

!   Coupling matrix (pristine).
    V = (Ei-unitshift(ntypeunits))*S1unit - H1unit(:,:,ispin)

!   Compute the first step.
    Gbfr = (Ei-unitshift(utype))*Sunits(utype)%S                        &
           - Hunits(utype)%H(:,:,ispin)
    Gbfr(dim-NR+1:dim,dim-NR+1:dim) = Gbfr(dim-NR+1:dim,dim-NR+1:dim)   &
                                      - Sigma_R
    call CHECKzsytrf (dim, 'L', Gbfr, ipiv)
    call CHECKzsytri (dim, 'L', Gbfr, ipiv)
    Gbfr_Mp = Gbfr

!   Free memory.
    deallocate (ipiv)

!   Loop over the blocks from right to left.
    do I = nunits+1,2,-1

!      Assign auxiliary variables.
       utype = unit_type(I) ! current unit type
       dim = unitdimensions(utype) ! current type dimension

!      Store matrix if required.
       if (ephIndic(utype)) then
          GR_pp(ueph)%G = Gbfr
          GR_Mp(ueph)%G = Gbfr_Mp
          ueph = ueph - 1
       endif

!      ('aux2 = V*Gbfr')
       aux1 = Gbfr(1:n,1:n)
       call zsymm ('R', 'L', n, n, (1.d0,0.d0), aux1, n,                &
                   V, n, (0.d0,0.d0), aux2, n)

!      ('aux1 = aux2*V^dagger')
       call zgemm ('N', 'C', n, n, n, (1.d0,0.d0), aux2, n,             &
                   V, n, (0.d0,0.d0), aux1, n)

!      (Re-)Allocate matrices and array.
       allocate (Gaft(dim,dim))
       allocate (ipiv(dim))

!      ('Gaft = (E*S - H - V^T*Gbfr*V)^-1')
       Gaft = (Ei-unitshift(utype))*Sunits(utype)%S                     &
              - Hunits(utype)%H(:,:,ispin)
       Gaft(dim-n+1:dim,dim-n+1:dim) = Gaft(dim-n+1:dim,dim-n+1:dim)    &
                                       - aux1
       call CHECKzsytrf (dim, 'L', Gaft, ipiv)
       call CHECKzsytri (dim, 'L', Gaft, ipiv)

!      (Re-)Allocate matrices.
       allocate (Gaft_Mp(NR,dim))
       allocate (foo3(n,dim))

!      ('foo2 = Gbfr_Mp*V^dagger')
       foo1 = Gbfr_Mp(1:NR,1:n)
       call zgemm ('N', 'C', NR, n, n, (1.d0,0.d0), foo1, NR,           &
                   V, n, (0.d0,0.d0), foo2, NR)

!      ('Gaft_Mp = - foo2 * Gaft')
       foo3 = Gaft(dim-n+1:dim,:)
       call zgemm ('N', 'N', NR, dim, n, (-1.d0,0.d0), foo2, NR,        &
                   foo3, n, (0.d0,0.d0), Gaft_Mp, NR)

!      Reallocate matrix.
       deallocate (Gbfr)
       deallocate (Gbfr_Mp)
       allocate (Gbfr(dim,dim))
       allocate (Gbfr_Mp(NR,dim))

       Gbfr = Gaft
       Gbfr_Mp = Gaft_Mp

!      (Re-)Free memory.
       deallocate (Gaft)
       deallocate (Gaft_Mp)
       deallocate (foo3)
       deallocate (ipiv)

    enddo

!   Store matrix if required.
    if (ephIndic(ntypeunits+1)) then
       GR_pp(ueph)%G = Gbfr
       GR_Mp(ueph)%G = Gbfr_Mp
    endif

!   Free memory.
    deallocate (V)
    deallocate (Gbfr)
    deallocate (Gbfr_Mp)
    deallocate (aux1)
    deallocate (aux2)
    deallocate (foo1)
    deallocate (foo2)

    if (IOnode) write(6,'(a)') " ok!"


  end subroutine RLsweep


!  *******************************************************************  !
!                                GFfull                                 !
!  *******************************************************************  !
!  Description: computes the full Green's function of required units.   !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode                       : True if it is the I/O node    !
!  integer ntypeunits                   : Number of unit types          !
!  integer nunits                       : Total number of units         !
!  integer unit_type(nunits+2)          : Units types                   !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  real*8 unitshift(ntypeunits+2)       : Units shift                   !
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
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Sigma_L(NL,NL)           : Left-lead self-energy          !
!  complex(8) Sigma_R(NR,NR)           : Right-lead self-energy         !
!  integer ephIdx(ntypeunits+2)        : Unit index (those with e-ph)   !
!  integer norbDyn(neph)               : Number of orbitals from        !
!                                        dynamic atoms                  !
!  integer idxF(neph)                  : First dynamic atom orbital     !
!  integer idxL(neph)                  : Last dynamic atom orbital      !
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 Ei                            : Energy grid point             !
!  ***************************** OUTPUT ******************************  !
!  TYPE(green) Gr_nn(neph)%G(unitdimensions,unitdimensions) :           !
!                                                  [complex] G^r_{n,n}  !
!  TYPE(green) Gr_1n(neph)%G(NL,unitdimensions) :  [complex] G^r_{1,n}  !
!  TYPE(green) Gr_Mn(neph)%G(NR,unitdimensions) :  [complex] G^r_{M,n}  !
!  complex(8) Gr_1M(NL,NR)                      : G^r_{1,M}             !
!  *******************************************************************  !
  subroutine GFfull (Ei, ispin)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_options,  only: ntypeunits, nunits
    use idsrdr_units,    only: unit_type, unitdimensions, unitshift,    &
                               S1unit, H1unit, Sunits, Hunits
    use idsrdr_leads,    only: NL, NR, Sigma_L, Sigma_R
    use idsrdr_ephcoupl, only: ephIdx, norbDyn, idxF, idxL
    use idsrdr_check,    only: CHECKzsytrf, CHECKzsytri

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    integer :: I, n, utype, dim, dimbfr, idx, ueph
    integer, allocatable, dimension (:) :: ipiv
    complex(8), allocatable, dimension (:,:) :: V ! pristine coupling
    complex(8), allocatable, dimension (:,:) :: aux1, aux2, aux3
    complex(8), allocatable, dimension (:,:) :: foo1L, foo2L,           &
                                                foo1R, foo2R, foo3
    external :: zsymm, zgemm

    if (IOnode) write (6,'(a)', advance='no')                           &
            '      computing full Greens functions... '

!   Initialize variables.
    utype = ntypeunits + 1 ! current unit type (first unit)
    dim = unitdimensions(utype) ! current type dimension
    n = unitdimensions(ntypeunits) ! pristine dimension
    idx = ephIdx(utype) ! e-ph unit type
    ueph = 1 ! eph units indexing

!   Allocate matrices.
    allocate (V(n,n))
    allocate (aux1(n,n))
    allocate (aux2(n,n))
    allocate (foo1L(NL,n))
    allocate (foo2L(NL,n))
    allocate (foo1R(NR,n))
    allocate (foo2R(NR,n))

!   Coupling matrix (pristine).
    V = (Ei-unitshift(ntypeunits))*S1unit - H1unit(:,:,ispin)

    if (idx /= 0) then

!      Allocate auxiliary matrix.
       allocate (aux3(dim,dim))

!      ('aux3 = E*S - H - Sigma_L')
       aux3 = (0.d0,0.d0)
       aux3 = (Ei-unitshift(utype))*Sunits(utype)%S                     &
              - Hunits(utype)%H(:,:,ispin)
       aux3(1:NL,1:NL) = aux3(1:NL,1:NL) - Sigma_L

!      ('aux2 = V*GR_pp')
       aux1 = GR_pp(ueph)%G(1:n,1:n)
       call zsymm ('R', 'L', n, n, (1.d0,0.d0), aux1, n,                &
                   V, n, (0.d0,0.d0), aux2, n)

!      ('aux1 = aux2*V^dagger')
       call zgemm ('N', 'C', n, n, n, (1.d0,0.d0), aux2, n,             &
                   V, n, (0.d0,0.d0), aux1, n)

!      ('aux3 = aux3 - aux1')
       aux3(dim-n+1:dim,dim-n+1:dim) = aux3(dim-n+1:dim,dim-n+1:dim)    &
                                       - aux1

!      ('aux3 = aux3^-1')
       allocate (ipiv(dim))
       call CHECKzsytrf (dim, 'L', aux3, ipiv)
       call CHECKzsytri (dim, 'L', aux3, ipiv)
       deallocate (ipiv)

!      ('Gr_nn = aux3')
       Gr_nn(ueph)%G = aux3(idxF(idx):idxL(idx),idxF(idx):idxL(idx))

!      ('Gr_1n = aux3')
       Gr_1n(ueph)%G = aux3(1:NL,idxF(idx):idxL(idx))

!      Allocate auxiliary matrix.
       allocate (foo3(n,norbDyn(idx)))

!      ('foo2R = GR_Mp*V^dagger')
       foo1R = GR_Mp(ueph)%G(1:NR,1:n)
       call zgemm ('N', 'C', NR, n, n, (1.d0,0.d0), foo1R, NR,          &
                   V, n, (0.d0,0.d0), foo2R, NR)

!      ('Gr_Mn = - foo2R * Gr_nn')
       foo3 = aux3(dim-n+1:dim,idxF(idx):idxL(idx))
       call zgemm ('N', 'N', NR, norbDyn(idx), n, (-1.d0,0.d0),         &
                   foo2R, NR, foo3, n, (0.d0,0.d0), Gr_Mn(ueph)%G, NR)

!      Free memory.
       deallocate (foo3)
       deallocate (aux3)

       ueph = ueph + 1

    endif

!   Loop over the blocks from left to right.
    do I = 2,nunits+1

!      Assign auxiliary variables.
       dimbfr = dim ! last unit dimension
       utype = unit_type(I) ! current unit type
       dim = unitdimensions(utype) ! current type dimension
       idx = ephIdx(utype) ! e-ph unit type

       if (idx /= 0) then

!         Allocate auxiliary matrix.
          allocate (aux3(dim,dim))

!         ('aux3 = E*S - H')
          aux3 = (Ei-unitshift(utype))*Sunits(utype)%S                  &
                 - Hunits(utype)%H(:,:,ispin)

!         ('aux2 = GL_mm*V')
          aux1 = GL_mm(ueph)%G(dimbfr-n+1:dimbfr,dimbfr-n+1:dimbfr)
          call zsymm ('L', 'L', n, n, (1.d0,0.d0), aux1, n,             &
                      V, n, (0.d0,0.d0), aux2, n)

!         ('aux1 = V^dagger*aux2')
          call zgemm ('C', 'N', n, n, n, (1.d0,0.d0), V, n,             &
                      aux2, n, (0.d0,0.d0), aux1, n)

!         ('aux3 = aux3 - aux1')
          aux3(1:n,1:n) = aux3(1:n,1:n) - aux1

!         ('aux2 = V*GR_pp')
          aux1 = GR_pp(ueph)%G(1:n,1:n)
          call zsymm ('R', 'L', n, n, (1.d0,0.d0), aux1, n,             &
                      V, n, (0.d0,0.d0), aux2, n)

!         ('aux1 = aux2*V^dagger')
          call zgemm ('N', 'C', n, n, n, (1.d0,0.d0), aux2, n,          &
                      V, n, (0.d0,0.d0), aux1, n)

!         ('aux3 = aux3 - aux1')
          aux3(dim-n+1:dim,dim-n+1:dim) = aux3(dim-n+1:dim,dim-n+1:dim) &
                                          - aux1

!         ('aux3 = aux3^-1')
          allocate (ipiv(dim))
          call CHECKzsytrf (dim, 'L', aux3, ipiv)
          call CHECKzsytri (dim, 'L', aux3, ipiv)
          deallocate (ipiv)

!         ('Gr_nn = aux3')
          Gr_nn(ueph)%G = aux3(idxF(idx):idxL(idx),idxF(idx):idxL(idx))

!         Allocate auxiliary matrix.
          allocate (foo3(n,norbDyn(idx)))

!         ('foo2L = GL_1m*V')
          foo1L = GL_1m(ueph)%G(1:NL,dimbfr-n+1:dimbfr)
          call zgemm ('N', 'N', NL, n, n, (1.d0,0.d0), foo1L, NL,       &
                      V, n, (0.d0,0.d0), foo2L, NL)

!         ('Gr_1n = - foo2L * Gr_nn')
          foo3 = aux3(1:n,idxF(idx):idxL(idx))
          call zgemm ('N', 'N', NL, norbDyn(idx), n, (-1.d0,0.d0),      &
                      foo2L, NL, foo3, n, (0.d0,0.d0), Gr_1n(ueph)%G, NL)

!         ('foo2R = GR_Mp*V^dagger')
          foo1R = GR_Mp(ueph)%G(1:NR,1:n)
          call zgemm ('N', 'C', NR, n, n, (1.d0,0.d0), foo1R, NR,       &
                      V, n, (0.d0,0.d0), foo2R, NR)

!         ('Gr_Mn = - foo2R * Gr_nn')
          foo3 = aux3(dim-n+1:dim,idxF(idx):idxL(idx))
          call zgemm ('N', 'N', NR, norbDyn(idx), n, (-1.d0,0.d0),      &
                      foo2R, NR, foo3, n, (0.d0,0.d0), Gr_Mn(ueph)%G, NR)

!         Free memory.
          deallocate (foo3)
          deallocate (aux3)

          ueph = ueph + 1

       endif
    enddo

!   Assign auxiliary variables.
    dimbfr = dim ! last unit dimension
    utype = ntypeunits + 2 ! current unit type (last unit)
    dim = unitdimensions(utype) ! current type dimension
    idx = ephIdx(utype) ! e-ph unit type

    if (idx /= 0) then

!      Allocate auxiliary matrix.
       allocate (aux3(dim,dim))

!      ('aux3 = E*S - H')
       aux3 = (Ei-unitshift(utype))*Sunits(utype)%S                     &
              - Hunits(utype)%H(:,:,ispin)
       aux3(dim-NR+1:dim,dim-NR+1:dim) =                                &
                                        aux3(dim-NR+1:dim,dim-NR+1:dim) &
                                        - Sigma_R

!      ('aux2 = GL_mm*V')
       aux1 = GL_mm(ueph)%G(dimbfr-n+1:dimbfr,dimbfr-n+1:dimbfr)
       call zsymm ('L', 'L', n, n, (1.d0,0.d0), aux1, n,                &
                   V, n, (0.d0,0.d0), aux2, n)

!      ('aux1 = V^dagger*aux2')
       call zgemm ('C', 'N', n, n, n, (1.d0,0.d0), V, n,                &
                   aux2, n, (0.d0,0.d0), aux1, n)

!      ('aux3 = aux3 - aux1')
       aux3(1:n,1:n) = aux3(1:n,1:n) - aux1

!      ('aux3 = aux3^-1')
       allocate (ipiv(dim))
       call CHECKzsytrf (dim, 'L', aux3, ipiv)
       call CHECKzsytri (dim, 'L', aux3, ipiv)
       deallocate (ipiv)

!      ('Gr_nn = aux3')
       Gr_nn(ueph)%G = aux3(idxF(idx):idxL(idx),idxF(idx):idxL(idx))

!      ('Gr_Mn = aux3')
       Gr_Mn(ueph)%G = aux3(dim-NR+1:dim,idxF(idx):idxL(idx))

!      Allocate auxiliary matrix.
       allocate (foo3(n,dim))
       deallocate (aux1)
       allocate (aux1(NL,dim))
 
!      ('foo2L = GL_1m*V')
       foo1L = GL_1m(ueph)%G(1:NL,dimbfr-n+1:dimbfr)
       call zgemm ('N', 'N', NL, n, n, (1.d0,0.d0), foo1L, NL,          &
                   V, n, (0.d0,0.d0), foo2L, NL)

!      ('Gr_1n = - foo2L * Gr_nn')
       foo3 = aux3(1:n,1:dim)
       call zgemm ('N', 'N', NL, dim, n, (-1.d0,0.d0), foo2L, NL,       &
                   foo3, n, (0.d0,0.d0), aux1, NL)
       Gr_1n(ueph)%G = aux1(1:NL,idxF(idx):idxL(idx))
       Gr_1M = aux1(1:NL,dim-NR+1:dim)

!      Free memory.
       deallocate (foo3)
       deallocate (aux3)

    else

!      Allocate auxiliary matrix.
       allocate (aux3(dim,dim))

!      ('aux3 = E*S - H')
       aux3 = (Ei-unitshift(utype))*Sunits(utype)%S                     &
              - Hunits(utype)%H(:,:,ispin)
       aux3(dim-NR+1:dim,dim-NR+1:dim) =                                &
                                        aux3(dim-NR+1:dim,dim-NR+1:dim) &
                                        - Sigma_R

!      ('aux2 = GL_nn*V')
       aux1 = GL_nn(dimbfr-n+1:dimbfr,dimbfr-n+1:dimbfr)
       call zsymm ('L', 'L', n, n, (1.d0,0.d0), aux1, n,                &
                   V, n, (0.d0,0.d0), aux2, n)

!      ('aux1 = V^dagger*aux2')
       call zgemm ('C', 'N', n, n, n, (1.d0,0.d0), V, n,                &
                   aux2, n, (0.d0,0.d0), aux1, n)

!      ('aux3 = aux3 - aux1')
       aux3(1:n,1:n) = aux3(1:n,1:n) - aux1

!      ('aux3 = aux3^-1')
       allocate (ipiv(dim))
       call CHECKzsytrf (dim, 'L', aux3, ipiv)
       call CHECKzsytri (dim, 'L', aux3, ipiv)
       deallocate (ipiv)

!      Allocate auxiliary matrix.
       allocate (foo3(n,NR))

!      ('foo2L = GL_1m*V')
       foo1L = GL_1N(1:NL,dimbfr-n+1:dimbfr)
       call zgemm ('N', 'N', NL, n, n, (1.d0,0.d0), foo1L, NL,          &
                   V, n, (0.d0,0.d0), foo2L, NL)

!      ('Gr_1n = - foo2L * G_nn')
       foo3 = aux3(1:n,dim-NR+1:dim)
       call zgemm ('N', 'N', NL, NR, n, (-1.d0,0.d0), foo2L, NL,        &
                   foo3, NR, (0.d0,0.d0), Gr_1M, NL)

!      Free memory.
       deallocate (foo3)
       deallocate (aux3)

    endif

!   Free memory.
    deallocate (V)
    deallocate (aux1)
    deallocate (aux2)
    deallocate (foo1L)
    deallocate (foo2L)
    deallocate (foo1R)
    deallocate (foo2R)

    if (IOnode) write(6,'(a)') " ok!"


  end subroutine GFfull


!  *******************************************************************  !
!                              LRsweepDisk                              !
!  *******************************************************************  !
!  Description: left-to-right sweep for computing 'G^L_{n-1,n-1}' and   !
!  'G^L_{1,n-1}'.                                                       !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode                       : True if it is the I/O node    !
!  integer ntypeunits                   : Number of unit types          !
!  integer nunits                       : Total number of units         !
!  integer unit_type(nunits+2)          : Units types                   !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  real*8 unitshift(ntypeunits+2)       : Units shift                   !
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
!  integer ephIndic(ntypeunits+2)        : E-ph interaction in          !
!                                          unit type (0 or 1)           !
!  integer NL                           : Number of left lead orbitals  !
!  complex(8) Sigma_L(NL,NL)            : Left-lead self-energy         !
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 Ei                            : Energy grid point             !
!  ***************************** OUTPUT ******************************  !
!  TYPE(green) GL_mm(neph)%G(unitdimensions,unitdimensions) :           !
!                                              [complex] G^L_{n-1,n-1}  !
!  TYPE(green) GL_1m(neph)%G(NL,unitdimensions) :                       !
!                                              [complex] G^L_{1,n-1}    !
!  *******************************************************************  !
  subroutine LRsweepDisk (Ei, ispin)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_options,  only: ntypeunits, nunits
    use idsrdr_units,    only: unit_type, unitdimensions, unitshift,    &
                               S1unit, H1unit, Sunits, Hunits, ephIndic
    use idsrdr_leads,    only: NL, Sigma_L
    use idsrdr_check,    only: CHECKzsytrf, CHECKzsytri
    use idsrdr_io,       only: IOopenStream, IOcloseStream

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    integer :: I, n, utype, dim, dimbfr
    integer, allocatable, dimension (:) :: ipiv
    complex(8), allocatable, dimension (:,:) :: V ! pristine coupling
    complex(8), allocatable, dimension (:,:) :: Gbfr ! GF before
    complex(8), allocatable, dimension (:,:) :: Gaft ! GF after
    complex(8), allocatable, dimension (:,:) :: Gbfr_1m ! GF_1m before
    complex(8), allocatable, dimension (:,:) :: Gaft_1m ! GF_1m after
    complex(8), allocatable, dimension (:,:) :: aux1, aux2
    complex(8), allocatable, dimension (:,:) :: foo1, foo2, foo3
    external :: zsymm, zgemm

    if (IOnode) write (6,'(a)', advance='no')                           &
            '      computing left-to-right sweep... '

!   Open Green's functions files.
    call IOopenStream (GL_mm_disk%fname, GL_mm_disk%lun)
    call IOopenStream (GL_1m_disk%fname, GL_1m_disk%lun)

!   Initialize variables.
    utype = ntypeunits + 1 ! current unit type (first unit)
    dim = unitdimensions(utype) ! current type dimension
    n = unitdimensions(ntypeunits) ! pristine dimension

!   Allocate matrices.
    allocate (V(n,n))
    allocate (Gbfr(dim,dim))
    allocate (Gbfr_1m(dim,dim))
    allocate (ipiv(dim))
    allocate (aux1(n,n))
    allocate (aux2(n,n))
    allocate (foo1(NL,n))
    allocate (foo2(NL,n))

!   Coupling matrix (pristine).
    V = (Ei-unitshift(ntypeunits))*S1unit - H1unit(:,:,ispin)

!   Compute the first step.
    Gbfr = (Ei-unitshift(utype))*Sunits(utype)%S                        &
           - Hunits(utype)%H(:,:,ispin)
    Gbfr(1:NL,1:NL) = Gbfr(1:NL,1:NL) - Sigma_L
    call CHECKzsytrf (dim, 'L', Gbfr, ipiv)
    call CHECKzsytri (dim, 'L', Gbfr, ipiv)
    Gbfr_1m = Gbfr

!   Free memory.
    deallocate (ipiv)

!   Loop over the blocks from left to right.
    do I = 2,nunits+1

!      Assign auxiliary variables.
       dimbfr = dim ! old dimension
       utype = unit_type(I) ! current unit type
       dim = unitdimensions(utype) ! current type dimension

!      Store matrix if required.
       if (ephIndic(utype)) then
          write (GL_mm_disk%lun) Gbfr
          write (GL_1m_disk%lun) Gbfr_1m
       endif

!      ('aux2 = Gbfr*V')
       aux1 = Gbfr(dimbfr-n+1:dimbfr,dimbfr-n+1:dimbfr)
       call zsymm ('L', 'L', n, n, (1.d0,0.d0), aux1, n,                &
                   V, n, (0.d0,0.d0), aux2, n)

!      ('aux1 = V^dagger*aux2')
       call zgemm ('C', 'N', n, n, n, (1.d0,0.d0), V, n,                &
                   aux2, n, (0.d0,0.d0), aux1, n)

!      (Re-)Allocate matrices and array.
       allocate (Gaft(dim,dim))
       allocate (ipiv(dim))

!      ('Gaft = (E*S - H - V^T*Gbfr*V)^-1')
       Gaft = (Ei-unitshift(utype))*Sunits(utype)%S                     &
              - Hunits(utype)%H(:,:,ispin)
       Gaft(1:n,1:n) = Gaft(1:n,1:n) - aux1
       call CHECKzsytrf (dim, 'L', Gaft, ipiv)
       call CHECKzsytri (dim, 'L', Gaft, ipiv)

!      (Re-)Allocate matrices.
       allocate (Gaft_1m(NL,dim))
       allocate (foo3(n,dim))

!      ('foo2 = Gbfr_1m*V')
       foo1 = Gbfr_1m(1:NL,dimbfr-n+1:dimbfr)
       call zgemm ('N', 'N', NL, n, n, (1.d0,0.d0), foo1, NL,           &
                   V, n, (0.d0,0.d0), foo2, NL)

!      ('Gaft_1m = - foo2 * Gaft')
       foo3 = Gaft(1:n,:)
       call zgemm ('N', 'N', NL, dim, n, (-1.d0,0.d0), foo2, NL,        &
                   foo3, n, (0.d0,0.d0), Gaft_1m, NL)

!      Reallocate matrix.
       deallocate (Gbfr)
       deallocate (Gbfr_1m)
       allocate (Gbfr(dim,dim))
       allocate (Gbfr_1m(NL,dim))

       Gbfr = Gaft
       Gbfr_1m = Gaft_1m

!      (Re-)Free memory.
       deallocate (Gaft)
       deallocate (Gaft_1m)
       deallocate (foo3)
       deallocate (ipiv)

    enddo

!   Store matrix if required.
    if (ephIndic(ntypeunits+2)) then
       write (GL_mm_disk%lun) Gbfr
       write (GL_1m_disk%lun) Gbfr_1m
    else
       GL_nn = Gbfr
       GL_1N = Gbfr_1m
    endif

!   Close Green's functions files.
    call IOcloseStream (GL_mm_disk%fname, GL_mm_disk%lun)
    call IOcloseStream (GL_1m_disk%fname, GL_1m_disk%lun)

!   Free memory.
    deallocate (V)
    deallocate (Gbfr)
    deallocate (Gbfr_1m)
    deallocate (aux1)
    deallocate (aux2)
    deallocate (foo1)
    deallocate (foo2)

    if (IOnode) write(6,'(a)') " ok!"


  end subroutine LRsweepDisk


!  *******************************************************************  !
!                              RLsweepDisk                              !
!  *******************************************************************  !
!  Description: right-to-left sweep for computing 'G^R_{n+1,n+1}' and   !
!  'G^R_{M,n+1}'.                                                       !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode                       : True if it is the I/O node    !
!  integer ntypeunits                   : Number of unit types          !
!  integer nunits                       : Total number of units         !
!  integer unit_type(nunits+2)          : Units types                   !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  real*8 unitshift(ntypeunits+2)       : Units shift                   !
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
!  integer ephIndic(ntypeunits+2)        : E-ph interaction in          !
!                                          unit type (0 or 1)           !
!  integer nunitseph                   : Number of units with eph       !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Sigma_R(NR,NR)           : Right-lead self-energy         !
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 Ei                            : Energy grid point             !
!  ***************************** OUTPUT ******************************  !
!  TYPE(green) GR_pp(neph)%G(unitdimensions,unitdimensions) :           !
!                                              [complex] G^R_{n+1,n+1}  !
!  TYPE(green) GR_Mp(neph)%G(NR,unitdimensions) :                       !
!                                              [complex] G^R_{M,n+1}    !
!  *******************************************************************  !
  subroutine RLsweepDisk (Ei, ispin)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_options,  only: ntypeunits, nunits
    use idsrdr_units,    only: unit_type, unitdimensions, unitshift,    &
                               S1unit, H1unit, Sunits, Hunits,          &
                               ephIndic, nunitseph
    use idsrdr_leads,    only: NR, Sigma_R
    use idsrdr_check,    only: CHECKzsytrf, CHECKzsytri
    use idsrdr_io,       only: IOopenStream, IOcloseStream

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    integer :: I, n, utype, dim, ueph, nbytes_pp, nbytes_Mp
    integer, allocatable, dimension (:) :: ipiv
    complex(8), allocatable, dimension (:,:) :: V ! pristine coupling
    complex(8), allocatable, dimension (:,:) :: Gbfr ! GF before
    complex(8), allocatable, dimension (:,:) :: Gaft ! GF after
    complex(8), allocatable, dimension (:,:) :: Gbfr_Mp ! GF_Mp before
    complex(8), allocatable, dimension (:,:) :: Gaft_Mp ! GF_Mp after
    complex(8), allocatable, dimension (:,:) :: aux1, aux2
    complex(8), allocatable, dimension (:,:) :: foo1, foo2, foo3
    external :: zsymm, zgemm

    if (IOnode) write (6,'(a)', advance='no')                           &
            '      computing right-to-left sweep... '

!   Open Green's functions files.
    call IOopenStream (GR_pp_disk%fname, GR_pp_disk%lun)
    call IOopenStream (GR_Mp_disk%fname, GR_Mp_disk%lun)

!   Initialize variables.
    utype = ntypeunits + 2 ! current unit type (last unit)
    dim = unitdimensions(utype) ! current type dimension
    n = unitdimensions(ntypeunits) ! pristine dimension
    if (ephIndic(utype)) then ! eph units indexing
       ueph = nunitseph - 1
    else
       ueph = nunitseph
    endif
    nbytes_pp = 1 ! position at 'GR_pp' file
    nbytes_Mp = 1 ! position at 'GR_Mp' file

!   Allocate matrices.
    allocate (V(n,n))
    allocate (Gbfr(dim,dim))
    allocate (Gbfr_Mp(dim,dim))
    allocate (ipiv(dim))
    allocate (aux1(n,n))
    allocate (aux2(n,n))
    allocate (foo1(NR,n))
    allocate (foo2(NR,n))

!   Coupling matrix (pristine).
    V = (Ei-unitshift(ntypeunits))*S1unit - H1unit(:,:,ispin)

!   Compute the first step.
    Gbfr = (Ei-unitshift(utype))*Sunits(utype)%S                        &
           - Hunits(utype)%H(:,:,ispin)
    Gbfr(dim-NR+1:dim,dim-NR+1:dim) = Gbfr(dim-NR+1:dim,dim-NR+1:dim)   &
                                      - Sigma_R
    call CHECKzsytrf (dim, 'L', Gbfr, ipiv)
    call CHECKzsytri (dim, 'L', Gbfr, ipiv)
    Gbfr_Mp = Gbfr

!   Free memory.
    deallocate (ipiv)

!   Loop over the blocks from right to left.
    do I = nunits+1,2,-1

!      Assign auxiliary variables.
       utype = unit_type(I) ! current unit type
       dim = unitdimensions(utype) ! current type dimension

!      Store matrix if required.
       if (ephIndic(utype)) then
          write (GR_pp_disk%lun) Gbfr
          pos_pp(ueph) = nbytes_pp
          nbytes_pp = nbytes_pp + sizeof(Gbfr)
          write (GR_Mp_disk%lun) Gbfr_Mp
          pos_Mp(ueph) = nbytes_Mp
          nbytes_Mp = nbytes_Mp + sizeof(Gbfr_Mp)
          ueph = ueph - 1
       endif

!      ('aux2 = V*Gbfr')
       aux1 = Gbfr(1:n,1:n)
       call zsymm ('R', 'L', n, n, (1.d0,0.d0), aux1, n,                &
                   V, n, (0.d0,0.d0), aux2, n)

!      ('aux1 = aux2*V^dagger')
       call zgemm ('N', 'C', n, n, n, (1.d0,0.d0), aux2, n,             &
                   V, n, (0.d0,0.d0), aux1, n)

!      (Re-)Allocate matrices and array.
       allocate (Gaft(dim,dim))
       allocate (ipiv(dim))

!      ('Gaft = (E*S - H - V^T*Gbfr*V)^-1')
       Gaft = (Ei-unitshift(utype))*Sunits(utype)%S                     &
              - Hunits(utype)%H(:,:,ispin)
       Gaft(dim-n+1:dim,dim-n+1:dim) = Gaft(dim-n+1:dim,dim-n+1:dim)    &
                                       - aux1
       call CHECKzsytrf (dim, 'L', Gaft, ipiv)
       call CHECKzsytri (dim, 'L', Gaft, ipiv)

!      (Re-)Allocate matrices.
       allocate (Gaft_Mp(NR,dim))
       allocate (foo3(n,dim))

!      ('foo2 = Gbfr_Mp*V^dagger')
       foo1 = Gbfr_Mp(1:NR,1:n)
       call zgemm ('N', 'C', NR, n, n, (1.d0,0.d0), foo1, NR,           &
                   V, n, (0.d0,0.d0), foo2, NR)

!      ('Gaft_Mp = - foo2 * Gaft')
       foo3 = Gaft(dim-n+1:dim,:)
       call zgemm ('N', 'N', NR, dim, n, (-1.d0,0.d0), foo2, NR,        &
                   foo3, n, (0.d0,0.d0), Gaft_Mp, NR)

!      Reallocate matrix.
       deallocate (Gbfr)
       deallocate (Gbfr_Mp)
       allocate (Gbfr(dim,dim))
       allocate (Gbfr_Mp(NR,dim))

       Gbfr = Gaft
       Gbfr_Mp = Gaft_Mp

!      (Re-)Free memory.
       deallocate (Gaft)
       deallocate (Gaft_Mp)
       deallocate (foo3)
       deallocate (ipiv)

    enddo

!   Store matrix if required.
    if (ephIndic(ntypeunits+1)) then
       write (GR_pp_disk%lun) Gbfr
       pos_pp(ueph) = nbytes_pp
       write (GR_Mp_disk%lun) Gbfr_Mp
       pos_Mp(ueph) = nbytes_Mp
    endif

!   Close Green's functions files.
    call IOcloseStream (GR_pp_disk%fname, GR_pp_disk%lun)
    call IOcloseStream (GR_Mp_disk%fname, GR_Mp_disk%lun)

!   Free memory.
    deallocate (V)
    deallocate (Gbfr)
    deallocate (Gbfr_Mp)
    deallocate (aux1)
    deallocate (aux2)
    deallocate (foo1)
    deallocate (foo2)

    if (IOnode) write(6,'(a)') " ok!"


  end subroutine RLsweepDisk


!  *******************************************************************  !
!                              GFfullDisk                               !
!  *******************************************************************  !
!  Description: computes the full Green's function of required units.   !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode                       : True if it is the I/O node    !
!  integer ntypeunits                   : Number of unit types          !
!  integer nunits                       : Total number of units         !
!  integer unit_type(nunits+2)          : Units types                   !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  real*8 unitshift(ntypeunits+2)       : Units shift                   !
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
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Sigma_L(NL,NL)           : Left-lead self-energy          !
!  complex(8) Sigma_R(NR,NR)           : Right-lead self-energy         !
!  integer ephIdx(ntypeunits+2)        : Unit index (those with e-ph)   !
!  integer norbDyn(neph)               : Number of orbitals from        !
!                                        dynamic atoms                  !
!  integer idxF(neph)                  : First dynamic atom orbital     !
!  integer idxL(neph)                  : Last dynamic atom orbital      !
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 Ei                            : Energy grid point             !
!  ***************************** OUTPUT ******************************  !
!  TYPE(green) Gr_nn(neph)%G(unitdimensions,unitdimensions) :           !
!                                                  [complex] G^r_{n,n}  !
!  TYPE(green) Gr_1n(neph)%G(NL,unitdimensions) :  [complex] G^r_{1,n}  !
!  TYPE(green) Gr_Mn(neph)%G(NR,unitdimensions) :  [complex] G^r_{M,n}  !
!  complex(8) Gr_1M(NL,NR)                      : G^r_{1,M}             !
!  *******************************************************************  !
  subroutine GFfullDisk (Ei, ispin)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_options,  only: ntypeunits, nunits
    use idsrdr_units,    only: unit_type, unitdimensions, unitshift,    &
                               S1unit, H1unit, Sunits, Hunits
    use idsrdr_leads,    only: NL, NR, Sigma_L, Sigma_R
    use idsrdr_ephcoupl, only: ephIdx, norbDyn, idxF, idxL
    use idsrdr_check,    only: CHECKzsytrf, CHECKzsytri
    use idsrdr_io,       only: IOopenStream, IOcloseStream

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    integer :: I, n, utype, dim, dimbfr, dimaft, idx, ueph
    integer, allocatable, dimension (:) :: ipiv
    complex(8), allocatable, dimension (:,:) :: V ! pristine coupling
    complex(8), allocatable, dimension (:,:) :: aux1, aux2, aux3, aux4
    complex(8), allocatable, dimension (:,:) :: foo1L, foo2L,           &
                                                foo1R, foo2R, foo3
    external :: zsymm, zgemm

    if (IOnode) write (6,'(a)', advance='no')                           &
            '      computing full Greens functions... '

!   Open Green's functions files.
    call IOopenStream (GL_mm_disk%fname, GL_mm_disk%lun)
    call IOopenStream (GL_1m_disk%fname, GL_1m_disk%lun)
    call IOopenStream (GR_pp_disk%fname, GR_pp_disk%lun)
    call IOopenStream (GR_Mp_disk%fname, GR_Mp_disk%lun)
    call IOopenStream (Gr_nn_disk%fname, Gr_nn_disk%lun)
    call IOopenStream (Gr_1n_disk%fname, Gr_1n_disk%lun)
    call IOopenStream (Gr_Mn_disk%fname, Gr_Mn_disk%lun)

!   Initialize variables.
    utype = ntypeunits + 1 ! current unit type (first unit)
    dim = unitdimensions(utype) ! current type dimension
    dimaft = unitdimensions(unit_type(2)) ! next unit dimension
    n = unitdimensions(ntypeunits) ! pristine dimension
    idx = ephIdx(utype) ! e-ph unit type
    ueph = 1 ! eph units indexing

!   Allocate matrices.
    allocate (V(n,n))
    allocate (aux1(n,n))
    allocate (aux2(n,n))
    allocate (foo1L(NL,n))
    allocate (foo2L(NL,n))
    allocate (foo1R(NR,n))
    allocate (foo2R(NR,n))

!   Coupling matrix (pristine).
    V = (Ei-unitshift(ntypeunits))*S1unit - H1unit(:,:,ispin)

    if (idx /= 0) then

!      Allocate auxiliary matrix.
       allocate (aux3(dim,dim))

!      ('aux3 = E*S - H - Sigma_L')
       aux3 = (0.d0,0.d0)
       aux3 = (Ei-unitshift(utype))*Sunits(utype)%S                     &
              - Hunits(utype)%H(:,:,ispin)
       aux3(1:NL,1:NL) = aux3(1:NL,1:NL) - Sigma_L

!      ('aux2 = V*GR_pp')
       call greenloadR (GR_pp_disk%lun, dimaft, dimaft, 1, n,           &
                        1, n, aux1, n, n, pos_pp(ueph))
       call zsymm ('R', 'L', n, n, (1.d0,0.d0), aux1, n,                &
                   V, n, (0.d0,0.d0), aux2, n)

!      ('aux1 = aux2*V^dagger')
       call zgemm ('N', 'C', n, n, n, (1.d0,0.d0), aux2, n,             &
                   V, n, (0.d0,0.d0), aux1, n)

!      ('aux3 = aux3 - aux1')
       aux3(dim-n+1:dim,dim-n+1:dim) = aux3(dim-n+1:dim,dim-n+1:dim)    &
                                       - aux1

!      ('aux3 = aux3^-1')
       allocate (ipiv(dim))
       call CHECKzsytrf (dim, 'L', aux3, ipiv)
       call CHECKzsytri (dim, 'L', aux3, ipiv)
       deallocate (ipiv)

!      ('Gr_nn = aux3')
       allocate (aux4(norbDyn(idx),norbDyn(idx)))
       aux4 = aux3(idxF(idx):idxL(idx),idxF(idx):idxL(idx))
       write (Gr_nn_disk%lun) aux4
       deallocate (aux4)

!      ('Gr_1n = aux3')
       allocate (aux4(NL,norbDyn(idx)))
       aux4 = aux3(1:NL,idxF(idx):idxL(idx))
       write (Gr_1n_disk%lun) aux4
       deallocate (aux4)

!      Allocate auxiliary matrix.
       allocate (foo3(n,norbDyn(idx)))

!      ('foo2R = GR_Mp*V^dagger')
       call greenloadR (GR_Mp_disk%lun, NR, dimaft, 1, NR,              &
                        1, n, foo1R, NR, n, pos_Mp(ueph))
       call zgemm ('N', 'C', NR, n, n, (1.d0,0.d0), foo1R, NR,          &
                   V, n, (0.d0,0.d0), foo2R, NR)

!      ('Gr_Mn = - foo2R * Gr_nn')
       allocate (aux4(NR,norbDyn(idx)))
       foo3 = aux3(dim-n+1:dim,idxF(idx):idxL(idx))
       call zgemm ('N', 'N', NR, norbDyn(idx), n, (-1.d0,0.d0),         &
                   foo2R, NR, foo3, n, (0.d0,0.d0), aux4, NR)
       write (Gr_Mn_disk%lun) aux4
       deallocate (aux4)

!      Free memory.
       deallocate (foo3)
       deallocate (aux3)

       ueph = ueph + 1

    endif

!   Loop over the blocks from left to right.
    do I = 2,nunits

!      Assign auxiliary variables.
       dimbfr = dim ! last unit dimension
       utype = unit_type(I) ! current unit type
       dim = unitdimensions(utype) ! current type dimension
       dimaft = unitdimensions(unit_type(I+1)) ! next unit dimension
       idx = ephIdx(utype) ! e-ph unit type

       if (idx /= 0) then

!         Allocate auxiliary matrix.
          allocate (aux3(dim,dim))

!         ('aux3 = E*S - H')
          aux3 = (Ei-unitshift(utype))*Sunits(utype)%S                  &
                 - Hunits(utype)%H(:,:,ispin)

!         ('aux2 = GL_mm*V')
          call greenload (GL_mm_disk%lun, dimbfr, dimbfr, dimbfr-n+1,   &
                          dimbfr, dimbfr-n+1, dimbfr, aux1, n, n)
          call zsymm ('L', 'L', n, n, (1.d0,0.d0), aux1, n,             &
                      V, n, (0.d0,0.d0), aux2, n)

!         ('aux1 = V^dagger*aux2')
          call zgemm ('C', 'N', n, n, n, (1.d0,0.d0), V, n,             &
                      aux2, n, (0.d0,0.d0), aux1, n)

!         ('aux3 = aux3 - aux1')
          aux3(1:n,1:n) = aux3(1:n,1:n) - aux1

!         ('aux2 = V*GR_pp')
          call greenloadR (GR_pp_disk%lun, dimaft, dimaft, 1, n,        &
                           1, n, aux1, n, n, pos_pp(ueph))
          call zsymm ('R', 'L', n, n, (1.d0,0.d0), aux1, n,             &
                      V, n, (0.d0,0.d0), aux2, n)

!         ('aux1 = aux2*V^dagger')
          call zgemm ('N', 'C', n, n, n, (1.d0,0.d0), aux2, n,          &
                      V, n, (0.d0,0.d0), aux1, n)

!         ('aux3 = aux3 - aux1')
          aux3(dim-n+1:dim,dim-n+1:dim) = aux3(dim-n+1:dim,dim-n+1:dim) &
                                          - aux1

!         ('aux3 = aux3^-1')
          allocate (ipiv(dim))
          call CHECKzsytrf (dim, 'L', aux3, ipiv)
          call CHECKzsytri (dim, 'L', aux3, ipiv)
          deallocate (ipiv)

!         ('Gr_nn = aux3')
          allocate (aux4(norbDyn(idx),norbDyn(idx)))
          aux4 = aux3(idxF(idx):idxL(idx),idxF(idx):idxL(idx))
          write (Gr_nn_disk%lun) aux4
          deallocate (aux4)

!         Allocate auxiliary matrix.
          allocate (foo3(n,norbDyn(idx)))

!         ('foo2L = GL_1m*V')
          call greenload (GL_1m_disk%lun, NL, dimbfr, 1, NL,            &
                          dimbfr-n+1, dimbfr, foo1L, NL, n)
          call zgemm ('N', 'N', NL, n, n, (1.d0,0.d0), foo1L, NL,       &
                      V, n, (0.d0,0.d0), foo2L, NL)

!         ('Gr_1n = - foo2L * Gr_nn')
          allocate (aux4(NL,norbDyn(idx)))
          foo3 = aux3(1:n,idxF(idx):idxL(idx))
          call zgemm ('N', 'N', NL, norbDyn(idx), n, (-1.d0,0.d0),      &
                      foo2L, NL, foo3, n, (0.d0,0.d0), aux4, NL)
          write (Gr_1n_disk%lun) aux4
          deallocate (aux4)

!         ('foo2R = GR_Mp*V^dagger')
          call greenloadR (GR_Mp_disk%lun, NR, dimaft, 1, NR,           &
                           1, n, foo1R, NR, n, pos_Mp(ueph))
          call zgemm ('N', 'C', NR, n, n, (1.d0,0.d0), foo1R, NR,       &
                      V, n, (0.d0,0.d0), foo2R, NR)

!         ('Gr_Mn = - foo2R * Gr_nn')
          allocate (aux4(NR,norbDyn(idx)))
          foo3 = aux3(dim-n+1:dim,idxF(idx):idxL(idx))
          call zgemm ('N', 'N', NR, norbDyn(idx), n, (-1.d0,0.d0),      &
                      foo2R, NR, foo3, n, (0.d0,0.d0), aux4, NR)
          write (Gr_Mn_disk%lun) aux4
          deallocate (aux4)

!         Free memory.
          deallocate (foo3)
          deallocate (aux3)

          ueph = ueph + 1

       endif
    enddo

!   Assign auxiliary variables.
    dimbfr = dim ! last unit dimension
    utype = unit_type(I) ! current unit type
    dim = unitdimensions(utype) ! current type dimension
    dimaft = unitdimensions(ntypeunits+2) ! next unit dimension
    idx = ephIdx(utype) ! e-ph unit type

    if (idx /= 0) then

!      Allocate auxiliary matrix.
       allocate (aux3(dim,dim))

!      ('aux3 = E*S - H')
       aux3 = (Ei-unitshift(utype))*Sunits(utype)%S                     &
              - Hunits(utype)%H(:,:,ispin)

!      ('aux2 = GL_mm*V')
       call greenload (GL_mm_disk%lun, dimbfr, dimbfr, dimbfr-n+1,      &
                       dimbfr, dimbfr-n+1, dimbfr, aux1, n, n)
       call zsymm ('L', 'L', n, n, (1.d0,0.d0), aux1, n,                &
                   V, n, (0.d0,0.d0), aux2, n)

!      ('aux1 = V^dagger*aux2')
       call zgemm ('C', 'N', n, n, n, (1.d0,0.d0), V, n,                &
                   aux2, n, (0.d0,0.d0), aux1, n)

!      ('aux3 = aux3 - aux1')
       aux3(1:n,1:n) = aux3(1:n,1:n) - aux1

!      ('aux2 = V*GR_pp')
       call greenloadR (GR_pp_disk%lun, dimaft, dimaft, 1, n,           &
                        1, n, aux1, n, n, pos_pp(ueph))
       call zsymm ('R', 'L', n, n, (1.d0,0.d0), aux1, n,                &
                   V, n, (0.d0,0.d0), aux2, n)

!      ('aux1 = aux2*V^dagger')
       call zgemm ('N', 'C', n, n, n, (1.d0,0.d0), aux2, n,             &
                   V, n, (0.d0,0.d0), aux1, n)

!      ('aux3 = aux3 - aux1')
       aux3(dim-n+1:dim,dim-n+1:dim) = aux3(dim-n+1:dim,dim-n+1:dim)    &
                                       - aux1

!      ('aux3 = aux3^-1')
       allocate (ipiv(dim))
       call CHECKzsytrf (dim, 'L', aux3, ipiv)
       call CHECKzsytri (dim, 'L', aux3, ipiv)
       deallocate (ipiv)

!      ('Gr_nn = aux3')
       allocate (aux4(norbDyn(idx),norbDyn(idx)))
       aux4 = aux3(idxF(idx):idxL(idx),idxF(idx):idxL(idx))
       write (Gr_nn_disk%lun) aux4
       deallocate (aux4)

!      Allocate auxiliary matrix.
       allocate (foo3(n,norbDyn(idx)))

!      ('foo2L = GL_1m*V')
       call greenload (GL_1m_disk%lun, NL, dimbfr, 1, NL,               &
                       dimbfr-n+1, dimbfr, foo1L, NL, n)
       call zgemm ('N', 'N', NL, n, n, (1.d0,0.d0), foo1L, NL,          &
                   V, n, (0.d0,0.d0), foo2L, NL)

!      ('Gr_1n = - foo2L * Gr_nn')
       allocate (aux4(NL,norbDyn(idx)))
       foo3 = aux3(1:n,idxF(idx):idxL(idx))
       call zgemm ('N', 'N', NL, norbDyn(idx), n, (-1.d0,0.d0),         &
                   foo2L, NL, foo3, n, (0.d0,0.d0), aux4, NL)
       write (Gr_1n_disk%lun) aux4
       deallocate (aux4)

!      ('foo2R = GR_Mp*V^dagger')
       call greenloadR (GR_Mp_disk%lun, NR, dimaft, 1, NR,              &
                        1, n, foo1R, NR, n, pos_Mp(ueph))
       call zgemm ('N', 'C', NR, n, n, (1.d0,0.d0), foo1R, NR,          &
                   V, n, (0.d0,0.d0), foo2R, NR)

!      ('Gr_Mn = - foo2R * Gr_nn')
       allocate (aux4(NR,norbDyn(idx)))
       foo3 = aux3(dim-n+1:dim,idxF(idx):idxL(idx))
       call zgemm ('N', 'N', NR, norbDyn(idx), n, (-1.d0,0.d0),         &
                   foo2R, NR, foo3, n, (0.d0,0.d0), aux4, NR)
       write (Gr_Mn_disk%lun) aux4
       deallocate (aux4)

!      Free memory.
       deallocate (foo3)
       deallocate (aux3)

       ueph = ueph + 1

    endif

!   Assign auxiliary variables.
    dimbfr = dim ! last unit dimension
    utype = ntypeunits + 2 ! current unit type (last unit)
    dim = unitdimensions(utype) ! current type dimension
    idx = ephIdx(utype) ! e-ph unit type

    if (idx /= 0) then

!      Allocate auxiliary matrix.
       allocate (aux3(dim,dim))

!      ('aux3 = E*S - H')
       aux3 = (Ei-unitshift(utype))*Sunits(utype)%S                     &
              - Hunits(utype)%H(:,:,ispin)
       aux3(dim-NR+1:dim,dim-NR+1:dim) =                                &
                                        aux3(dim-NR+1:dim,dim-NR+1:dim) &
                                        - Sigma_R

!      ('aux2 = GL_mm*V')
       call greenload (GL_mm_disk%lun, dimbfr, dimbfr, dimbfr-n+1,      &
                       dimbfr, dimbfr-n+1, dimbfr, aux1, n, n)
       call zsymm ('L', 'L', n, n, (1.d0,0.d0), aux1, n,                &
                   V, n, (0.d0,0.d0), aux2, n)

!      ('aux1 = V^dagger*aux2')
       call zgemm ('C', 'N', n, n, n, (1.d0,0.d0), V, n,                &
                   aux2, n, (0.d0,0.d0), aux1, n)

!      ('aux3 = aux3 - aux1')
       aux3(1:n,1:n) = aux3(1:n,1:n) - aux1

!      ('aux3 = aux3^-1')
       allocate (ipiv(dim))
       call CHECKzsytrf (dim, 'L', aux3, ipiv)
       call CHECKzsytri (dim, 'L', aux3, ipiv)
       deallocate (ipiv)

!      ('Gr_nn = aux3')
       allocate (aux4(norbDyn(idx),norbDyn(idx)))
       aux4 = aux3(idxF(idx):idxL(idx),idxF(idx):idxL(idx))
       write (Gr_nn_disk%lun) aux4
       deallocate (aux4)

!      ('Gr_Mn = aux3')
       allocate (aux4(NR,norbDyn(idx)))
       aux4 = aux3(dim-NR+1:dim,idxF(idx):idxL(idx))
       write (Gr_Mn_disk%lun) aux4
       deallocate (aux4)

!      Allocate auxiliary matrix.
       allocate (foo3(n,dim))
       deallocate (aux1)
       allocate (aux1(NL,dim))
 
!      ('foo2L = GL_1m*V')
       call greenload (GL_1m_disk%lun, NL, dimbfr, 1, NL,               &
                       dimbfr-n+1, dimbfr, foo1L, NL, n)
       call zgemm ('N', 'N', NL, n, n, (1.d0,0.d0), foo1L, NL,          &
                   V, n, (0.d0,0.d0), foo2L, NL)

!      ('Gr_1n = - foo2L * Gr_nn')
       foo3 = aux3(1:n,1:dim)
       call zgemm ('N', 'N', NL, dim, n, (-1.d0,0.d0), foo2L, NL,       &
                   foo3, n, (0.d0,0.d0), aux1, NL)
       allocate (aux4(NL,norbDyn(idx)))
       aux4 = aux1(1:NL,idxF(idx):idxL(idx))
       write (Gr_1n_disk%lun) aux4
       deallocate (aux4)
       Gr_1M = aux1(1:NL,dim-NR+1:dim)

!      Free memory.
       deallocate (foo3)
       deallocate (aux3)

    else

!      Allocate auxiliary matrix.
       allocate (aux3(dim,dim))

!      ('aux3 = E*S - H')
       aux3 = (Ei-unitshift(utype))*Sunits(utype)%S                     &
              - Hunits(utype)%H(:,:,ispin)
       aux3(dim-NR+1:dim,dim-NR+1:dim) =                                &
                                        aux3(dim-NR+1:dim,dim-NR+1:dim) &
                                        - Sigma_R

!      ('aux2 = GL_nn*V')
       aux1 = GL_nn(dimbfr-n+1:dimbfr,dimbfr-n+1:dimbfr)
       call zsymm ('L', 'L', n, n, (1.d0,0.d0), aux1, n,                &
                   V, n, (0.d0,0.d0), aux2, n)

!      ('aux1 = V^dagger*aux2')
       call zgemm ('C', 'N', n, n, n, (1.d0,0.d0), V, n,                &
                   aux2, n, (0.d0,0.d0), aux1, n)

!      ('aux3 = aux3 - aux1')
       aux3(1:n,1:n) = aux3(1:n,1:n) - aux1

!      ('aux3 = aux3^-1')
       allocate (ipiv(dim))
       call CHECKzsytrf (dim, 'L', aux3, ipiv)
       call CHECKzsytri (dim, 'L', aux3, ipiv)
       deallocate (ipiv)

!      Allocate auxiliary matrix.
       allocate (foo3(n,NR))

!      ('foo2L = GL_1m*V')
       foo1L = GL_1N(1:NL,dimbfr-n+1:dimbfr)
       call zgemm ('N', 'N', NL, n, n, (1.d0,0.d0), foo1L, NL,          &
                   V, n, (0.d0,0.d0), foo2L, NL)

!      ('Gr_1n = - foo2L * G_nn')
       foo3 = aux3(1:n,dim-NR+1:dim)
       call zgemm ('N', 'N', NL, NR, n, (-1.d0,0.d0), foo2L, NL,        &
                   foo3, NR, (0.d0,0.d0), Gr_1M, NL)

!      Free memory.
       deallocate (foo3)
       deallocate (aux3)

    endif

!   Close Green's functions files.
    call IOcloseStream (GL_mm_disk%fname, GL_mm_disk%lun)
    call IOcloseStream (GL_1m_disk%fname, GL_1m_disk%lun)
    call IOcloseStream (GR_pp_disk%fname, GR_pp_disk%lun)
    call IOcloseStream (GR_Mp_disk%fname, GR_Mp_disk%lun)
    call IOcloseStream (Gr_nn_disk%fname, Gr_nn_disk%lun)
    call IOcloseStream (Gr_1n_disk%fname, Gr_1n_disk%lun)
    call IOcloseStream (Gr_Mn_disk%fname, Gr_Mn_disk%lun)

!   Free memory.
    deallocate (V)
    deallocate (aux1)
    deallocate (aux2)
    deallocate (foo1L)
    deallocate (foo2L)
    deallocate (foo1R)
    deallocate (foo2R)

    if (IOnode) write(6,'(a)') " ok!"


  end subroutine GFfullDisk


!  *******************************************************************  !
!                             greenFilesSet                             !
!  *******************************************************************  !
!  Description: assign the Green's function data structure (logical     !
!  unit number and system name file.                                    !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer Node                : Actual node (MPI_Comm_rank)            !
!  character(60) directory     : Working directory                      !
!  ****************************** INPUT ******************************  !
!  integer lun                 : Logical unit number                    !
!  character(*) name           : File name                              !
!  ***************************** OUTPUT ******************************  !
!  TYPE(greenfile) GF          : Green's function data structure        !
!  *******************************************************************  !
  subroutine greenFilesSet (lun, name, GF)

!
!   Modules
!
    use parallel,        only: Node
    use idsrdr_options,  only: directory
    use idsrdr_io,       only: IOopenStreamnew, IOcloseStream
    use idsrdr_string,   only: STRconcat, STRpaste

!   Input variables.
    integer, intent(in) :: lun
    character(*), intent(in) :: name
    TYPE(greenfile), intent(out) :: GF

!   Local variables.
    character(len=7) :: suffix

!   Set logical unit number.
    GF%lun = lun

!   Set file name.
    write (suffix,'(i3)') Node
    call STRconcat (suffix, '.GF', suffix)
    call STRpaste ('_', suffix, suffix)
    call STRpaste (name, suffix, GF%fname)
    call STRpaste (directory, GF%fname, GF%fname)

!   Create a new file.
    call IOopenStreamnew (GF%fname, GF%lun)
    call IOcloseStream (GF%fname, GF%lun)


  end subroutine greenFilesSet


!  *******************************************************************  !
!                               greenload                               !
!  *******************************************************************  !
!  Description: reads the required part of Green's function from file.  !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  ****************************** INPUT ******************************  !
!  integer GFlun               : Green's function logical unit number   !
!  integer rdim                : Green's function number of rows        !
!  integer cdim                : Green's function number of columns     !
!  integer rIni                : Required initial row index             !
!  integer rEnd                : Required final row index               !
!  integer cIni                : Required initial colunm index          !
!  integer cEnd                : Required final colunm index            !
!  integer rgf                 : Number of rows of output matrix        !
!  integer cgf                 : Number of columns of output matrix     !
!  ***************************** OUTPUT ******************************  !
!  complex*8 GFout(rgf,cgf)    : Required part of Green's function      !
!                                read from file                         !
!  *******************************************************************  !
  subroutine greenload (GFlun, rdim, cdim, rIni, rEnd, cIni, cEnd,      &
                        GFout, rgf, cgf)

!   Input variables.
    integer, intent(in) :: GFlun, rdim, cdim, rIni, rEnd, cIni, cEnd,   &
                           rgf, cgf
    complex(8), dimension (rgf,cgf), intent(out) :: GFout

!   Local variables.
    complex(8), allocatable, dimension (:,:) :: aux

!   Allocate auxiliary matrix.
    allocate (aux(rdim,cdim))

!   Read Green's function from file.
    read (GFlun) aux

!   Copy the required part.
    GFout = aux(rIni:rEnd,cIni:cEnd)

!   Free memory.
    deallocate (aux)


  end subroutine greenload


!  *******************************************************************  !
!                              greenloadR                               !
!  *******************************************************************  !
!  Description: reads the required part of Green's function from file   !
!  (for reading 'GR_pp' and 'GR_Mp' which are written backwards in the  !
!  file, so the 'position' parameter is needed).                        !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  ****************************** INPUT ******************************  !
!  integer GFlun               : Green's function logical unit number   !
!  integer rdim                : Green's function number of rows        !
!  integer cdim                : Green's function number of columns     !
!  integer rIni                : Required initial row index             !
!  integer rEnd                : Required final row index               !
!  integer cIni                : Required initial colunm index          !
!  integer cEnd                : Required final colunm index            !
!  integer rgf                 : Number of rows of output matrix        !
!  integer cgf                 : Number of columns of output matrix     !
!  integer position            : Position in the file                   !
!  ***************************** OUTPUT ******************************  !
!  complex*8 GFout(rgf,cgf)    : Required part of Green's function      !
!                                read from file                         !
!  *******************************************************************  !
  subroutine greenloadR (GFlun, rdim, cdim, rIni, rEnd, cIni, cEnd,     &
                         GFout, rgf, cgf, position)

!   Input variables.
    integer, intent(in) :: GFlun, rdim, cdim, rIni, rEnd, cIni, cEnd,   &
                           rgf, cgf, position
    complex(8), dimension (rgf,cgf), intent(out) :: GFout

!   Local variables.
    complex(8), allocatable, dimension (:,:) :: aux

!   Allocate auxiliary matrix.
    allocate (aux(rdim,cdim))

!   Read Green's function from file.
    read (GFlun, POS=position) aux

!   Copy the required part.
    GFout = aux(rIni:rEnd,cIni:cEnd)

!   Free memory.
    deallocate (aux)


  end subroutine greenloadR


!  *******************************************************************  !
!                                GFtest                                 !
!  *******************************************************************  !
!  Description: build the full hamiltonian of the scattering region     !
!  and invert it to obtain the entire Green's function in order to      !
!  compare with the Green's function obtained with the recursive        !
!  method.                                                              !
!                                                                       !
!  ATENTION: the writting part only works in serial mode!               !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode                       : True if it is the I/O node    !
!  integer ntypeunits                   : Number of unit types          !
!  integer nunits                       : Total number of units         !
!  integer unit_type(nunits+2)          : Units types                   !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  real*8 unitshift(ntypeunits+2)       : Units shift                   !
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
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Sigma_L(NL,NL)           : Left-lead self-energy          !
!  complex(8) Sigma_R(NR,NR)           : Right-lead self-energy         !
!  integer ephIdx(ntypeunits+2)        : Unit index (those with e-ph)   !
!  integer idxF(neph)                  : First dynamic atom orbital     !
!  integer idxL(neph)                  : Last dynamic atom orbital      !
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 Ei                            : Energy grid point             !
!  *******************************************************************  !
  subroutine GFtest (Ei, ispin)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_options,  only: ntypeunits, nunits
    use idsrdr_units,    only: unit_type, unitdimensions, unitshift,    &
                               S1unit, H1unit, Sunits, Hunits
    use idsrdr_leads,    only: NL, NR, Sigma_L, Sigma_R
    use idsrdr_ephcoupl, only: ephIdx, idxF, idxL
    use idsrdr_check,    only: CHECKzgetrf, CHECKzgetri

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    integer :: i, j, k, utype, dim, dimTot, dimCpl, idxAnt, idx, ueph
    integer, allocatable, dimension (:) :: ipiv
    real(8) :: dosTot
    real(8), parameter :: pi = 3.1415926535897932384626433832795028841D0
    real(8), allocatable, dimension (:,:) :: Stot
    complex(8), allocatable, dimension (:,:) :: Htot, Gtot

    if (IOnode) write (6,'(a)', advance='no')                           &
            '      computing TOTAL Greens function... '

!   Get the total dimension.
    dimTot = unitdimensions(ntypeunits+1) + unitdimensions(ntypeunits+2)
    do I = 2,nunits+1
       utype = unit_type(I) ! current unit type
       dimTot = dimTot + unitdimensions(utype)
    enddo

!   Allocate and initialize matrices.
    allocate (Stot(dimTot,dimTot))
    allocate (Htot(dimTot,dimTot))
    allocate (Gtot(dimTot,dimTot))
    allocate (ipiv(dimTot))
    Stot = 0.d0
    Htot = 0.d0

!   Build total hamiltonian and overlap matrices.
    dimCpl = unitdimensions(ntypeunits) ! coupling unit dimensions
    utype = ntypeunits+1 ! current unit type (first unit)
    dim = unitdimensions(utype) ! current unit dimensions
    idx = ephIdx(utype) ! e-ph unit type
    ueph = 1 ! eph units indexing

    Stot(1:dim,1:dim) = (Ei-unitshift(utype))*Sunits(utype)%S(:,:)
    Htot(1:dim,1:dim) = Hunits(utype)%H(:,:,ispin)
    Stot(dim-dimCpl+1:dim,dim+1:dim+dimCpl) =                           &
         (Ei-unitshift(ntypeunits))*S1unit
    Htot(dim-dimCpl+1:dim,dim+1:dim+dimCpl) = H1unit(:,:,ispin)
    Stot(dim+1:dim+dimCpl,dim-dimCpl+1:dim) =                           &
         (Ei-unitshift(ntypeunits))*TRANSPOSE(S1unit)
    Htot(dim+1:dim+dimCpl,dim-dimCpl+1:dim) =                           &
         TRANSPOSE(H1unit(:,:,ispin))

    idxAnt = dim + 1
    do k = 2,nunits+1
       utype = unit_type(k) ! current unit type
       dim = unitdimensions(utype) ! current unit dimensions

       Stot(idxAnt:idxAnt+dim-1,idxAnt:idxAnt+dim-1) =                  &
            (Ei-unitshift(utype))*Sunits(utype)%S(:,:)
       Htot(idxAnt:idxAnt+dim-1,idxAnt:idxAnt+dim-1) =                  &
            Hunits(utype)%H(:,:,ispin)
       Stot(idxAnt+dim-dimCpl:idxAnt+dim-1,                             &
            idxAnt+dim:idxAnt+dim+dimCpl-1) =                           &
            (Ei-unitshift(ntypeunits))*S1unit
       Htot(idxAnt+dim-dimCpl:idxAnt+dim-1,                             &
            idxAnt+dim:idxAnt+dim+dimCpl-1) = H1unit(:,:,ispin)
       Stot(idxAnt+dim:idxAnt+dim+dimCpl-1,                             &
            idxAnt+dim-dimCpl:idxAnt+dim-1) =                           &
            (Ei-unitshift(ntypeunits))*TRANSPOSE(S1unit)
       Htot(idxAnt+dim:idxAnt+dim+dimCpl-1,                             &
            idxAnt+dim-dimCpl:idxAnt+dim-1) =                           &
            TRANSPOSE(H1unit(:,:,ispin))

       idxAnt = idxAnt + dim

    enddo

    utype = ntypeunits + 2 ! current unit type (last unit)
    dim = unitdimensions(utype) ! current unit dimensions

    Stot(idxAnt:idxAnt+dim-1,idxAnt:idxAnt+dim-1) =                     &
         (Ei-unitshift(utype))*Sunits(utype)%S(:,:)
    Htot(idxAnt:idxAnt+dim-1,idxAnt:idxAnt+dim-1) =                     &
         Hunits(utype)%H(:,:,ispin)

!   Green's function.
    Gtot = Stot - Htot
    Gtot(1:NL,1:NL) = Gtot(1:NL,1:NL) - Sigma_L
    Gtot(dimTot-NR+1:dimTot,dimTot-NR+1:dimTot) =                       &
         Gtot(dimTot-NR+1:dimTot,dimTot-NR+1:dimTot) - Sigma_R
    call CHECKzgetrf (dimTot, Gtot, ipiv)
    call CHECKzgetri (dimTot, Gtot, ipiv)

!   Free memory.
    deallocate (Htot)
    deallocate (Stot)
    deallocate (ipiv)

    if (IOnode) write(6,'(a)') " ok!"

!   Compute the total density of states.
    dosTot = 0.d0
    do i = 1,dimTot
       dosTot = dosTot - DIMAG(Gtot(i,i))
    enddo
!!$    dosTot = dosTot / pi
    write (4102,'(e17.8e3,e17.8e3)') Ei, dosTot

!   Write everything...
    do i = 1,dimTot
       do j = 1,dimTot
          write (102,'(i5,e17.8e3,e17.8e3)') j, DREAL(Gtot(i,j)),       &
               DIMAG(Gtot(i,j)) 
          write (202,'(i5,e17.8e3)') j, CDABS(Gtot(i,j))
       enddo
    enddo

!   First unit.
    if (idx /= 0) then

       do j = idxF(idx),idxL(idx)
          do i = idxF(idx),idxL(idx)
             write (3102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Gtot(i,j)),                                     &
                  DREAL(Gr_nn(ueph)%G(i-idxF(idx)+1,j-idxF(idx)+1)),    &
                  DIMAG(Gtot(i,j)),                                     &
                  DIMAG(Gr_nn(ueph)%G(i-idxF(idx)+1,j-idxF(idx)+1))
             write (2102,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Gtot(i,j)) -                                    &
                  DREAL(Gr_nn(ueph)%G(i-idxF(idx)+1,j-idxF(idx)+1)),    &
                  DIMAG(Gtot(i,j)) -                                    &
                  DIMAG(Gr_nn(ueph)%G(i-idxF(idx)+1,j-idxF(idx)+1))
          enddo
       enddo
       do j = idxF(idx),idxL(idx)
          do i = 1,NL
             write (6102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Gtot(i,j)),                                     &
                  DREAL(Gr_1n(ueph)%G(i,j-idxF(idx)+1)),                &
                  DIMAG(Gtot(i,j)),                                     &
                  DIMAG(Gr_1n(ueph)%G(i,j-idxF(idx)+1))
             write (5102,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Gtot(i,j)) -                                    &
                  DREAL(Gr_1n(ueph)%G(i,j-idxF(idx)+1)),                &
                  DIMAG(Gtot(i,j)) -                                    &
                  DIMAG(Gr_1n(ueph)%G(i,j-idxF(idx)+1))
          enddo
       enddo
       do j = idxF(idx),idxL(idx)
          do i = 1,NR
             write (8102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Gtot(dimTot-NR+i,j)),                           &
                  DREAL(Gr_Mn(ueph)%G(i,j-idxF(idx)+1)),                &
                  DIMAG(Gtot(dimTot-NR+i,j)),                           &
                  DIMAG(Gr_Mn(ueph)%G(i,j-idxF(idx)+1))
             write (7102,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Gtot(dimTot-NR+i,j)) -                          &
                  DREAL(Gr_Mn(ueph)%G(i,j-idxF(idx)+1)),                &
                  DIMAG(Gtot(dimTot-NR+i,j)) -                          &
                  DIMAG(Gr_Mn(ueph)%G(i,j-idxF(idx)+1))
          enddo
       enddo

       ueph = ueph + 1

    endif

!   Units in the middle.
    idxAnt = unitdimensions(ntypeunits+1) ! last unit dimensions
    do k = 2,nunits+1

       utype = unit_type(k) ! current unit type
       dim = unitdimensions(utype) ! current unit dimensions
       idx = ephIdx(utype) ! e-ph unit type

       if (idx /= 0) then

          do j = idxF(idx),idxL(idx)
             do i = idxF(idx),idxL(idx)
                write (3102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')        &
                     DREAL(Gtot(idxAnt+i,idxAnt+j)),                    &
                     DREAL(Gr_nn(ueph)%G(i-idxF(idx)+1,j-idxF(idx)+1)), &
                     DIMAG(Gtot(idxAnt+i,idxAnt+j)),                    &
                     DIMAG(Gr_nn(ueph)%G(i-idxF(idx)+1,j-idxF(idx)+1))
                write (2102,'(e17.8e3,e17.8e3)')                        &
                     DREAL(Gtot(idxAnt+i,idxAnt+j)) -                   &
                     DREAL(Gr_nn(ueph)%G(i-idxF(idx)+1,j-idxF(idx)+1)), &
                     DIMAG(Gtot(idxAnt+i,idxAnt+j)) -                   &
                     DIMAG(Gr_nn(ueph)%G(i-idxF(idx)+1,j-idxF(idx)+1))
             enddo
          enddo
          do j = idxF(idx),idxL(idx)
             do i = 1,NL
                write (6102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')        &
                     DREAL(Gtot(i,idxAnt+j)),                           &
                     DREAL(Gr_1n(ueph)%G(i,j-idxF(idx)+1)),             &
                     DIMAG(Gtot(i,idxAnt+j)),                           &
                     DIMAG(Gr_1n(ueph)%G(i,j-idxF(idx)+1))
                write (5102,'(e17.8e3,e17.8e3)')                        &
                     DREAL(Gtot(i,idxAnt+j)) -                          &
                     DREAL(Gr_1n(ueph)%G(i,j-idxF(idx)+1)),             &
                     DIMAG(Gtot(i,idxAnt+j)) -                          &
                     DIMAG(Gr_1n(ueph)%G(i,j-idxF(idx)+1))
             enddo
          enddo
          do j = idxF(idx),idxL(idx)
             do i = 1,NR
                write (8102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')        &
                     DREAL(Gtot(dimTot-NR+i,idxAnt+j)),                 &
                     DREAL(Gr_Mn(ueph)%G(i,j-idxF(idx)+1)),             &
                     DIMAG(Gtot(dimTot-NR+i,idxAnt+j)),                 &
                     DIMAG(Gr_Mn(ueph)%G(i,j-idxF(idx)+1))
                write (7102,'(e17.8e3,e17.8e3)')                        &
                     DREAL(Gtot(dimTot-NR+i,idxAnt+j)) -                &
                     DREAL(Gr_Mn(ueph)%G(i,j-idxF(idx)+1)),             &
                     DIMAG(Gtot(dimTot-NR+i,idxAnt+j)) -                &
                     DIMAG(Gr_Mn(ueph)%G(i,j-idxF(idx)+1))
             enddo
          enddo

          ueph = ueph + 1

       endif

       idxAnt = idxAnt + dim

    enddo

!   Last unit.
    idx = ephIdx(ntypeunits+2) ! e-ph unit type
    if (idx /= 0) then

       do j = idxF(idx),idxL(idx)
          do i = idxF(idx),idxL(idx)
             write (3102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Gtot(idxAnt+i,idxAnt+j)),                       &
                  DREAL(Gr_nn(ueph)%G(i-idxF(idx)+1,j-idxF(idx)+1)),    &
                  DIMAG(Gtot(idxAnt+i,idxAnt+j)),                       &
                  DIMAG(Gr_nn(ueph)%G(i-idxF(idx)+1,j-idxF(idx)+1))
             write (2102,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Gtot(idxAnt+i,idxAnt+j)) -                      &
                  DREAL(Gr_nn(ueph)%G(i-idxF(idx)+1,j-idxF(idx)+1)),    &
                  DIMAG(Gtot(idxAnt+i,idxAnt+j)) -                      &
                  DIMAG(Gr_nn(ueph)%G(i-idxF(idx)+1,j-idxF(idx)+1))
          enddo
       enddo
       do j = idxF(idx),idxL(idx)
          do i = 1,NL
             write (6102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Gtot(i,idxAnt+j)),                              &
                  DREAL(Gr_1n(ueph)%G(i,j-idxF(idx)+1)),                &
                  DIMAG(Gtot(i,idxAnt+j)),                              &
                  DIMAG(Gr_1n(ueph)%G(i,j-idxF(idx)+1))
             write (5102,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Gtot(i,idxAnt+j)) -                             &
                  DREAL(Gr_1n(ueph)%G(i,j-idxF(idx)+1)),                &
                  DIMAG(Gtot(i,idxAnt+j)) -                             &
                  DIMAG(Gr_1n(ueph)%G(i,j-idxF(idx)+1))
          enddo
       enddo
       do j = idxF(idx),idxL(idx)
          do i = 1,NR
             write (8102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Gtot(dimTot-NR+i,idxAnt+j)),                    &
                  DREAL(Gr_Mn(ueph)%G(i,j-idxF(idx)+1)),                &
                  DIMAG(Gtot(dimTot-NR+i,idxAnt+j)),                    &
                  DIMAG(Gr_Mn(ueph)%G(i,j-idxF(idx)+1))
             write (7102,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Gtot(dimTot-NR+i,idxAnt+j)) -                   &
                  DREAL(Gr_Mn(ueph)%G(i,j-idxF(idx)+1)),                &
                  DIMAG(Gtot(dimTot-NR+i,idxAnt+j)) -                   &
                  DIMAG(Gr_Mn(ueph)%G(i,j-idxF(idx)+1))
          enddo
       enddo

    endif

    do j = 1,NR
       do i = 1,NL
          write (310,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')               &
               DREAL(Gtot(i,dimTot-NR+j)), DREAL(Gr_1M(i,j)),           &
               DIMAG(Gtot(i,dimTot-NR+j)), DIMAG(Gr_1M(i,j))
          write (210,'(e17.8e3,e17.8e3)')                               &
               DREAL(Gtot(i,dimTot-NR+j)) - DREAL(Gr_1M(i,j)),          &
               DIMAG(Gtot(i,dimTot-NR+j)) - DIMAG(Gr_1M(i,j))
       enddo
    enddo

!   Free memory.
    deallocate (Gtot)


  end subroutine GFtest


!  *******************************************************************  !
!                              GFtestDisk                               !
!  *******************************************************************  !
!  Description: build the full hamiltonian of the scattering region     !
!  and invert it to obtain the entire Green's function in order to      !
!  compare with the Green's function obtained with the recursive        !
!  method. (Green's functions are read from disk).                      !
!                                                                       !
!  ATENTION: the writting part only works in serial mode!               !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode                       : True if it is the I/O node    !
!  integer ntypeunits                   : Number of unit types          !
!  integer nunits                       : Total number of units         !
!  integer unit_type(nunits+2)          : Units types                   !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  real*8 unitshift(ntypeunits+2)       : Units shift                   !
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
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Sigma_L(NL,NL)           : Left-lead self-energy          !
!  complex(8) Sigma_R(NR,NR)           : Right-lead self-energy         !
!  integer ephIdx(ntypeunits+2)        : Unit index (those with e-ph)   !
!  integer norbDyn(neph)               : Number of orbitals from        !
!                                        dynamic atoms                  !
!  integer idxF(neph)                  : First dynamic atom orbital     !
!  integer idxL(neph)                  : Last dynamic atom orbital      !
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 Ei                            : Energy grid point             !
!  *******************************************************************  !
  subroutine GFtestDisk (Ei, ispin)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_options,  only: ntypeunits, nunits
    use idsrdr_units,    only: unit_type, unitdimensions, unitshift,    &
                               S1unit, H1unit, Sunits, Hunits
    use idsrdr_leads,    only: NL, NR, Sigma_L, Sigma_R
    use idsrdr_ephcoupl, only: ephIdx, norbDyn, idxF, idxL
    use idsrdr_check,    only: CHECKzgetrf, CHECKzgetri
    use idsrdr_io,       only: IOopenStream, IOcloseStream

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    integer :: i, j, k, utype, dim, dimTot, dimCpl, idxAnt, idx, ueph
    integer, allocatable, dimension (:) :: ipiv
    real(8) :: dosTot
    real(8), parameter :: pi = 3.1415926535897932384626433832795028841D0
    real(8), allocatable, dimension (:,:) :: Stot
    complex(8), allocatable, dimension (:,:) :: Htot, Gtot, aux

    if (IOnode) write (6,'(a)', advance='no')                           &
            '      computing TOTAL Greens function... '

!   Get the total dimension.
    dimTot = unitdimensions(ntypeunits+1) + unitdimensions(ntypeunits+2)
    do I = 2,nunits+1
       utype = unit_type(I) ! current unit type
       dimTot = dimTot + unitdimensions(utype)
    enddo

!   Allocate and initialize matrices.
    allocate (Stot(dimTot,dimTot))
    allocate (Htot(dimTot,dimTot))
    allocate (Gtot(dimTot,dimTot))
    allocate (ipiv(dimTot))
    Stot = 0.d0
    Htot = 0.d0

!   Build total hamiltonian and overlap matrices.
    dimCpl = unitdimensions(ntypeunits) ! coupling unit dimensions
    utype = ntypeunits+1 ! current unit type (first unit)
    dim = unitdimensions(utype) ! current unit dimensions
    idx = ephIdx(utype) ! e-ph unit type
    ueph = 1 ! eph units indexing

    Stot(1:dim,1:dim) = (Ei-unitshift(utype))*Sunits(utype)%S(:,:)
    Htot(1:dim,1:dim) = Hunits(utype)%H(:,:,ispin)
    Stot(dim-dimCpl+1:dim,dim+1:dim+dimCpl) =                           &
         (Ei-unitshift(ntypeunits))*S1unit
    Htot(dim-dimCpl+1:dim,dim+1:dim+dimCpl) = H1unit(:,:,ispin)
    Stot(dim+1:dim+dimCpl,dim-dimCpl+1:dim) =                           &
         (Ei-unitshift(ntypeunits))*TRANSPOSE(S1unit)
    Htot(dim+1:dim+dimCpl,dim-dimCpl+1:dim) =                           &
         TRANSPOSE(H1unit(:,:,ispin))

    idxAnt = dim + 1
    do k = 2,nunits+1
       utype = unit_type(k) ! current unit type
       dim = unitdimensions(utype) ! current unit dimensions

       Stot(idxAnt:idxAnt+dim-1,idxAnt:idxAnt+dim-1) =                  &
            (Ei-unitshift(utype))*Sunits(utype)%S(:,:)
       Htot(idxAnt:idxAnt+dim-1,idxAnt:idxAnt+dim-1) =                  &
            Hunits(utype)%H(:,:,ispin)
       Stot(idxAnt+dim-dimCpl:idxAnt+dim-1,                             &
            idxAnt+dim:idxAnt+dim+dimCpl-1) =                           &
            (Ei-unitshift(ntypeunits))*S1unit
       Htot(idxAnt+dim-dimCpl:idxAnt+dim-1,                             &
            idxAnt+dim:idxAnt+dim+dimCpl-1) = H1unit(:,:,ispin)
       Stot(idxAnt+dim:idxAnt+dim+dimCpl-1,                             &
            idxAnt+dim-dimCpl:idxAnt+dim-1) =                           &
            (Ei-unitshift(ntypeunits))*TRANSPOSE(S1unit)
       Htot(idxAnt+dim:idxAnt+dim+dimCpl-1,                             &
            idxAnt+dim-dimCpl:idxAnt+dim-1) =                           &
            TRANSPOSE(H1unit(:,:,ispin))

       idxAnt = idxAnt + dim

    enddo

    utype = ntypeunits + 2 ! current unit type (last unit)
    dim = unitdimensions(utype) ! current unit dimensions

    Stot(idxAnt:idxAnt+dim-1,idxAnt:idxAnt+dim-1) =                     &
         (Ei-unitshift(utype))*Sunits(utype)%S(:,:)
    Htot(idxAnt:idxAnt+dim-1,idxAnt:idxAnt+dim-1) =                     &
         Hunits(utype)%H(:,:,ispin)

!   Green's function.
    Gtot = Stot - Htot
    Gtot(1:NL,1:NL) = Gtot(1:NL,1:NL) - Sigma_L
    Gtot(dimTot-NR+1:dimTot,dimTot-NR+1:dimTot) =                       &
         Gtot(dimTot-NR+1:dimTot,dimTot-NR+1:dimTot) - Sigma_R
    call CHECKzgetrf (dimTot, Gtot, ipiv)
    call CHECKzgetri (dimTot, Gtot, ipiv)

!   Free memory.
    deallocate (Htot)
    deallocate (Stot)
    deallocate (ipiv)

    if (IOnode) write(6,'(a)') " ok!"

!   Compute the total density of states.
    dosTot = 0.d0
    do i = 1,dimTot
       dosTot = dosTot - DIMAG(Gtot(i,i))
    enddo
!!$    dosTot = dosTot / pi
    write (4102,'(e17.8e3,e17.8e3)') Ei, dosTot

!   Write everything...
    do i = 1,dimTot
       do j = 1,dimTot
          write (102,'(i5,e17.8e3,e17.8e3)') j, DREAL(Gtot(i,j)),       &
               DIMAG(Gtot(i,j)) 
          write (202,'(i5,e17.8e3)') j, CDABS(Gtot(i,j))
       enddo
    enddo

!   Open Green's functions files.
    call IOopenStream (Gr_nn_disk%fname, Gr_nn_disk%lun)
    call IOopenStream (Gr_1n_disk%fname, Gr_1n_disk%lun)
    call IOopenStream (Gr_Mn_disk%fname, Gr_Mn_disk%lun)

!   First unit.
    if (idx /= 0) then

!      Auxiliary matrix to copy 'Gr_nn' from file.
       allocate (aux(norbDyn(idx),norbDyn(idx)))
       call greenload (Gr_nn_disk%lun, norbDyn(idx), norbDyn(idx),      &
                       1, norbDyn(idx), 1, norbDyn(idx), aux,           &
                       norbDyn(idx), norbDyn(idx))

       do j = idxF(idx),idxL(idx)
          do i = idxF(idx),idxL(idx)
             write (3102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Gtot(i,j)),                                     &
                  DREAL(aux(i-idxF(idx)+1,j-idxF(idx)+1)),              &
                  DIMAG(Gtot(i,j)),                                     &
                  DIMAG(aux(i-idxF(idx)+1,j-idxF(idx)+1))
             write (2102,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Gtot(i,j)) -                                    &
                  DREAL(aux(i-idxF(idx)+1,j-idxF(idx)+1)),              &
                  DIMAG(Gtot(i,j)) -                                    &
                  DIMAG(aux(i-idxF(idx)+1,j-idxF(idx)+1))
          enddo
       enddo

!      Free memory.
       deallocate (aux)

!      Auxiliary matrix to copy 'Gr_1n' from file.
       allocate (aux(NL,norbDyn(idx)))
       call greenload (Gr_1n_disk%lun, NL, norbDyn(idx), 1, NL,         &
                       1, norbDyn(idx), aux, NL, norbDyn(idx))

       do j = idxF(idx),idxL(idx)
          do i = 1,NL
             write (6102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Gtot(i,j)), DREAL(aux(i,j-idxF(idx)+1)),        &
                  DIMAG(Gtot(i,j)), DIMAG(aux(i,j-idxF(idx)+1))
             write (5102,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Gtot(i,j)) - DREAL(aux(i,j-idxF(idx)+1)),       &
                  DIMAG(Gtot(i,j)) - DIMAG(aux(i,j-idxF(idx)+1))
          enddo
       enddo

!      Free memory.
       deallocate (aux)

!      Auxiliary matrix to copy 'Gr_Mn' from file.
       allocate (aux(NR,norbDyn(idx)))
       call greenload (Gr_Mn_disk%lun, NR, norbDyn(idx), 1, NR,         &
                       1, norbDyn(idx), aux, NR, norbDyn(idx))

       do j = idxF(idx),idxL(idx)
          do i = 1,NR
             write (8102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Gtot(dimTot-NR+i,j)),                           &
                  DREAL(aux(i,j-idxF(idx)+1)),                          &
                  DIMAG(Gtot(dimTot-NR+i,j)),                           &
                  DIMAG(aux(i,j-idxF(idx)+1))
             write (7102,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Gtot(dimTot-NR+i,j)) -                          &
                  DREAL(aux(i,j-idxF(idx)+1)),                          &
                  DIMAG(Gtot(dimTot-NR+i,j)) -                          &
                  DIMAG(aux(i,j-idxF(idx)+1))
          enddo
       enddo

!      Free memory.
       deallocate (aux)

       ueph = ueph + 1

    endif

!   Units in the middle.
    idxAnt = unitdimensions(ntypeunits+1) ! last unit dimensions
    do k = 2,nunits+1

       utype = unit_type(k) ! current unit type
       dim = unitdimensions(utype) ! current unit dimensions
       idx = ephIdx(utype) ! e-ph unit type

       if (idx /= 0) then

!         Auxiliary matrix to copy 'Gr_nn' from file.
          allocate (aux(norbDyn(idx),norbDyn(idx)))
          call greenload (Gr_nn_disk%lun, norbDyn(idx), norbDyn(idx),   &
                          1, norbDyn(idx), 1, norbDyn(idx), aux,        &
                          norbDyn(idx), norbDyn(idx))

          do j = idxF(idx),idxL(idx)
             do i = idxF(idx),idxL(idx)
                write (3102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')        &
                     DREAL(Gtot(idxAnt+i,idxAnt+j)),                    &
                     DREAL(aux(i-idxF(idx)+1,j-idxF(idx)+1)),           &
                     DIMAG(Gtot(idxAnt+i,idxAnt+j)),                    &
                     DIMAG(aux(i-idxF(idx)+1,j-idxF(idx)+1))
                write (2102,'(e17.8e3,e17.8e3)')                        &
                     DREAL(Gtot(idxAnt+i,idxAnt+j)) -                   &
                     DREAL(aux(i-idxF(idx)+1,j-idxF(idx)+1)),           &
                     DIMAG(Gtot(idxAnt+i,idxAnt+j)) -                   &
                     DIMAG(aux(i-idxF(idx)+1,j-idxF(idx)+1))
             enddo
          enddo

!         Free memory.
          deallocate (aux)

!         Auxiliary matrix to copy 'Gr_1n' from file.
          allocate (aux(NL,norbDyn(idx)))
          call greenload (Gr_1n_disk%lun, NL, norbDyn(idx), 1, NL,      &
                          1, norbDyn(idx), aux, NL, norbDyn(idx))

          do j = idxF(idx),idxL(idx)
             do i = 1,NL
                write (6102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')        &
                     DREAL(Gtot(i,idxAnt+j)),                           &
                     DREAL(aux(i,j-idxF(idx)+1)),                       &
                     DIMAG(Gtot(i,idxAnt+j)),                           &
                     DIMAG(aux(i,j-idxF(idx)+1))
                write (5102,'(e17.8e3,e17.8e3)')                        &
                     DREAL(Gtot(i,idxAnt+j)) -                          &
                     DREAL(aux(i,j-idxF(idx)+1)),                       &
                     DIMAG(Gtot(i,idxAnt+j)) -                          &
                     DIMAG(aux(i,j-idxF(idx)+1))
             enddo
          enddo

!         Free memory.
          deallocate (aux)

!         Auxiliary matrix to copy 'Gr_Mn' from file.
          allocate (aux(NR,norbDyn(idx)))
          call greenload (Gr_Mn_disk%lun, NR, norbDyn(idx), 1, NR,      &
                          1, norbDyn(idx), aux, NR, norbDyn(idx))

          do j = idxF(idx),idxL(idx)
             do i = 1,NR
                write (8102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')        &
                     DREAL(Gtot(dimTot-NR+i,idxAnt+j)),                 &
                     DREAL(aux(i,j-idxF(idx)+1)),                       &
                     DIMAG(Gtot(dimTot-NR+i,idxAnt+j)),                 &
                     DIMAG(aux(i,j-idxF(idx)+1))
                write (7102,'(e17.8e3,e17.8e3)')                        &
                     DREAL(Gtot(dimTot-NR+i,idxAnt+j)) -                &
                     DREAL(aux(i,j-idxF(idx)+1)),                       &
                     DIMAG(Gtot(dimTot-NR+i,idxAnt+j)) -                &
                     DIMAG(aux(i,j-idxF(idx)+1))
             enddo
          enddo

!         Free memory.
          deallocate (aux)

          ueph = ueph + 1

       endif

       idxAnt = idxAnt + dim

    enddo

!   Last unit.
    idx = ephIdx(ntypeunits+2) ! e-ph unit type
    if (idx /= 0) then

!      Auxiliary matrix to copy 'Gr_nn' from file.
       allocate (aux(norbDyn(idx),norbDyn(idx)))
       call greenload (Gr_nn_disk%lun, norbDyn(idx), norbDyn(idx),      &
                       1, norbDyn(idx), 1, norbDyn(idx), aux,           &
                       norbDyn(idx), norbDyn(idx))

       do j = idxF(idx),idxL(idx)
          do i = idxF(idx),idxL(idx)
             write (3102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Gtot(idxAnt+i,idxAnt+j)),                       &
                  DREAL(aux(i-idxF(idx)+1,j-idxF(idx)+1)),    &
                  DIMAG(Gtot(idxAnt+i,idxAnt+j)),                       &
                  DIMAG(aux(i-idxF(idx)+1,j-idxF(idx)+1))
             write (2102,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Gtot(idxAnt+i,idxAnt+j)) -                      &
                  DREAL(aux(i-idxF(idx)+1,j-idxF(idx)+1)),    &
                  DIMAG(Gtot(idxAnt+i,idxAnt+j)) -                      &
                  DIMAG(aux(i-idxF(idx)+1,j-idxF(idx)+1))
          enddo
       enddo

!      Free memory.
       deallocate (aux)

!      Auxiliary matrix to copy 'Gr_1n' from file.
       allocate (aux(NL,norbDyn(idx)))
       call greenload (Gr_1n_disk%lun, NL, norbDyn(idx), 1, NL,         &
                       1, norbDyn(idx), aux, NL, norbDyn(idx))

       do j = idxF(idx),idxL(idx)
          do i = 1,NL
             write (6102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Gtot(i,idxAnt+j)), DREAL(aux(i,j-idxF(idx)+1)), &
                  DIMAG(Gtot(i,idxAnt+j)), DIMAG(aux(i,j-idxF(idx)+1))
             write (5102,'(e17.8e3,e17.8e3)')                           &
                 DREAL(Gtot(i,idxAnt+j)) - DREAL(aux(i,j-idxF(idx)+1)), &
                 DIMAG(Gtot(i,idxAnt+j)) - DIMAG(aux(i,j-idxF(idx)+1))
          enddo
       enddo

!      Free memory.
       deallocate (aux)

!      Auxiliary matrix to copy 'Gr_Mn' from file.
       allocate (aux(NR,norbDyn(idx)))
       call greenload (Gr_Mn_disk%lun, NR, norbDyn(idx), 1, NR,         &
                       1, norbDyn(idx), aux, NR, norbDyn(idx))

       do j = idxF(idx),idxL(idx)
          do i = 1,NR
             write (8102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Gtot(dimTot-NR+i,idxAnt+j)),                    &
                  DREAL(aux(i,j-idxF(idx)+1)),                          &
                  DIMAG(Gtot(dimTot-NR+i,idxAnt+j)),                    &
                  DIMAG(aux(i,j-idxF(idx)+1))
             write (7102,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Gtot(dimTot-NR+i,idxAnt+j)) -                   &
                  DREAL(aux(i,j-idxF(idx)+1)),                          &
                  DIMAG(Gtot(dimTot-NR+i,idxAnt+j)) -                   &
                  DIMAG(aux(i,j-idxF(idx)+1))
          enddo
       enddo

!      Free memory.
       deallocate (aux)

    endif

    do j = 1,NR
       do i = 1,NL
          write (310,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')               &
               DREAL(Gtot(i,dimTot-NR+j)), DREAL(Gr_1M(i,j)),           &
               DIMAG(Gtot(i,dimTot-NR+j)), DIMAG(Gr_1M(i,j))
          write (210,'(e17.8e3,e17.8e3)')                               &
               DREAL(Gtot(i,dimTot-NR+j)) - DREAL(Gr_1M(i,j)),          &
               DIMAG(Gtot(i,dimTot-NR+j)) - DIMAG(Gr_1M(i,j))
       enddo
    enddo

!   Close Green's functions files.
    call IOcloseStream (Gr_nn_disk%fname, Gr_nn_disk%lun)
    call IOcloseStream (Gr_1n_disk%fname, Gr_1n_disk%lun)
    call IOcloseStream (Gr_Mn_disk%fname, Gr_Mn_disk%lun)

!   Free memory.
    deallocate (Gtot)


  end subroutine GFtestDisk


!  *******************************************************************  !
!                               freegreen                               !
!  *******************************************************************  !
!  Description: free allocated vectors.                                 !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer nunitseph                   : Number of units with eph       !
!  *******************************************************************  !
  subroutine freegreen

!
!   Modules
!
    use idsrdr_units,    only: nunitseph
    use idsrdr_options,  only: writeondisk

!   Local variables.
    integer :: I

    IF (writeondisk) THEN

       deallocate (pos_pp, pos_Mp)

    ELSE ! wrote on memory...

!      First deallocates pointed matrices.
       do I = 1,nunitseph
          deallocate (GL_mm(I)%G)
          deallocate (GL_1m(I)%G)
          deallocate (GR_pp(I)%G)
          deallocate (GR_Mp(I)%G)
          deallocate (Gr_nn(I)%G)
          deallocate (Gr_1n(I)%G)
          deallocate (Gr_Mn(I)%G)
       enddo
       deallocate (GL_mm)
       deallocate (GL_1m)
       deallocate (GR_pp)
       deallocate (GR_Mp)
       deallocate (Gr_nn)
       deallocate (Gr_1n)
       deallocate (Gr_Mn)

    ENDIF ! IF (writeondisk)

    if (allocated(GL_nn)) then
       deallocate (GL_nn)
       deallocate (GL_1N)
    endif
    deallocate (Gr_1M)


  end subroutine freegreen


!  *******************************************************************  !


END MODULE idsrdr_green
