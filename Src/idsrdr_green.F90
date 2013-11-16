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

  implicit none
  
  PUBLIC  :: greeninit, greenfunctions, freegreen
  PRIVATE :: LRsweep , RLsweep, GFfull, GFtest

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

  complex(8), allocatable, dimension (:,:) :: GL_nn ! G^L_{n,n}
  complex(8), allocatable, dimension (:,:) :: GL_1N ! G^L_{1,M}
  complex(8), allocatable, dimension (:,:) :: Gr_1M ! G^r_{1,M}


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
!  integer unit_type(nunits+2)          : Units types                   !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  integer neph                        : Number of units with e-ph      !
!                                        interaction                    !
!  integer ephIdx(ntypeunits+2)        : Unit index (those with e-ph)   !
!  integer norbDyn(neph)               : Number of orbitals from        !
!                                        dynamic atoms                  !
!  *******************************************************************  !
  subroutine greeninit

!
!   Modules
!
    use idsrdr_options,  only: ntypeunits, nunits
    use idsrdr_units,    only: unit_type, unitdimensions
    use idsrdr_leads,    only: NL, NR
    use idsrdr_ephcoupl, only: neph, ephIdx, norbDyn

!   Local variables.
    integer :: I, J, dim, dimbfr, dimaft, ephType

!   Allocate Green's functions pointer array.
    allocate (GL_mm(neph))
    allocate (GL_1m(neph))
    allocate (GR_pp(neph))
    allocate (GR_Mp(neph))
    allocate (Gr_nn(neph))
    allocate (Gr_1n(neph))
    allocate (Gr_Mn(neph))

! DUVIDA: GL_1m tem que ter dimensão NL mesmo?
! E GR_Mp tem que ter dimensão NL mesmo?

!   Allocate Green's functions matrices.
    dim = unitdimensions(ntypeunits+1) ! current unit dimension
    dimbfr = NL ! last unit dimension
    dimaft = unitdimensions(unit_type(2)) ! next unit dimension
    ephType = ephIdx(ntypeunits+1) ! e-ph unit type
    J = 1
    if (ephType /= 0) then
       allocate (GL_mm(J)%G(dimbfr,dimbfr))
       allocate (GL_1m(J)%G(NL,dimbfr))
       allocate (GR_pp(J)%G(dimaft,dimaft))
       allocate (GR_Mp(J)%G(NR,dimaft))
       allocate (Gr_nn(J)%G(norbDyn(ephType),norbDyn(ephType)))
       allocate (Gr_1n(J)%G(NL,norbDyn(ephType)))
       allocate (Gr_Mn(J)%G(NR,norbDyn(ephType)))
       J = J + 1
    endif

    do I = 2,nunits
       dimbfr = dim ! last unit dimension
       dim = unitdimensions(unit_type(I)) ! current unit dimension
       dimaft = unitdimensions(unit_type(I+1)) ! next unit dimension
       ephType = ephIdx(unit_type(I)) ! e-ph unit type
       if (ephType /= 0) then
          allocate (GL_mm(J)%G(dimbfr,dimbfr))
          allocate (GL_1m(J)%G(NL,dimbfr))
          allocate (GR_pp(J)%G(dimaft,dimaft))
          allocate (GR_Mp(J)%G(NR,dimaft))
          allocate (Gr_nn(J)%G(norbDyn(ephType),norbDyn(ephType)))
          allocate (Gr_1n(J)%G(NL,norbDyn(ephType)))
          allocate (Gr_Mn(J)%G(NR,norbDyn(ephType)))
          J = J + 1
       endif
    enddo

    dimbfr = dim ! last unit dimension
    dim = unitdimensions(unit_type(I)) ! current unit dimension
    dimaft = unitdimensions(ntypeunits+2) ! next unit dimension
    ephType = ephIdx(unit_type(I)) ! e-ph unit type
    if (ephType /= 0) then
       allocate (GL_mm(J)%G(dimbfr,dimbfr))
       allocate (GL_1m(J)%G(NL,dimbfr))
       allocate (GR_pp(J)%G(dimaft,dimaft))
       allocate (GR_Mp(J)%G(NR,dimaft))
       allocate (Gr_nn(J)%G(norbDyn(ephType),norbDyn(ephType)))
       allocate (Gr_1n(J)%G(NL,norbDyn(ephType)))
       allocate (Gr_Mn(J)%G(NR,norbDyn(ephType)))
       J = J + 1
    endif

    dimbfr = dim ! last unit dimension
    dim = unitdimensions(ntypeunits+2) ! current unit dimension
    dimaft = NR ! next unit dimension
    ephType = ephIdx(ntypeunits+2) ! e-ph unit type
    if (ephType /= 0) then
       allocate (GL_mm(J)%G(dimbfr,dimbfr))
       allocate (GL_1m(J)%G(NL,dimbfr))
       allocate (GR_pp(J)%G(dimaft,dimaft))
       allocate (GR_Mp(J)%G(NR,dimaft))
       allocate (Gr_nn(J)%G(norbDyn(ephType),norbDyn(ephType)))
       allocate (Gr_1n(J)%G(NL,norbDyn(ephType)))
       allocate (Gr_Mn(J)%G(NR,norbDyn(ephType)))
    else
       allocate (GL_nn(dimbfr,dimbfr))
       allocate (GL_1N(NL,dimbfr))
    endif
    allocate (Gr_1M(NL,NR))


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
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 Ei                            : Energy grid point             !
!  *******************************************************************  !
  subroutine greenfunctions (Ei, ispin)

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Perform the left-to-right sweep.
    call LRsweep (Ei, ispin)

!   Perform the right-to-left sweep.
    call RLsweep (Ei, ispin)

!   Compute the required full Green's functions.
    call GFfull (Ei, ispin)

!   [test] Compute the entire Green's function of scattering region.
!!$    call GFtest (Ei, ispin)


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
!  integer NL                           : Number of left lead orbitals  !
!  complex(8) Sigma_L(NL,NL)            : Left-lead self-energy         !
!  integer ephIdx(ntypeunits+2)         : Unit index (those with e-ph)  !
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
                               S1unit, H1unit, Sunits, Hunits
    use idsrdr_leads,    only: NL, Sigma_L
    use idsrdr_ephcoupl, only: ephIdx
    use idsrdr_check,    only: CHECKzsytrf, CHECKzsytri

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    integer :: I, J, n, utype, dim, dimbfr, ephType
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

    if (nunits == 1) go to 203

!   Initialize variables.
    utype = ntypeunits + 1 ! current unit type (first unit)
    dim = unitdimensions(utype) ! current type dimension
    n = unitdimensions(ntypeunits) ! pristine dimension
    J = ephIdx(utype) + 1 ! first e-ph index

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
       ephType = ephIdx(utype) ! e-ph unit type

!      Store matrix if required.
       if (ephType /= 0) then
          GL_mm(J)%G = Gbfr
          GL_1m(J)%G = Gbfr_1m
          J = J + 1
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

!      ('Gaft = (E*S - H + V^T*Gbfr*V)^-1')
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

!      ('Gaft_1m = foo2 * Gaft')
       foo3 = Gaft(1:n,:)
       call zgemm ('N', 'N', NL, dim, n, (1.d0,0.d0), foo2, NL,         &
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
    ephType = ephIdx(ntypeunits+2) ! e-ph unit type
    if (ephType /= 0) then
       GL_mm(J)%G = Gbfr
       GL_1m(J)%G = Gbfr_1m
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

203 if (IOnode) write(6,'(a)') " ok!"


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
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Sigma_R(NR,NR)            : Right-lead self-energy        !
!  integer ephIdx(ntypeunits+2)         : Unit index (those with e-ph)  !
!  integer neph                         : Number of units with e-ph     !
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
                               S1unit, H1unit, Sunits, Hunits
    use idsrdr_leads,    only: NR, Sigma_R
    use idsrdr_ephcoupl, only: ephIdx, neph
    use idsrdr_check,    only: CHECKzsytrf, CHECKzsytri

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    integer :: I, J, n, utype, dim, ephType
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

    if (nunits == 1) go to 302

!   Initialize variables.
    utype = ntypeunits + 2 ! current unit type (last unit)
    dim = unitdimensions(utype) ! current type dimension
    n = unitdimensions(ntypeunits) ! pristine dimension

!   First e-ph unit from the right.
    if (ephIdx(utype) /= 0) then
       J = neph - 1
    else
       J = neph
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
       ephType = ephIdx(utype) ! e-ph unit type

!      Store matrix if required.
       if (ephType /= 0) then
          GR_pp(J)%G = Gbfr
          GR_Mp(J)%G = Gbfr_Mp
          J = J - 1
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

!      ('Gaft = (E*S - H + V^T*Gbfr*V)^-1')
       Gaft = (Ei-unitshift(utype))*Sunits(utype)%S                     &
              - Hunits(utype)%H(:,:,ispin)
       Gaft(dim-n+1:dim,dim-n+1:dim) = Gaft(dim-n+1:dim,dim-n+1:dim)    &
                                       - aux1
       call CHECKzsytrf (dim, 'L', Gaft, ipiv)
       call CHECKzsytri (dim, 'L', Gaft, ipiv)

!      (Re-)Allocate matrices.
       allocate (Gaft_Mp(NR,dim))
       allocate (foo3(n,dim))

!      ('foo2 = Gbfr_1m*V^dagger')
       foo1 = Gbfr_Mp(1:NR,1:n)
       call zgemm ('N', 'C', NR, n, n, (1.d0,0.d0), foo1, NR,           &
                   V, n, (0.d0,0.d0), foo2, NR)

!      ('Gaft_1m = foo2 * Gaft')
       foo3 = Gaft(dim-n+1:dim,:)
       call zgemm ('N', 'N', NR, dim, n, (1.d0,0.d0), foo2, NR,         &
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
    ephType = ephIdx(ntypeunits+1) ! e-ph unit type
    if (ephType /= 0) then
       GR_pp(J)%G = Gbfr
       GR_Mp(J)%G = Gbfr_Mp
    endif

!   Free memory.
    deallocate (V)
    deallocate (Gbfr)
    deallocate (Gbfr_Mp)
    deallocate (aux1)
    deallocate (aux2)
    deallocate (foo1)
    deallocate (foo2)

302 if (IOnode) write(6,'(a)') " ok!"


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
!  TYPE(green) Gr_Mn(neph)%G(NR,unitdimensions) :  [complex] G^r_{1,n}  !
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
    integer :: I, J, n, utype, dim, dimbfr, ephType
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
    ephType = ephIdx(utype) ! e-ph unit type
    J = 1

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

    if (ephType /= 0) then

!      Allocate auxiliary matrix.
       allocate (aux3(dim,dim))

!      ('aux3 = E*S - H - Sigma_L')
       aux3 = (0.d0,0.d0)
       aux3 = (Ei-unitshift(utype))*Sunits(utype)%S                     &
              - Hunits(utype)%H(:,:,ispin)
       aux3(1:NL,1:NL) = aux3(1:NL,1:NL) - Sigma_L

!      ('aux2 = V*GR_pp')
       aux1 = GR_pp(J)%G(1:n,1:n)
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
       Gr_nn(J)%G = aux3(idxF(ephType):idxL(ephType),                   &
                         idxF(ephType):idxL(ephType))

!      ('Gr_1n = aux3')
       Gr_1n(J)%G = aux3(1:NL,idxF(ephType):idxL(ephType))

!      Allocate auxiliary matrix.
       allocate (foo3(n,norbDyn(ephType)))

!      ('foo2R = GR_Mp*V^dagger')
       foo1R = GR_Mp(J)%G(1:NR,1:n)
       call zgemm ('N', 'C', NR, n, n, (1.d0,0.d0), foo1R, NR,          &
                   V, n, (0.d0,0.d0), foo2R, NR)

!      ('Gr_Mn = foo2R * Gr_nn')
       foo3 = aux3(dim-n+1:dim,idxF(ephType):idxL(ephType))
       call zgemm ('N', 'N', NR, norbDyn(ephType), n, (1.d0,0.d0),      &
                   foo2R, NR, foo3, n, (0.d0,0.d0), Gr_Mn(J)%G, NR)

!      Free memory.
       deallocate (foo3)
       deallocate (aux3)
       J = J + 1

    endif

!   Loop over the blocks from left to right.
    do I = 2,nunits+1

!      Assign auxiliary variables.
       dimbfr = dim ! last unit dimension
       utype = unit_type(I) ! current unit type
       dim = unitdimensions(utype) ! current type dimension
       ephType = ephIdx(utype) ! e-ph unit type

       if (ephType /= 0) then

!         Allocate auxiliary matrix.
          allocate (aux3(dim,dim))

!         ('aux3 = E*S - H')
          aux3 = (Ei-unitshift(utype))*Sunits(utype)%S                  &
                 - Hunits(utype)%H(:,:,ispin)

!         ('aux2 = GL_mm*V')
          aux1 = GL_mm(J)%G(dimbfr-n+1:dimbfr,dimbfr-n+1:dimbfr)
          call zsymm ('L', 'L', n, n, (1.d0,0.d0), aux1, n,             &
                      V, n, (0.d0,0.d0), aux2, n)

!         ('aux1 = V^dagger*aux2')
          call zgemm ('C', 'N', n, n, n, (1.d0,0.d0), V, n,             &
                      aux2, n, (0.d0,0.d0), aux1, n)

!         ('aux3 = aux3 - aux1')
          aux3(1:n,1:n) = aux3(1:n,1:n) - aux1

!         ('aux2 = V*GR_pp')
          aux1 = GR_pp(J)%G(1:n,1:n)
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
          Gr_nn(J)%G = aux3(idxF(ephType):idxL(ephType),                &
                            idxF(ephType):idxL(ephType))

!         Allocate auxiliary matrix.
          allocate (foo3(n,norbDyn(ephType)))

!         ('foo2L = GL_1m*V')
          foo1L = GL_1m(J)%G(1:NL,dimbfr-n+1:dimbfr)
          call zgemm ('N', 'N', NL, n, n, (1.d0,0.d0), foo1L, NL,       &
                      V, n, (0.d0,0.d0), foo2L, NL)

!         ('Gr_1n = foo2L * Gr_nn')
          foo3 = aux3(1:n,idxF(ephType):idxL(ephType))
          call zgemm ('N', 'N', NL, norbDyn(ephType), n, (1.d0,0.d0),   &
                      foo2L, NL, foo3, n, (0.d0,0.d0), Gr_1n(J)%G, NL)

!         ('foo2R = GR_Mp*V^dagger')
          foo1R = GR_Mp(J)%G(1:NR,1:n)
          call zgemm ('N', 'C', NR, n, n, (1.d0,0.d0), foo1R, NR,       &
                      V, n, (0.d0,0.d0), foo2R, NR)

!         ('Gr_Mn = foo2R * Gr_nn')
          foo3 = aux3(dim-n+1:dim,idxF(ephType):idxL(ephType))
          call zgemm ('N', 'N', NR, norbDyn(ephType), n, (1.d0,0.d0),   &
                      foo2R, NR, foo3, n, (0.d0,0.d0), Gr_Mn(J)%G, NR)

!         Free memory.
          deallocate (foo3)
          deallocate (aux3)

          J = J + 1

       endif
    enddo

!   Assign auxiliary variables.
    dimbfr = dim ! last unit dimension
    utype = ntypeunits + 2 ! current unit type (last unit)
    dim = unitdimensions(utype) ! current type dimension
    ephType = ephIdx(utype) ! e-ph unit type

    if (ephType /= 0) then

!      Allocate auxiliary matrix.
       allocate (aux3(dim,dim))

!      ('aux3 = E*S - H')
       aux3 = (Ei-unitshift(utype))*Sunits(utype)%S                     &
              - Hunits(utype)%H(:,:,ispin)
       aux3(dim-NR+1:dim,dim-NR+1:dim) =                                &
                                        aux3(dim-NR+1:dim,dim-NR+1:dim) &
                                        - Sigma_R

!      ('aux2 = GL_mm*V')
       aux1 = GL_mm(J)%G(dimbfr-n+1:dimbfr,dimbfr-n+1:dimbfr)
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
       Gr_nn(J)%G = aux3(idxF(ephType):idxL(ephType),                   &
                         idxF(ephType):idxL(ephType))

!      ('Gr_Mn = aux3')
       Gr_Mn(J)%G = aux3(dim-NR+1:dim,idxF(ephType):idxL(ephType))

!      Allocate auxiliary matrix.
       allocate (foo3(n,norbDyn(ephType)))

!      ('foo2L = GL_1m*V')
       foo1L = GL_1m(J)%G(1:NL,dimbfr-n+1:dimbfr)
       call zgemm ('N', 'N', NL, n, n, (1.d0,0.d0), foo1L, NL,          &
                   V, n, (0.d0,0.d0), foo2L, NL)

!      ('Gr_1n = foo2L * Gr_nn')
       foo3 = aux3(1:n,idxF(ephType):idxL(ephType))
       call zgemm ('N', 'N', NL, norbDyn(ephType), n, (1.d0,0.d0),      &
                   foo2L, NL, foo3, n, (0.d0,0.d0), Gr_1n(J)%G, NL)
       Gr_1M = Gr_1n(J)%G

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

!      ('aux2 = GL_mm*V')
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

!      ('Gr_1n = foo2L * aux3')
       foo3 = aux3(1:n,dim-NR+1:dim)
       call zgemm ('N', 'N', NL, NR, n, (1.d0,0.d0), foo2L, NL,        &
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
!                                GFtest                                 !
!  *******************************************************************  !
!  Description: build the full hamiltonian of the scattering region     !
!  and invert it to obtain the entire Green's function in order to      !
!  compare with the Green's function obtained with the recursive        !
!  method.                                                              !
!                                                                       !
!  ATENTION: the writting part only works in serial mode!
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
    use idsrdr_check,    only: CHECKzsytrf, CHECKzsytri

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    integer :: i, j, k, w, utype, dim, dimTot, dimCpl, idxAnt, ephType
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
    ephType = ephIdx(utype) ! e-ph unit type

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
    call CHECKzsytrf (dimTot, 'L', Gtot, ipiv)
    call CHECKzsytri (dimTot, 'L', Gtot, ipiv)

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

!   First unit.
    w = 1
    if (ephType /= 0) then

       do j = idxF(ephType),idxL(ephType)
          do i = j,idxL(ephType)
             write (3102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Gtot(i,j)),                                     &
                  DREAL(Gr_nn(w)%G(i-idxF(ephType)+1,                   &
                                   j-idxF(ephType)+1)),                 &
                  DIMAG(Gtot(i,j)),                                     &
                  DIMAG(Gr_nn(w)%G(i-idxF(ephType)+1,                   &
                                   j-idxF(ephType)+1))
             write (2102,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Gtot(i,j)) -                                    &
                  DREAL(Gr_nn(w)%G(i-idxF(ephType)+1,                   &
                                   j-idxF(ephType)+1)),                 &
                  DIMAG(Gtot(i,j)) -                                    &
                  DIMAG(Gr_nn(w)%G(i-idxF(ephType)+1,                   &
                                   j-idxF(ephType)+1))
          enddo
       enddo
       w = w + 1

    endif

!   Units in the middle.
    idxAnt = unitdimensions(ntypeunits+1) ! last unit dimensions
    do k = 2,nunits+1

       utype = unit_type(k) ! current unit type
       dim = unitdimensions(utype) ! current unit dimensions
       ephType = ephIdx(utype) ! e-ph unit type

       if (ephType /= 0) then

          do j = idxF(ephType),idxL(ephType)
             do i = j,idxL(ephType)
                write (3102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')        &
                     DREAL(Gtot(idxAnt+i,idxAnt+j)),                    &
                     DREAL(Gr_nn(w)%G(i-idxF(ephType)+1,                &
                                      j-idxF(ephType)+1)),              &
                     DIMAG(Gtot(idxAnt+i,idxAnt+j)),                    &
                     DIMAG(Gr_nn(w)%G(i-idxF(ephType)+1,                &
                                      j-idxF(ephType)+1))
                write (2102,'(e17.8e3,e17.8e3)')                        &
                     DREAL(Gtot(idxAnt+i,idxAnt+j)) -                   &
                     DREAL(Gr_nn(w)%G(i-idxF(ephType)+1,                &
                                      j-idxF(ephType)+1)),              &
                     DIMAG(Gtot(idxAnt+i,idxAnt+j)) -                   &
                     DIMAG(Gr_nn(w)%G(i-idxF(ephType)+1,                &
                                      j-idxF(ephType)+1))
             enddo
          enddo
          w = w + 1

       endif

       idxAnt = idxAnt + dim

    enddo

!   Last unit.
    ephType = ephIdx(ntypeunits+2) ! e-ph unit type
    if (ephType /= 0) then

       do j = idxF(ephType),idxL(ephType)
          do i = j,idxL(ephType)
             write (3102,'(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')           &
                  DREAL(Gtot(idxAnt+i,idxAnt+j)),                       &
                  DREAL(Gr_nn(w)%G(i-idxF(ephType)+1,                   &
                                   j-idxF(ephType)+1)),                 &
                  DIMAG(Gtot(idxAnt+i,idxAnt+j)),                       &
                  DIMAG(Gr_nn(w)%G(i-idxF(ephType)+1,                   &
                                   j-idxF(ephType)+1))
             write (2102,'(e17.8e3,e17.8e3)')                           &
                  DREAL(Gtot(idxAnt+i,idxAnt+j)) -                      &
                  DREAL(Gr_nn(w)%G(i-idxF(ephType)+1,                   &
                                   j-idxF(ephType)+1)),                 &
                  DIMAG(Gtot(idxAnt+i,idxAnt+j)) -                      &
                  DIMAG(Gr_nn(w)%G(i-idxF(ephType)+1,                   &
                                   j-idxF(ephType)+1))
          enddo
       enddo

    endif

!   Free memory.
    deallocate (Gtot)


  end subroutine GFtest


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
!  integer neph                : Number of units with e-ph interaction  !
!  *******************************************************************  !
  subroutine freegreen

!
!   Modules
!
    use idsrdr_ephcoupl, only: neph

!   Local variables.
    integer :: I

!   First deallocates pointed matrices.
    do I = 1,neph
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
    if (allocated(GL_nn)) then
       deallocate (GL_nn)
       deallocate (GL_1N)
    endif
    deallocate (Gr_1M)


  end subroutine freegreen


!  *******************************************************************  !


END MODULE idsrdr_green
