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
  PRIVATE :: LRsweep , RLsweep, GFfull

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
!  integer nunits                       : Total number of units         !
!  integer unit_type(nunits+2)          : Units types                   !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  integer ephIndic(ntypeunits+2)       : E-ph interaction indicator    !
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  integer neph                        : Number of units with e-ph      !
!                                        interaction                    !
!  integer norbDyn(neph)               : Number of orbitals from        !
!                                        dynamic atoms                  !
!  *******************************************************************  !
  subroutine greeninit

!
!   Modules
!
    use idsrdr_options,  only: nunits
    use idsrdr_units,    only: unit_type, unitdimensions, ephIndic
    use idsrdr_leads,    only: NL, NR
    use idsrdr_ephcoupl, only: neph, norbDyn

!   Local variables.
    integer :: I, J, dim

!   Allocate Green's functions pointer array.
    allocate (GL_mm(neph))
    allocate (GL_1m(neph))
    allocate (GR_pp(neph))
    allocate (GR_Mp(neph))
    allocate (Gr_nn(neph))
    allocate (Gr_1n(neph))
    allocate (Gr_Mn(neph))

!   Allocate Green's functions matrices.
    J = 1
    do I = 2,nunits+1
       if (ephIndic(unit_type(I)) == 1) then
          dim = unitdimensions(unit_type(I)) ! current unit dimension
          allocate (GL_mm(J)%G(dim,dim))
          allocate (GL_1m(J)%G(NL,dim))
          allocate (GR_pp(J)%G(dim,dim))
          allocate (GR_Mp(J)%G(NR,dim))
          dim = norbDyn(J) ! dynamic atoms dimension
          allocate (Gr_nn(J)%G(dim,dim))
          allocate (Gr_1n(J)%G(NL,dim))
          allocate (Gr_Mn(J)%G(NR,dim))
          J = J + 1
       endif
    enddo


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
!  integer ephIndic(ntypeunits+2)       : E-ph interaction indicator    !
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
                               ephIndic, S1unit, H1unit, Sunits, Hunits
    use idsrdr_leads,    only: NL, Sigma_L
    use idsrdr_check,    only: CHECKzsytrf, CHECKzsytri

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    integer :: I, n, J, utype, dim
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
    utype = ntypeunits+1 ! current unit type (left lead)
    dim = unitdimensions(utype) ! current type dimension
    n = unitdimensions(ntypeunits) ! pristine dimension
    J = 1 ! Green's function store indexing

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
           - Hunits(utype)%H(:,:,ispin) - Sigma_L
    call CHECKzsytrf (dim, 'L', Gbfr, ipiv)
    call CHECKzsytri (dim, 'L', Gbfr, ipiv)
    Gbfr_1m = Gbfr

!   Free memory.
    deallocate (ipiv)

!   Store matrix if required.
    if (ephIndic(unit_type(2)) == 1) then
       GL_mm(J)%G = Gbfr
       GL_1m(J)%G = Gbfr_1m
       J = J + 1
    endif

!   Loop over the blocks from left to right.
    do I = 2,nunits

!      ('aux2 = Gbfr*V')
       aux1 = Gbfr(dim-n+1:dim,dim-n+1:dim)
       call zsymm ('L', 'L', n, n, (1.d0,0.d0), aux1, n,                &
                   V, n, (0.d0,0.d0), aux2, n)

!      ('aux1 = V^dagger*aux2')
       call zgemm ('C', 'N', n, n, n, (1.d0,0.d0), V, n,                &
                   aux2, n, (0.d0,0.d0), aux1, n)

!      Assign auxiliary variables.
       utype = unit_type(I) ! current unit type
       dim = unitdimensions(utype) ! current type dimension

!      (Re-)Allocate matrices and array.
       allocate (Gaft(dim,dim))
       allocate (ipiv(dim))

!      ('Gaft = (E*S - H + V^T*Gbfr*V)^-1')
       Gaft = (Ei-unitshift(utype))*Sunits(utype)%S                     &
              - Hunits(utype)%H(:,:,ispin)
       Gaft(dim-n+1:dim,dim-n+1:dim) = Gaft(dim-n+1:dim,dim-n+1:dim)    &
                                       + aux1
       call CHECKzsytrf (dim, 'L', Gaft, ipiv)
       call CHECKzsytri (dim, 'L', Gaft, ipiv)

!      (Re-)Allocate matrices.
       allocate (Gaft_1m(NL,dim))
       allocate (foo3(n,dim))

!      ('foo2 = Gbfr_1m*V')
       foo1 = Gbfr_1m(1:NL,dim-n+1:dim)
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

!      Store matrix if required.
       if (ephIndic(unit_type(I+1)) == 1) then
          GL_mm(J)%G = Gbfr
          GL_1m(J)%G = Gbfr_1m
          J = J + 1
       endif

!      (Re-)Free memory.
       deallocate (Gaft)
       deallocate (Gaft_1m)
       deallocate (foo3)
       deallocate (ipiv)

    enddo

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
!  integer ephIndic(ntypeunits+2)       : E-ph interaction indicator    !
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
!  integer neph                         : Number of units with e-ph     !
!                                         interaction                   !
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
                               ephIndic, S1unit, H1unit, Sunits, Hunits
    use idsrdr_leads,    only: NR, Sigma_R
    use idsrdr_ephcoupl, only: neph
    use idsrdr_check,    only: CHECKzsytrf, CHECKzsytri

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    integer :: I, n, J, utype, dim
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
    utype = ntypeunits+2 ! current unit type (right lead)
    dim = unitdimensions(utype) ! current type dimension
    n = unitdimensions(ntypeunits) ! pristine dimension
    J = neph ! Green's function store indexing

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
           - Hunits(utype)%H(:,:,ispin) - Sigma_R
    call CHECKzsytrf (dim, 'L', Gbfr, ipiv)
    call CHECKzsytri (dim, 'L', Gbfr, ipiv)
    Gbfr_Mp = Gbfr

!   Free memory.
    deallocate (ipiv)

!   Loop over the blocks from left to right.
    do I = nunits+1,3,-1

!      Store matrix if required.
       if (ephIndic(unit_type(I)) == 1) then
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

!      Assign auxiliary variables.
       utype = unit_type(I) ! current unit type
       dim = unitdimensions(utype) ! current type dimension

!      (Re-)Allocate matrices and array.
       allocate (Gaft(dim,dim))
       allocate (ipiv(dim))

!      ('Gaft = (E*S - H + V^T*Gbfr*V)^-1')
       Gaft = (Ei-unitshift(utype))*Sunits(utype)%S                     &
              - Hunits(utype)%H(:,:,ispin)
       Gaft(dim-n+1:dim,dim-n+1:dim) = Gaft(dim-n+1:dim,dim-n+1:dim)    &
                                       + aux1
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
    if (ephIndic(unit_type(I)) == 1) then
       GR_pp(J)%G = Gbfr
       GR_Mp(J)%G = Gbfr_Mp
       J = J - 1
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
!  integer ephIndic(ntypeunits+2)       : E-ph interaction indicator    !
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
!  integer neph                        : Number of units with e-ph      !
!                                        interaction                    !
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
                               ephIndic, S1unit, H1unit, Sunits, Hunits
    use idsrdr_leads,    only: NL, NR
    use idsrdr_ephcoupl, only: neph, norbDyn, idxF, idxL
    use idsrdr_check,    only: CHECKzsytrf, CHECKzsytri

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    integer :: I, n, J, utype, dim
    integer, allocatable, dimension (:) :: ipiv
    complex(8), allocatable, dimension (:,:) :: V ! pristine coupling
    complex(8), allocatable, dimension (:,:) :: aux1, aux2, aux3
    complex(8), allocatable, dimension (:,:) :: foo1, foo2, foo3
    external :: zsymm, zgemm

    if (IOnode) write (6,'(a)', advance='no')                           &
            '      computing full Greens functions... '

!   Initialize variables.
    utype = ntypeunits+1 ! current unit type (left lead)
    dim = unitdimensions(utype) ! current type dimension
    n = unitdimensions(ntypeunits) ! pristine dimension
    J = 1 ! Green's function store indexing

!   Allocate matrices.
    allocate (V(n,n))
    allocate (aux1(n,n))
    allocate (aux2(n,n))
    allocate (foo1(NL,n))
    allocate (foo2(NL,n))

!   Coupling matrix (pristine).
    V = (Ei-unitshift(ntypeunits))*S1unit - H1unit(:,:,ispin)

!   Loop over the blocks from left to right.
    do I = 2,nunits+1
       if (ephIndic(unit_type(I)) == 1) then

!         Assign auxiliary variables.
          utype = unit_type(I) ! current unit type
          dim = unitdimensions(utype) ! current type dimension

!         Allocate auxiliary matrix.
          allocate (aux3(dim,dim))

!         ('aux3 = E*S - H + V^T*Gbfr*V')
          aux3 = (Ei-unitshift(utype))*Sunits(utype)%S                  &
                 - Hunits(utype)%H(:,:,ispin)

!         ('aux2 = GL_mm*V')
          aux1 = GL_mm(J)%G(dim-n+1:dim,dim-n+1:dim)
          call zsymm ('L', 'L', n, n, (1.d0,0.d0), aux1, n,             &
                      V, n, (0.d0,0.d0), aux2, n)

!         ('aux1 = V^dagger*aux2')
          call zgemm ('C', 'N', n, n, n, (1.d0,0.d0), V, n,             &
               aux2, n, (0.d0,0.d0), aux1, n)

!         ('aux3 = aux3 + aux1')
          aux3(1:n,1:n) = aux3(1:n,1:n) + aux1

!         ('aux2 = V*GR_pp')
          aux1 = GR_pp(J)%G(1:n,1:n)
          call zsymm ('R', 'L', n, n, (1.d0,0.d0), aux1, n,             &
                      V, n, (0.d0,0.d0), aux2, n)

!         ('aux1 = aux2*V^dagger')
          call zgemm ('N', 'C', n, n, n, (1.d0,0.d0), aux2, n,          &
                      V, n, (0.d0,0.d0), aux1, n)

!         ('aux3 = aux3 + aux1')
          aux3(dim-n+1:dim,dim-n+1:dim) = aux3(dim-n+1:dim,dim-n+1:dim) &
                                          + aux1

!         ('aux3 = aux3^-1')
          allocate (ipiv(dim))
          call CHECKzsytrf (dim, 'L', aux3, ipiv)
          call CHECKzsytri (dim, 'L', aux3, ipiv)
          deallocate (ipiv)

!         ('Gr_nn = aux3')
          Gr_nn(J)%G = aux3(idxF(J):idxL(J),idxF(J):idxL(J))

!         Allocate auxiliary matrix.
          allocate (foo3(n,norbDyn(J)))

!         ('foo2 = GL_1m*V')
          foo1 = GL_1m(J)%G(1:NL,dim-n+1:dim)
          call zgemm ('N', 'N', NL, n, n, (1.d0,0.d0), foo1, NL,        &
                      V, n, (0.d0,0.d0), foo2, NL)

!         ('Gr_1n = foo2 * Gr_nn')
          foo3 = aux3(1:n,idxF(J):idxL(J))
          call zgemm ('N', 'N', NL, norbDyn(J), n, (1.d0,0.d0),         &
                      foo2, NL, foo3, n, (0.d0,0.d0), Gr_1n(J)%G, NL)

!         ('foo2 = GR_Mp*V^dagger')
          foo1 = GR_Mp(J)%G(1:NR,1:n)
          call zgemm ('N', 'C', NR, n, n, (1.d0,0.d0), foo1, NR,        &
                      V, n, (0.d0,0.d0), foo2, NR)

!         ('Gr_Mn = foo2 * Gr_nn')
          foo3 = aux3(dim-n+1:dim,idxF(J):idxL(J))
          call zgemm ('N', 'N', NR, norbDyn(J), n, (1.d0,0.d0),         &
                      foo2, NR, foo3, n, (0.d0,0.d0), Gr_Mn(J)%G, NR)

!         Free memory.
          deallocate (foo3)
          deallocate (aux3)

          J = J + 1

       endif
    enddo

!   Free memory.
    deallocate (V)
    deallocate (aux1)
    deallocate (aux2)
    deallocate (foo1)
    deallocate (foo2)

    if (IOnode) write(6,'(a)') " ok!"


  end subroutine GFfull


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


  end subroutine freegreen


!  *******************************************************************  !


END MODULE idsrdr_green
