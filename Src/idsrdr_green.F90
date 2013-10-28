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

  implicit none
  
  PUBLIC ! default is public
!!$  PRIVATE :: LRsweep, RLsweep


CONTAINS


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
!  logical IOnode                       : True if it is the I/O node    !
!  integer ntypeunits                   : Number of unit types          !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  real*8 unitshift(ntypeunits+2)       : Units shift                   !
!  TYPE(unitS) Sunits(ntypeunits+2)%S(unitdimensions,unitdimensions) :  !
!                                         [real*8] Units overlap        !
!  TYPE(unitH)                                                          !
!        Hunits(ntypeunits+2)%H(unitdimensions,unitdimensions,nspin) :  !
!                                         [real*8] Units hamiltonian    !
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 Ei                            : Energy grid point             !
!  ***************************** OUTPUT ******************************  !
!  integer N1                           : Leads total number of         !
!                                         orbitals (left + right)       !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  real*8 theta(NDeffects+1)            :                               !
!  real*8 unitshift(ntypeunits+2)       : Units shift                   !
!  *******************************************************************  !
  subroutine greenfunctions (Ei, ispin)

!
!   Modules
!
    use idsrdr_options,  only: ntypeunits
    use idsrdr_units,    only: unitdimensions, unitshift, Sunits, Hunits

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    real(8), allocatable, dimension (:,:) :: V0 ! pristine coupling

!   Allocate arrays.
    allocate (V0(unitdimensions(ntypeunits),unitdimensions(ntypeunits)))

!   Initialize arrays.
    V0 = (Ei-unitshift(ntypeunits))*Sunits(ntypeunits)%S                &
         - Hunits(ntypeunits)%H(:,:,ispin)

!   Free memory.
    deallocate (V0)


  end subroutine greenfunctions


!  *******************************************************************  !
!                                LRsweep                                !
!  *******************************************************************  !
!  Description: Left-to-right sweep for computing 'G^L_{n,n}',          !
!  'G^L_{1,n}' and 'G^L_{N,n}'.                                         !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode               : True if it is the I/O node            !
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
!  *******************************************************************  !
!!$  subroutine LRsweep
!!$
!!$!
!!$!   Modules
!!$!
!!$
!!$
!!$  end subroutine LRsweep


!  *******************************************************************  !
!                                RLsweep                                !
!  *******************************************************************  !
!  Description:                                                         !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
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
!  *******************************************************************  !
!!$  subroutine RLsweep
!!$
!!$!
!!$!   Modules
!!$!
!!$
!!$
!!$  end subroutine RLsweep


!  *******************************************************************  !
!                                gfHead                                 !
!  *******************************************************************  !
!  Description: Prints an "begin" message.                              !
!  Description: free allocated vectors.                                 !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *******************************************************************  !
  subroutine gfHead

!
!   Modules
!
    use parallel,        only: IOnode

    if (IOnode) write (6,'(/,25("*"),a,26("*"),/)')                     &
         ' Computing Greens functions '

  end subroutine gfHead


!  *******************************************************************  !
!                                gfTail                                 !
!  *******************************************************************  !
!  Description: Prints an "finish" message.                             !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *******************************************************************  !
  subroutine gfTail

!
!   Modules
!
    use parallel,        only: IOnode

    if (IOnode) write (6,'(/,2a)') 'greenfunctions: ', repeat('*', 63)


  end subroutine gfTail


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
!  *******************************************************************  !
!!$  subroutine freegreen
!!$
!!$!
!!$!   Modules
!!$!
!!$    use idsrdr_options,  only: ntypeunits
!!$
!!$!   Local variables.
!!$    integer :: I
!!$
!!$!   Free memory.
!!$    deallocate (theta)
!!$    deallocate (unitshift)
!!$    deallocate (S1unit)
!!$    deallocate (H1unit)
!!$!   First deallocates pointed matrices.
!!$    do I = 1,ntypeunits
!!$       deallocate (Sunits(I)%S)
!!$       deallocate (Hunits(I)%H)
!!$    enddo
!!$    deallocate (Sunits)
!!$    deallocate (Hunits)
!!$    deallocate (unit_type)
!!$    deallocate (unitdimensions)
!!$
!!$
!!$  end subroutine freegreen


!  *******************************************************************  !


END MODULE idsrdr_green
