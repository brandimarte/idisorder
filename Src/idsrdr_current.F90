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
!                         MODULE idsrdr_current                         !
!  *******************************************************************  !
!  Description: compute the eletronic current.                          !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *******************************************************************  !

MODULE idsrdr_current

!
!   Modules
!
  use parallel,        only: 
  use idsrdr_leads,    only: 
  use idsrdr_green,    only: 

  implicit none
  
  PUBLIC  :: current
  PRIVATE :: elastic, inelastic


CONTAINS


!  *******************************************************************  !
!                                current                                !
!  *******************************************************************  !
!  Description: main subroutine for computing the current.              !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Sigma_L(NL,NL)           : Left-lead self-energy          !
!  complex(8) Sigma_R(NR,NR)           : Right-lead self-energy         !
!  ****************************** INPUT ******************************  !
!  real*8 Ei                           : Energy grid point              !
!  integer ispin                       : Spin component index           !
!  *******************************************************************  !
  subroutine current (Ei, ispin)

!
!   Modules
!
    use idsrdr_leads,    only: NL, NR, Sigma_L, Sigma_R

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.
    integer :: i, j
    complex(8), parameter :: zi = (0.D0,1.D0) ! complex i
    complex(8), allocatable, dimension (:,:) :: Gamma_L
    complex(8), allocatable, dimension (:,:) :: Gamma_R

!   Allocate coupling and auxiliary matrices.
    allocate (Gamma_L(NL,NL))
    allocate (Gamma_R(NR,NR))

!   Sets the lead's coupling matricesx (triangular inferior part).
    if (NL >= NR) then
       Do j = 1,NR
          do i = j,NR
             Gamma_L(i,j) = zi * (Sigma_L(i,j) - DCONJG(Sigma_L(j,i)))
             Gamma_R(i,j) = zi * (Sigma_R(i,j) - DCONJG(Sigma_R(j,i)))
          enddo
          do i = NR+1,NL
             Gamma_L(i,j) = zi * (Sigma_L(i,j) - DCONJG(Sigma_L(j,i)))
          enddo
       Enddo
       Do j = NR+1,NL
          do i = j,NL
             Gamma_L(i,j) = zi * (Sigma_L(i,j) - DCONJG(Sigma_L(j,i)))
          enddo
       Enddo
    else ! NL < NR
       Do j = 1,NL
          do i = j,NL
             Gamma_L(i,j) = zi * (Sigma_L(i,j) - DCONJG(Sigma_L(j,i)))
             Gamma_R(i,j) = zi * (Sigma_R(i,j) - DCONJG(Sigma_R(j,i)))
          enddo
          do i = NL+1,NR
             Gamma_R(i,j) = zi * (Sigma_R(i,j) - DCONJG(Sigma_R(j,i)))
          enddo
       Enddo
       Do j = NL+1,NR
          do i = j,NR
             Gamma_R(i,j) = zi * (Sigma_R(i,j) - DCONJG(Sigma_R(j,i)))
          enddo
       Enddo
    endif

!   Compute elastic contribution.
    call elastic (Ei, ispin, NL, Gamma_L, NR, Gamma_R)

!   Compute inelastic contribution.
    call inelastic (Ei, ispin, NL, Gamma_L, NR, Gamma_R)

!   Free memory.
    deallocate (Gamma_L)
    deallocate (Gamma_R)


  end subroutine current


!  *******************************************************************  !
!                                elastic                                !
!  *******************************************************************  !
!  Description: compute the elastic part from current expression.       !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode                      : True if it is the I/O node     !
!  ****************************** INPUT ******************************  !
!  real*8 Ei                           : Energy grid point              !
!  integer ispin                       : Spin component index           !
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Gamma_L(NL,NL)           : Left-lead coupling matrix      !
!  complex(8) Gamma_R(NR,NR)           : Right-lead coupling matrix     !
!  ***************************** OUTPUT ******************************  !
!  *******************************************************************  !
  subroutine elastic (Ei, ispin, NL, Gamma_L, NR, Gamma_R)

!
!   Modules
!
    use parallel,        only: IOnode

!   Input variables.
    integer, intent(in) :: ispin, NL, NR
    real(8), intent(in) :: Ei
    complex(8), dimension (NL,NL), intent(in) :: Gamma_L
    complex(8), dimension (NR,NR), intent(in) :: Gamma_R

!   Local variables.

    if (IOnode) write (6,'(a)', advance='no')                           &
            '      computing elastic current... '


    if (IOnode) write(6,'(a)') " ok!"


  end subroutine elastic


!  *******************************************************************  !
!                               inelastic                               !
!  *******************************************************************  !
!  Description: compute the inelastic part from current expression.     !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer neph                : Number of units with e-ph interaction  !
!  integer nModes(neph)        : Number of vibrational modes            !
!  integer norbDyn(neph)       : Number of orbitals from dynamic atoms  !
!  TYPE(ephCplng) Meph(neph)%M(norbDyn,norbDyn,nspin,nModes) :          !
!                                   [complex*8] E-ph coupling matrices  !
!  TYPE(green) Gr_Mn(neph)%G(NR,unitdimensions) :  [complex] G^r_{1,n}  !
!  ****************************** INPUT ******************************  !
!  real*8 Ei                           : Energy grid point              !
!  integer ispin                       : Spin component index           !
!  integer NL                          : Number of left lead orbitals   !
!  integer NR                          : Number of right lead orbitals  !
!  complex(8) Gamma_L(NL,NL)           : Left-lead coupling matrix      !
!  complex(8) Gamma_R(NR,NR)           : Right-lead coupling matrix     !
!  ***************************** OUTPUT ******************************  !
!  *******************************************************************  !
  subroutine inelastic (Ei, ispin, NL, Gamma_L, NR, Gamma_R)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_ephcoupl, only: neph, nModes, norbDyn, Meph
    use idsrdr_green,    only: Gr_Mn

!   Input variables.
    integer, intent(in) :: ispin, NL, NR
    real(8), intent(in) :: Ei
    complex(8), dimension (NL,NL), intent(in) :: Gamma_L
    complex(8), dimension (NR,NR), intent(in) :: Gamma_R

!   Local variables.
    integer :: j, w
    complex(8), dimension(:,:), allocatable :: GrT_Mn, Aux
    external :: zsymm, zgemm

    if (IOnode) write (6,'(a)', advance='no')                           &
            '      computing inelastic current... '

    do j = 1,neph ! over unit with e-ph

!      Allocate auxiliary matrices.
       allocate (GrT_Mn(norbDyn(j),NR))
       allocate (Aux(norbDyn(j),NR))

       GrT_Mn = TRANSPOSE(Gr_Mn(j)%G)

       do w = 1,nModes(j) ! over phonon modes

!         ('Aux = Meph*GrT_Mn')
          call zsymm ('L', 'L', norbDyn(j), NR, (1.d0,0.d0),            &
                      Meph(j)%M(:,:,ispin,w), norbDyn(j), GrT_Mn,       &
                      norbDyn(j), (0.d0,0.d0), Aux, norbDyn(j))

          
       enddo

!      Free memory.
       deallocate (GrT_Mn)
       deallocate (Aux)

    enddo


    if (IOnode) write(6,'(a)') " ok!"


  end subroutine inelastic


!  *******************************************************************  !
!                              freecurrent                              !
!  *******************************************************************  !
!  Description: free allocated vectors.                                 !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer neph                : Number of units with e-ph interaction  !
!  *******************************************************************  !
!!$  subroutine freecurrent
!!$
!!$!
!!$!   Modules
!!$!
!!$    use idsrdr_ephcoupl, only: neph
!!$
!!$!   Local variables.
!!$    integer :: I
!!$
!!$!   First deallocates pointed matrices.
!!$    do I = 1,neph
!!$       deallocate (GL_mm(I)%G)
!!$       deallocate (GL_1m(I)%G)
!!$       deallocate (GR_pp(I)%G)
!!$       deallocate (GR_Mp(I)%G)
!!$       deallocate (Gr_nn(I)%G)
!!$       deallocate (Gr_1n(I)%G)
!!$       deallocate (Gr_Mn(I)%G)
!!$    enddo
!!$    deallocate (GL_mm)
!!$    deallocate (GL_1m)
!!$    deallocate (GR_pp)
!!$    deallocate (GR_Mp)
!!$    deallocate (Gr_nn)
!!$    deallocate (Gr_1n)
!!$    deallocate (Gr_Mn)
!!$
!!$
!!$  end subroutine freecurrent


!  *******************************************************************  !


END MODULE idsrdr_current
