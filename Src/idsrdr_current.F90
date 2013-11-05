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
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 Ei                            : Energy grid point             !
!  *******************************************************************  !
  subroutine current (Ei, ispin)

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Compute elastic contribution.
    call elastic (Ei, ispin)

!   Compute inelastic contribution.
    call inelastic (Ei, ispin)


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
!  logical IOnode                       : True if it is the I/O node    !
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 Ei                            : Energy grid point             !
!  ***************************** OUTPUT ******************************  !
!  *******************************************************************  !
  subroutine elastic (Ei, ispin)

!
!   Modules
!
    use parallel,        only: IOnode

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

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
!  logical IOnode                       : True if it is the I/O node    !
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 Ei                            : Energy grid point             !
!  ***************************** OUTPUT ******************************  !
!  *******************************************************************  !
  subroutine inelastic (Ei, ispin)

!
!   Modules
!
    use parallel,        only: IOnode

!   Input variables.
    integer, intent(in) :: ispin
    real(8), intent(in) :: Ei

!   Local variables.

    if (IOnode) write (6,'(a)', advance='no')                           &
            '      computing inelastic current... '



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
