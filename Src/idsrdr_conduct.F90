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
!                         MODULE idsrdr_conduct                         !
!  *******************************************************************  !
!  Description: compute the differential conductance (dIdV) and its     !
!  derivative (d2IdV2).                                                 !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *******************************************************************  !

MODULE idsrdr_conduct

!
!   Modules
!
  use parallel,        only: 
  use idsrdr_options,  only: 
  use idsrdr_engrid,   only:
  use idsrdr_current,  only: 

  implicit none
  
  PUBLIC  :: conductinit, conduct, alldIdV, dIdV, d2IdV2, freedIdV,     &
             sumAlldIdV
  PRIVATE :: computedIdV, computed2IdV2

! Type for storing calculated differential conductance.
  TYPE alldIdV
     sequence
     real(8) :: el    ! elastic part
     real(8) :: symm  ! inelastic symmetric
     real(8) :: asymm ! inelastic asymmetric
  END TYPE alldIdV

! Differential conductance (dIdV) and its derivative (d2IdV2).
  TYPE(alldIdV), allocatable, dimension (:,:,:) :: dIdV, d2IdV2


CONTAINS


!  *******************************************************************  !
!                              conductinit                              !
!  *******************************************************************  !
!  Description: allocate array of type 'alldIdV' for storing            !
!  calculated differential conductance.                                 !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer nspin               : Number of spin components              !
!  integer NIVP                : Number of bias potential points        !
!  integer NTenerg_div         : Number of energy grid points per node  !
!  *******************************************************************  !
  subroutine conductinit

!
!   Modules
!
    use idsrdr_options,  only: nspin, NIVP
    use idsrdr_engrid,   only: NTenerg_div

    if (NIVP < 3) return

!   Allocate and initializes current array.
    allocate (dIdV(NTenerg_div,nspin,NIVP))
    dIdV%el = 0.d0
    dIdV%symm = 0.d0
    dIdV%asymm = 0.d0
    allocate (d2IdV2(NTenerg_div,nspin,NIVP))
    d2IdV2%el = 0.d0
    d2IdV2%symm = 0.d0
    d2IdV2%asymm = 0.d0


  end subroutine conductinit


!  *******************************************************************  !
!                                conduct                                !
!  *******************************************************************  !
!  Description: interface subroutine for computing the differential     !
!  conductance.                                                         !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer NIVP                : Number of bias potential points        !
!  *******************************************************************  !
  subroutine conduct

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_options,  only: NIVP

    if (NIVP < 3) return

    if (IOnode) write (6,'(a)', advance='no')                           &
         '      computing differential conductance... '

!   Compute the differential conductance.
    call computedIdV

!   Compute the derivative of the differential conductance.
    call computed2IdV2

    if (IOnode) write(6,'(a)') " ok!"


  end subroutine conduct


!  *******************************************************************  !
!                              computedIdV                              !
!  *******************************************************************  !
!  Description: compute the differential conductance by applying        !
!  finite differences at calculated currents.                           !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer nspin               : Number of spin components              !
!  integer NIVP                : Number of bias potential points        !
!  real*8 dV                   : Bias potential step                    !
!  integer NTenerg_div         : Number of energy grid points per node  !
!  TYPE(calcCurr) allcurr(NTenerg_div,nspin,NIVP) : [real] calculated   !
!                                                   currents            !
!  *******************************************************************  !
  subroutine computedIdV

!
!   Modules
!
    use idsrdr_options,  only: nspin, NIVP, dV
    use idsrdr_engrid,   only: NTenerg_div
    use idsrdr_current,  only: allcurr

!   Local variables.
    integer :: e, s, v

    do e = 1,NTenerg_div
       do s = 1,nspin

!         First point.
          dIdV(e,s,1)%el = (allcurr(e,s,2)%el - allcurr(e,s,1)%el) / dV
          dIdV(e,s,1)%symm = (allcurr(e,s,2)%isymm -                    &
               allcurr(e,s,1)%isymm) / dV
          dIdV(e,s,1)%asymm = (allcurr(e,s,2)%iasymm -                  &
               allcurr(e,s,1)%iasymm) / dV

!         Intermediate points.
          do v = 2,NIVP-1
             dIdV(e,s,v)%el = (allcurr(e,s,v+1)%el -                    &
                  allcurr(e,s,v-1)%el) / (2.d0 * dV)
             dIdV(e,s,v)%symm = (allcurr(e,s,v+1)%isymm -               &
                  allcurr(e,s,v-1)%isymm) / (2.d0 * dV)
             dIdV(e,s,v)%asymm = (allcurr(e,s,v+1)%iasymm -             &
                  allcurr(e,s,v-1)%iasymm) / (2.d0 * dV)
          enddo

!         Last point.
          dIdV(e,s,NIVP)%el = (allcurr(e,s,NIVP)%el -                   &
               allcurr(e,s,NIVP-1)%el) / dV
          dIdV(e,s,NIVP)%symm = (allcurr(e,s,NIVP)%isymm -              &
               allcurr(e,s,NIVP-1)%isymm) / dV
          dIdV(e,s,NIVP)%asymm = (allcurr(e,s,NIVP)%iasymm -            &
               allcurr(e,s,NIVP-1)%iasymm) / dV

       enddo
    enddo


  end subroutine computedIdV


!  *******************************************************************  !
!                             computed2IdV2                             !
!  *******************************************************************  !
!  Description: compute the derivative of the differential conductance  !
!  by applying finite differences.                                      !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer nspin               : Number of spin components              !
!  integer NIVP                : Number of bias potential points        !
!  real*8 dV                   : Bias potential step                    !
!  integer NTenerg_div         : Number of energy grid points per node  !
!  *******************************************************************  !
  subroutine computed2IdV2

!
!   Modules
!
    use idsrdr_options,  only: nspin, NIVP, dV
    use idsrdr_engrid,   only: NTenerg_div

!   Local variables.
    integer :: e, s, v

    do e = 1,NTenerg_div
       do s = 1,nspin

!         First point.
          d2IdV2(e,s,1)%el = (dIdV(e,s,2)%el - dIdV(e,s,1)%el) / dV
          d2IdV2(e,s,1)%symm = (dIdV(e,s,2)%symm -                      &
               dIdV(e,s,1)%symm) / dV
          d2IdV2(e,s,1)%asymm = (dIdV(e,s,2)%asymm -                    &
               dIdV(e,s,1)%asymm) / dV

!         Intermediate points.
          do v = 2,NIVP-1
             d2IdV2(e,s,v-1)%el = (dIdV(e,s,v+1)%el -                   &
                  dIdV(e,s,v-1)%el) / (2.d0 * dV)
             d2IdV2(e,s,v-1)%symm = (dIdV(e,s,v+1)%symm -               &
                  dIdV(e,s,v-1)%symm) / (2.d0 * dV)
             d2IdV2(e,s,v-1)%asymm = (dIdV(e,s,v+1)%asymm -             &
                  dIdV(e,s,v-1)%asymm) / (2.d0 * dV)
          enddo

!         Last point.
          d2IdV2(e,s,NIVP)%el = (dIdV(e,s,NIVP)%el -                    &
               dIdV(e,s,NIVP-1)%el) / dV
          d2IdV2(e,s,NIVP)%symm = (dIdV(e,s,NIVP)%symm -                &
               dIdV(e,s,NIVP-1)%symm) / dV
          d2IdV2(e,s,NIVP)%asymm = (dIdV(e,s,NIVP)%asymm -              &
               dIdV(e,s,NIVP-1)%asymm) / dV

       enddo
    enddo


  end subroutine computed2IdV2


!  *******************************************************************  !
!                               freedIdV                                !
!  *******************************************************************  !
!  Description: free allocated memory.                                  !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer NIVP                : Number of bias potential points        !
!  *******************************************************************  !
  subroutine freedIdV

!
!   Modules
!
    use idsrdr_options,  only: NIVP

    if (NIVP < 3) return

!   Free memory.
    deallocate (dIdV)
    deallocate (d2IdV2)


  end subroutine freedIdV


!  *******************************************************************  !
!                              sumAlldIdV                               !
!  *******************************************************************  !
!  Description: function for suming two items of type 'alldIdV' (to be  !
!  used at MPI collective operations like 'MPI_Reduce').                !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  ****************************** INPUT ******************************  !
!  TYPE(alldIdV) indIdV(:,:,:)   : 1st item to be added                 !
!  integer len                   : Total length of the items            !
!                                  (1st dim x 2nd dim x 3rd dim)        !
!  integer type                  : Created MPI data type                !
!  ************************** INPUT/OUTPUT ***************************  !
!  TYPE(alldIdV) outdIdV(:,:,:) : 2nd item to be added                 !
!  *******************************************************************  !
  subroutine sumAlldIdV (indIdV, outdIdV, len, type)

!   Input variables.
    integer, intent(in) :: len, type
    TYPE(alldIdV), intent(in) :: indIdV(len)
    TYPE(alldIdV), intent(inout) :: outdIdV(len)

!   Local variables.
    integer :: i

    i = type ! I know that it seems useless...
    do i = 1,len
       outdIdV(i)%el = indIdV(i)%el + outdIdV(i)%el
       outdIdV(i)%symm = indIdV(i)%symm + outdIdV(i)%symm
       outdIdV(i)%asymm = indIdV(i)%asymm + outdIdV(i)%asymm
    enddo


  end subroutine sumAlldIdV


!  *******************************************************************  !


END MODULE idsrdr_conduct
