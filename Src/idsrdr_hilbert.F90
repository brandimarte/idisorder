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
!                         MODULE idsrdr_hilbert                         !
!  *******************************************************************  !
!  Description: calculate the required non-equilibrium green's          !
!  functions with recursive technique.                                  !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2013                                   !
!  *******************************************************************  !

MODULE idsrdr_hilbert

!
!   Modules
!
  use parallel,        only: 
  use idsrdr_options,  only: 
  use idsrdr_ephcoupl, only: 
  use idsrdr_recipes,  only: 

  implicit none
  
  PUBLIC  :: hilbertinit, hilbert, freehilb, hilbEn, hilbWe
  PRIVATE :: hilbertkernel, hilbertgrid, hilbertmodes, gridPts, kbTol

  integer :: gridPts ! number of energy grid points per node

  real(8) :: hghFreq ! highest vibrational mode energy
  real(8) :: kbTol = 18.d0 ! tolerance value to take into account

  real(8), allocatable, dimension (:) :: hilbEn ! energy grid points
  real(8), allocatable, dimension (:) :: hilbWe ! energy grid weights

  real(8), allocatable, dimension (:) :: ker ! Hilbert transform kernel

  TYPE hilbVal
     real(8), pointer :: HIL(:) ! pointer to Hilbert transform array
  END TYPE hilbVal

  TYPE(hilbVal), dimension(:), allocatable :: hilb ! Hilbert transforms


CONTAINS


!  *******************************************************************  !
!                              hilbertinit                              !
!  *******************************************************************  !
!  Description: allocates grid energy points and weights arrays, and    !
!  Hilbert transform structure. Also allocates and calculates the       !
!  kernel interpolation function, and computes the highest vibrational  !
!  mode energy.                                                         !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer Nodes               : Total number of nodes (MPI_Comm_size)  !
!  integer nAsymmPts           : Number of energy grid points           !
!                                for asymmetric term integral           !
!  integer neph                : Number of units with e-ph interaction  !
!  TYPE(ephFreq) freq(neph)%F(nModes) : [real*8] Vibrational mode's     !
!                                       frequencies                     !
!  *******************************************************************  !
  subroutine hilbertinit

!
!   Modules
!
    use parallel,        only: Nodes
    use idsrdr_options,  only: nAsymmPts
    use idsrdr_ephcoupl, only: neph, freq

!   Local variables.
    integer :: i

!   Set the number of grid points per node.
#ifdef MPI
    gridPts = nAsymmPts / Nodes
#else
    gridPts = nAsymmPts
#endif

!   Allocate the energy grid points and weights arrays.
    allocate (hilbEn(gridPts), hilbWe(gridPts))

!   Compute interpolation kernel function.
    call hilbertkernel

!   Allocate Hilbert transform array.
    allocate (hilb(neph))

!   Find the highest vibrational mode energy.
    hghFreq = freq(1)%F(1)
    do i = 2,neph
       if (hghFreq < freq(i)%F(1)) hghFreq = freq(i)%F(1)
    enddo

  end subroutine hilbertinit


!  *******************************************************************  !
!                                hilbert                                !
!  *******************************************************************  !
!  Description: main function that comput the Hilbert transform.        !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2013                                   !
!  ****************************** INPUT ******************************  !
!  real*8 Efermi               : Energy around where the integration    !
!                                will be performed                      !
!  *******************************************************************  !
  subroutine hilbert (Efermi)

!   Input variables.
    real(8), intent(in) :: Efermi

!   Calculate the energy grid.
    call hilbertgrid (Efermi)

!   Calculate the Hilbert transform for each vibrational mode.
    call hilbertmodes


  end subroutine hilbert


!  *******************************************************************  !
!                              hilbertgrid                              !
!  *******************************************************************  !
!  Description: create the energy grid for Hilbert transform at the     !
!  integral from asymmetric term of the current expression.             !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer Nodes               : Total number of nodes (MPI_Comm_size)  !
!  integer nAsymmPts           : Number of energy grid points           !
!                                for asymmetric term integral           !
!  real*8 temp                 : Electronic temperature                 !
!  real*8 VFinal               : Final value of the bias potential      !
!  integer neph                : Number of units with e-ph interaction  !
!  TYPE(ephFreq) freq(neph)%F(nModes) : [real*8] Vibrational mode's     !
!                                       frequencies                     !
!  ****************************** INPUT ******************************  !
!  real*8 Efermi               : Energy around where the integration    !
!                                will be performed                      !
!  *******************************************************************  !
  subroutine hilbertgrid (Efermi)

!
!   Modules
!
    use parallel,        only: IOnode, Nodes
    use idsrdr_options,  only: nAsymmPts, temp, VFinal
    use idsrdr_ephcoupl, only: neph, freq
    use idsrdr_recipes,  only: RECPSsimpson

!!$#ifdef MPI
!!$    include "mpif.h"
!!$#endif

!   Input variables.
    real(8), intent(in) :: Efermi

!   Local variables.
    real(8) :: enI, enF
!!$    real(8), allocatable, dimension (:) :: EX, EW
!!$#ifdef MPI
!!$    integer :: MPIerror ! Return error code in MPI routines
!!$    integer, dimension(MPI_Status_Size) :: MPIstatus
!!$#endif

    if (IOnode) write (6,'(/,a)')                                       &
         '         hilbertgrid: Computing Hilbert transform...'

!      Set lower and upper limit of energy integration.
    enI = Efermi - hghFreq - kbTol*temp
    enF = Efermi + hghFreq + kbTol*temp + VFinal
    if (IOnode) write (6,'(/,a,f12.4)')                                 &
         '         hilbertgrid: initial energy = ', enI
    if (IOnode) write (6,'(a,f12.4)')                                   &
         '         hilbertgrid: final energy = ', enF

!!$!   Allocate full size energy points and weights arrays.
!!$    allocate (EX(nAsymmPts), EW(nAsymmPts))

!   Compute energy points and weights on an equidistant grid.
    call RECPSsimpson (enI, enF, nAsymmPts, hilbEn, hilbWe)
!!$    call RECPSsimpson (enI, enF, nAsymmPts, EX, EW)

!!$!   Energy grid for node 0.
!!$    hilbEn = EX(1:gridPts)
!!$    hilbWe = EW(1:gridPts)

!!$#ifdef MPI
!!$!      Distribute 'EX' and 'EW' to the other nodes.
!!$       do i = 1,Nodes-1
!!$          call MPI_Send (EX(i*gridPts+1:(i+1)*gridPts), gridPts,        &
!!$                         MPI_Double_Precision, i, 1,                    &
!!$                         MPI_Comm_world, MPIerror)
!!$          call MPI_Send (EW(i*gridPts+1:(i+1)*gridPts), gridPts,        &
!!$                         MPI_Double_Precision, i, 2,                    &
!!$                         MPI_Comm_world, MPIerror)
!!$       enddo
!!$#endif
!!$
!!$!      Free memory.
!!$       deallocate (EX, EW)

    if (IOnode) write(6,'(/,a,/)') '         hilbertgrid: done!'

!!$#ifdef MPI
!!$    else
!!$
!!$!      Receive 'EX' and 'EW' from node 0.
!!$       call MPI_Recv (hilbEn, gridPts, MPI_Double_Precision, 0, 1,      &
!!$                      MPI_Comm_world, MPIstatus, MPIerror)
!!$       call MPI_Recv (hilbWe, gridPts, MPI_Double_Precision, 0, 2,      &
!!$                      MPI_Comm_world, MPIstatus, MPIerror)
!!$#endif
!!$
!!$    endif ! if (IOnode)


  end subroutine hilbertgrid


!  *******************************************************************  !
!                             hilbertmodes                              !
!  *******************************************************************  !
!  Description: compute the Hilbert transform for each vibrational      !
!  mode.                                                                !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer neph                : Number of units with e-ph interaction  !
!  integer nModes(neph)        : Number of vibrational modes            !
!  integer nAsymmPts           : Number of energy grid points           !
!                                for asymmetric term integral           !
!  *******************************************************************  !
  subroutine hilbertmodes

!
!   Modules
!
    use idsrdr_ephcoupl, only: neph, nModes
    use idsrdr_options,  only: nAsymmPts
    use MKL_DFTI

!   Local variables.
    integer :: i, j, k, e
    real(8) :: foo
    real(8), allocatable, dimension (:) :: aux

!   Intel MKL types.
    TYPE(DFTI_DESCRIPTOR), pointer :: MKLdesc
    integer :: MKLstatus

!   Allocate auxiliary array.
    allocate (aux(2*nAsymmPts))

!   Allocate and initialize the descriptor data structure.
    MKLstatus = DftiCreateDescriptor (MKLdesc, DFTI_DOUBLE,             &
                                      DFTI_REAL, 1, 2*nAsymmPts)

!   Set the scale factor for the backward transform
!   (to make it really the inverse of the forward).
    foo = 1.d0 / (2.d0 * nAsymmPts)
    MKLstatus = DftiSetValue (MKLdesc, DFTI_BACKWARD_SCALE, foo)

!   Complete initialization of the previously created descriptor.
    MKLstatus = DftiCommitDescriptor (MKLdesc)

    do i = 1,neph ! over unit with e-ph

!      Allocate Hilbert array for i-th unit.
       allocate (hilb(i)%HIL(nModes(i)))

       do j = 1,nModes(i) ! over phonon modes

!          init aux
          aux = 0.d0
          do k = 1,nAsymmPts
!!$             aux(k) = 
          enddo

!         Compute the forward FFT.
          MKLstatus = DftiComputeForward (MKLdesc, aux)

          do e = 1,2*nAsymmPts
             aux(e) = ker(e) * aux(e)
          enddo

!         Compute the forward FFT.
          MKLstatus = DftiComputeBackward (MKLdesc, aux)

       enddo

    enddo ! do i = 1,neph

!   Free memory.
    deallocate (aux)
    MKLstatus = DftiFreeDescriptor (MKLdesc)


  end subroutine hilbertmodes


!  *******************************************************************  !
!                             hilbertkernel                             !
!  *******************************************************************  !
!  Description: Computes the kernel function (and its Fourier           !
!  transform) associated with a linear interpolation on an equidistant  !
!  energy grid.                                                         !
!                                                                       !
!  Uses Intel MKL FFT algorithm.                                        !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer nAsymmPts           : Number of energy grid points           !
!                                for asymmetric term integral           !
!  ***************************** OUTPUT ******************************  !
!  real*8 ker(nAsymmPts)       : Kernel function                        !
!  *******************************************************************  !
  subroutine hilbertkernel

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_options,  only: nAsymmPts
    use MKL_DFTI

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    integer :: i
    real(8), parameter :: pi = 3.1415926535897932384626433832795028841D0
    real(8), dimension (:), allocatable :: foo
#ifdef MPI
    integer :: MPIerror
#endif

!   Intel MKL types.
    TYPE(DFTI_DESCRIPTOR), pointer :: MKLdesc
    integer :: MKLstatus

!   Allocate Hilbert transform kernel array.
    allocate (ker(2*nAsymmPts))

    if (IOnode) then

!      Allocate auxiliary vector.
       allocate (foo(nAsymmPts+1))

!      Build the kernel function.
       foo(1) = 0.d0
       do i = 1,nAsymmPts
          foo(i+1) = 1.d0 * i * DLOG (1.d0 * i)
       enddo

       ker(1) = 0.d0
       ker(nAsymmPts+1) = 0.d0
       do i = 1,nAsymmPts-1
          ker(i+1) = foo(i+2) - 2.d0*foo(i+1) + foo(i)
          ker(2*nAsymmPts-i+1) = - ker(i+1)
       enddo

!      Allocate and initialize the descriptor data structure.
       MKLstatus = DftiCreateDescriptor (MKLdesc, DFTI_DOUBLE,          &
                                         DFTI_REAL, 1, 2*nAsymmPts)

!      Complete initialization of the previously created descriptor.
       MKLstatus = DftiCommitDescriptor (MKLdesc)

!      Compute the forward FFT.
       MKLstatus = DftiComputeForward (MKLdesc, ker)

       ker = - ker / pi

!      Free memory
       MKLstatus = DftiFreeDescriptor (MKLdesc)
       deallocate(foo)

    endif

!   Distribute 'ker' to all nodes.
#ifdef MPI
    call MPI_Bcast (ker, 2*nAsymmPts, MPI_Double_Precision, 0,          &
                    MPI_Comm_world, MPIerror)
#endif


  end subroutine hilbertkernel


!  *******************************************************************  !
!                               freehilb                                !
!  *******************************************************************  !
!  Description: free allocated vectors.                                 !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer neph                : Number of units with e-ph interaction  !
!  *******************************************************************  !
  subroutine freehilb

!
!   Modules
!
    use idsrdr_ephcoupl, only: neph

!   Local variables.
    integer :: i

!   Free memory.
    deallocate (hilbEn, hilbWe)
    deallocate (ker)

!   First deallocates pointed arrays and matrices.
    do i = 1,neph
       deallocate (hilb(i)%HIL)
    enddo
    deallocate (hilb)


  end subroutine freehilb


!  *******************************************************************  !


END MODULE idsrdr_hilbert
