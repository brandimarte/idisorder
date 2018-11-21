!  *******************************************************************  !
!  I-Disorder Fortran Code 2007-2014                                    !
!                                                                       !
!  Written by Pedro Brandimarte (brandimarte@gmail.com),                !
!             Alberto Torres (alberto.trj@gmail.com) and                !
!             Alexandre Reily Rocha (reilya@ift.unesp.br).              !
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
!  Description: computes a Hilbert transform H{Sr(w')}(w) on an         !
!  equidistant grid by discrete convolution with FFT.                   !
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
  use MKL_DFTI

  implicit none
  
  PUBLIC  :: hilbertkernel, hilbert, freehilb
  PRIVATE :: ker

  complex(8), allocatable, dimension (:) :: ker ! kernel function


CONTAINS


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
!  complex*8 ker(2*nAsymmPts)  : Kernel function                        !
!  *******************************************************************  !
  subroutine hilbertkernel

!
!   Modules
!
#ifdef MPI
    use parallel,        only: IOnode, MPI_Comm_MyWorld
#else
    use parallel,        only: IOnode
#endif
    use idsrdr_options,  only: nAsymmPts

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
          ker(i+1) = DCMPLX(foo(i+2) - 2.d0*foo(i+1) + foo(i))
          ker(2*nAsymmPts-i+1) = - ker(i+1)
       enddo

!      Allocate and initialize the descriptor data structure.
       MKLstatus = DftiCreateDescriptor (MKLdesc, DFTI_DOUBLE,          &
                                         DFTI_COMPLEX, 1, 2*nAsymmPts)

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
    call MPI_Bcast (ker, 2*nAsymmPts, MPI_Double_Complex, 0,            &
                    MPI_Comm_MyWorld, MPIerror)
#endif


  end subroutine hilbertkernel


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
!  integer nAsymmPts           : Number of energy grid points           !
!                                for asymmetric term integral           !
!  complex*8 aux(2*nAsymmPts)  : Array where the Hibert transform will  !
!                                be performed                           !
!  *******************************************************************  !
  subroutine hilbert (nAsymmPts, aux)

!   Input variables.
    integer, intent(in) :: nAsymmPts
    complex(8), dimension (2*nAsymmPts), intent(inout) :: aux

!   Local variables.
    integer :: i
    real(8) :: foo

!   Intel MKL types.
    TYPE(DFTI_DESCRIPTOR), pointer :: MKLdesc
    integer :: MKLstatus

!   Allocate and initialize the descriptor data structure.
    MKLstatus = DftiCreateDescriptor (MKLdesc, DFTI_DOUBLE,             &
                                      DFTI_COMPLEX, 1, 2*nAsymmPts)

!   Set the scale factor for the backward transform
!   (to make it really the inverse of the forward).
    foo = 1.d0 / (2.d0 * nAsymmPts)
    MKLstatus = DftiSetValue (MKLdesc, DFTI_BACKWARD_SCALE, foo)

!   Complete initialization of the previously created descriptor.
    MKLstatus = DftiCommitDescriptor (MKLdesc)

!   Compute the forward FFT.
    MKLstatus = DftiComputeForward (MKLdesc, aux)

    do i = 1,2*nAsymmPts
       aux(i) = ker(i) * aux(i)
    enddo

!   Compute the backward FFT.
    MKLstatus = DftiComputeBackward (MKLdesc, aux)

!   Free memory.
    MKLstatus = DftiFreeDescriptor (MKLdesc)


  end subroutine hilbert


!  *******************************************************************  !
!                               freehilb                                !
!  *******************************************************************  !
!  Description: free allocated vector.                                  !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2013                                   !
!  *******************************************************************  !
  subroutine freehilb


!   Free memory.
    deallocate (ker)


  end subroutine freehilb


!  *******************************************************************  !


END MODULE idsrdr_hilbert
