!  *******************************************************************  !
!  I-Disorder Fortran Code                                              !
!                                                                       !
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
!                          MODULE idsrdr_check                          !
!  *******************************************************************  !
!  Description: call lapack subroutines and check if execution was      !
!  successful (the idea is to avoid code repetition in the code).       !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *******************************************************************  !

MODULE idsrdr_check

!
!   Modules
!
#ifdef MPI
  use parallel,        only: 
#endif

  implicit none
  
  PUBLIC  :: CHECKzsytrf, CHECKzsytri, CHECKzgetrf, CHECKzgetri,        &
             CHECKzhetrf, CHECKzhetri


CONTAINS


!  *******************************************************************  !
!                              CHECKzsytrf                              !
!  *******************************************************************  !
!  Description: (Lapack computational rotine) computes Bunch-Kaufman    !
!  factorization of a complex symmetric matrix.                         !
!                                                                       !
!  Written by Pedro Brandimarte, Apr 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    April 2013                                      !
!  ****************************** INPUT ******************************  !
!  integer n           : Order of the squared matrix A                  !
!  character(1) LorU   : Indicates whether the upper or lower           !
!                        triangular part of the matrix 'A' is stored    !
!  ************************** INPUT/OUTPUT ***************************  !
!  complex*8 A(n,n)    : Square symmetric matrix to be factorized as    !
!                        A = P*L*D*L^T*P^T and overwritten by details   !
!                        of the block-diagonal matrix D and the         !
!                        multipliers used to obtain the factor L        !
!  ***************************** OUTPUT ******************************  !
!  integer ipiv(n)     : Contains details of the interchanges and the   !
!                        block structure of D                           !
!  *******************************************************************  !
  subroutine CHECKzsytrf (n, LorU, A, ipiv)

!
!   Modules
!
#ifdef MPI
    use parallel, only: MPI_Comm_MyWorld

    include "mpif.h"
#endif

!   Input variables.
    integer, intent(in) :: n
    integer, intent(out) :: ipiv(n)
    complex(8), intent(inout) :: A(n,n)
    character(len=1), intent(in) :: LorU

!   Local variables.
    integer :: info, lwork
    complex(8) :: l
    complex(8), dimension(:), allocatable :: work
    external :: zsytrf
#ifdef MPI
    integer MPIerror
#endif

!   Workspace query: calculates the optimal size of 'work'.
    lwork = -1
    call zsytrf (LorU, n, A, n, ipiv, l, lwork, info)
    lwork = INT(l)
    allocate (work(lwork))

!   Computes the Bunch-Kaufman factorization of 'A'.
    call zsytrf (LorU, n, A, n, ipiv, work, lwork, info)

    if (info .lt. 0) then
       write(0,'(a,/,a,i,a)') "ERROR: In lapack zsytrf: ",              &
            " The ", -info, " argument had an illegal value!"
#ifdef MPI
       call MPI_Abort (MPI_Comm_MyWorld, 1, MPIerror)
#endif
       stop
    elseif (info .gt. 0) then
       write(0,'(a,/,a,/,a,/,a)') "ERROR: In lapack zsytrf: ",          &
            " The factorization has been completed, but D is exactly",  &
            " singular, so division by 0 will occur if the D is",       &
            " used for solving a system of linear equations!"
#ifdef MPI
       call MPI_Abort (MPI_Comm_MyWorld, 1, MPIerror)
#endif
       stop
    endif

!   Free memory.
    deallocate (work)


  end subroutine CHECKzsytrf


!  *******************************************************************  !
!                              CHECKzsytri                              !
!  *******************************************************************  !
!  Description: (Lapack computational rotine) computes the inverse of   !
!  an Bunch-Kaufman-factored complex symmetric matrix A.                !
!                                                                       !
!  Written by Pedro Brandimarte, Apr 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    April 2013                                      !
!  ****************************** INPUT ******************************  !
!  integer n           : Order of the squared matrix A                  !
!  character(1) LorU   : Indicates whether the upper or lower           !
!                        triangular part of the matrix 'A' is stored    !
!  integer ipiv(n)     : Contains details of the interchanges and the   !
!                        block structure of D as returned by 'zsytrf'   !
!  ************************** INPUT/OUTPUT ***************************  !
!  complex*8 A(n,n)    : Symmetric matrix factored as                   !
!                        A = P*L*D*L^T*P^T (returned by 'zsytrf') and   !
!                        overwritten by the inverse of 'A'              !
!  *******************************************************************  !
  subroutine CHECKzsytri (n, LorU, A, ipiv)

!
!   Modules
!
#ifdef MPI
    use parallel, only: MPI_Comm_MyWorld

    include "mpif.h"
#endif

!   Input variables.
    integer, intent(in) :: n
    integer, intent(in) :: ipiv(n)
    complex(8), intent(inout) :: A(n,n)
    character(len=1), intent(in) :: LorU

!   Local variables.
    integer :: info, i, j
    complex(8), dimension(:), allocatable :: work
    external :: zsytri
#ifdef MPI
    integer MPIerror
#endif

!   Allocate a workspace array.
    allocate (work(n))

!   Computes the inverse of 'A'.
    call zsytri (LorU, n, A, n, ipiv, work, info)

    if (info .lt. 0) then
       write(0,'(a,/,a,i,a)') "ERROR: In lapack zhetri: ",              &
            " The ", -info, " argument had an illegal value!"
#ifdef MPI
       call MPI_Abort (MPI_Comm_MyWorld, 1, MPIerror)
#endif
       stop
    elseif (info .gt. 0) then
       write(0,'(a,/,a,i,a,/,a)') "ERROR: In lapack zhetri: ",          &
            " The ", info, "-th diagonal element of D is zero, D is",   &
            " singular, and the inversion could not be completed!"
#ifdef MPI
       call MPI_Abort (MPI_Comm_MyWorld, 1, MPIerror)
#endif
       stop
    endif

!   Copy the computed triangular part to the symmetric one.
    if (LorU == 'L') then
       do i = 1,n-1
          do j = i+1,n
             A(i,j) = A(j,i)
          enddo
       enddo
    else
       do j = 1,n-1
          do i = j+1,n
             A(i,j) = A(j,i)             
          enddo
       enddo
    endif

!   Free memory.
    deallocate (work)


  end subroutine CHECKzsytri


!  *******************************************************************  !
!                              CHECKzgetrf                              !
!  *******************************************************************  !
!  Description: (Lapack computational rotine) computes the LU           !
!  factorization of a complex general squared matrix A.                 !
!                                                                       !
!  Written by Pedro Brandimarte, Feb 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    February 2013                                   !
!  ****************************** INPUT ******************************  !
!  integer n           : Order of the squared matrix A                  !
!  ************************** INPUT/OUTPUT ***************************  !
!  complex*8 A(n,n)    : Square general matrix to be factorized as      !
!                        A = P*L*U and overwritten by the factors L     !
!                        and U (the unit diagonal elements of L are     !
!                        not stored)                                    !
!  ***************************** OUTPUT ******************************  !
!  integer ipiv(n)     : The pivot indices (for 1 <= i <= n, row 'i'    !
!                        is interchanged with row 'ipiv(i)')            !
!  *******************************************************************  !
  subroutine CHECKzgetrf (n, A, ipiv)

!
!   Modules
!
#ifdef MPI
    use parallel, only: MPI_Comm_MyWorld

    include "mpif.h"
#endif

!   Input variables.
    integer, intent(in) :: n
    integer, intent(out) :: ipiv(n)
    complex(8), intent(inout) :: A(n,n)

!   Local variables.
    integer :: info
    external :: zgetrf
#ifdef MPI
    integer MPIerror
#endif

!   Computes the LU factorization of 'A'.
    call zgetrf (n, n, A, n, ipiv, info)

    if (info .lt. 0) then
       write(0,'(a,/,a,i,a)') "ERROR: In lapack zgetrf: ",              &
            " The ", -info, " argument had an illegal value!"
#ifdef MPI
       call MPI_Abort (MPI_Comm_MyWorld, 1, MPIerror)
#endif
       stop
    elseif (info .gt. 0) then
       write(0,'(a,/,a,/,a,/,a)') "ERROR: In lapack zgetrf: ",          &
            " The factorization has been completed, but U is exactly",  &
            " singular, so division by 0 will occur if the factor U",   &
            " is used for solving a system of linear equations!"
#ifdef MPI
       call MPI_Abort (MPI_Comm_MyWorld, 1, MPIerror)
#endif
       stop
    endif


  end subroutine CHECKzgetrf


!  *******************************************************************  !
!                              CHECKzgetri                              !
!  *******************************************************************  !
!  Description: (Lapack computational rotine) computes the inverse of   !
!  an LU-factored complex general squared matrix A.                     !
!                                                                       !
!  Written by Pedro Brandimarte, Feb 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    February 2013                                   !
!  ****************************** INPUT ******************************  !
!  integer n           : Order of the squared matrix A                  !
!  integer ipiv(n)     : The pivot indices (for 1 <= i <= n, row 'i'    !
!                        is interchanged with row 'ipiv(i)')            !
!  ************************** INPUT/OUTPUT ***************************  !
!  complex*8 A(n,n)    : Square general matrix factored as A = P*L*U    !
!                        (returned by 'zgetrf') and overwritten by the  !
!                        inverse of 'A'                                 !
!  *******************************************************************  !
  subroutine CHECKzgetri (n, A, ipiv)

!
!   Modules
!
#ifdef MPI
    use parallel, only: MPI_Comm_MyWorld

    include "mpif.h"
#endif

!   Input variables.
    integer, intent(in) :: n
    integer, intent(in) :: ipiv(n)
    complex(8), intent(inout) :: A(n,n)

!   Local variables.
    integer :: info, lwork
    complex(8) :: l
    complex(8), dimension(:), allocatable :: work
    external :: zgetri
#ifdef MPI
    integer MPIerror
#endif

!   Workspace query: calculates the optimal size of 'work'.
    lwork = -1
    call zgetri (n, A, n, ipiv, l, lwork, info)
    lwork = INT(l)
    allocate (work(lwork))

!   Computes the inverse of 'A'.
    call zgetri (n, A, n, ipiv, work, lwork, info)

    if (info .lt. 0) then
       write(0,'(a,/,a,i,a)') "ERROR: In lapack zgetri: ",              &
            " The ", -info, " argument had an illegal value!"
#ifdef MPI
       call MPI_Abort (MPI_Comm_MyWorld, 1, MPIerror)
#endif
       stop
    elseif (info .gt. 0) then
       write(0,'(a,/,a,i,a,/,a,/,a)') "ERROR: In lapack zgetri: ",      &
            " The ", info, "-th diagonal element of the factor U",      &
            " is zero, U is singular, and the inversion could not",     &
            " be completed!"
#ifdef MPI
       call MPI_Abort (MPI_Comm_MyWorld, 1, MPIerror)
#endif
       stop
    endif

!   Free memory.
    deallocate (work)


  end subroutine CHECKzgetri


!  *******************************************************************  !
!                              CHECKzhetrf                              !
!  *******************************************************************  !
!  Description: (Lapack computational rotine) computes Bunch-Kaufman    !
!  factorization of a complex Hermitian matrix.                         !
!                                                                       !
!  Written by Pedro Brandimarte, Apr 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    April 2013                                      !
!  ****************************** INPUT ******************************  !
!  integer n           : Order of the squared matrix A                  !
!  character(1) LorU   : Indicates whether the upper or lower           !
!                        triangular part of the matrix 'A' is stored    !
!  ************************** INPUT/OUTPUT ***************************  !
!  complex*8 A(n,n)    : Square Hermitian matrix to be factorized as    !
!                        A = P*L*D*L^H*P^T and overwritten by details   !
!                        of the block-diagonal matrix D and the         !
!                        multipliers used to obtain the factor L        !
!  ***************************** OUTPUT ******************************  !
!  integer ipiv(n)     : Contains details of the interchanges and the   !
!                        block structure of D                           !
!  *******************************************************************  !
  subroutine CHECKzhetrf (n, LorU, A, ipiv)

!
!   Modules
!
#ifdef MPI
    use parallel, only: MPI_Comm_MyWorld

    include "mpif.h"
#endif

!   Input variables.
    integer, intent(in) :: n
    integer, intent(out) :: ipiv(n)
    complex(8), intent(inout) :: A(n,n)
    character(len=1), intent(in) :: LorU

!   Local variables.
    integer :: info, lwork
    complex(8) :: l
    complex(8), dimension(:), allocatable :: work
    external :: zhetrf
#ifdef MPI
    integer MPIerror
#endif

!   Workspace query: calculates the optimal size of 'work'.
    lwork = -1
    call zhetrf (LorU, n, A, n, ipiv, l, lwork, info)
    lwork = INT(l)
    allocate (work(lwork))

!   Computes the Bunch-Kaufman factorization of 'A'.
    call zhetrf (LorU, n, A, n, ipiv, work, lwork, info)

    if (info .lt. 0) then
       write(0,'(a,/,a,i,a)') "ERROR: In lapack zhetrf: ",              &
            " The ", -info, " argument had an illegal value!"
#ifdef MPI
       call MPI_Abort (MPI_Comm_MyWorld, 1, MPIerror)
#endif
       stop
    elseif (info .gt. 0) then
       write(0,'(a,/,a,/,a,/,a)') "ERROR: In lapack zhetrf: ",          &
            " The factorization has been completed, but D is exactly",  &
            " singular, so division by 0 will occur if the D is",       &
            " used for solving a system of linear equations!"
#ifdef MPI
       call MPI_Abort (MPI_Comm_MyWorld, 1, MPIerror)
#endif
       stop
    endif

!   Free memory.
    deallocate (work)


  end subroutine CHECKzhetrf


!  *******************************************************************  !
!                              CHECKzhetri                              !
!  *******************************************************************  !
!  Description: (Lapack computational rotine) computes the inverse of   !
!  an Bunch-Kaufman-factored complex Hermitian matrix A.                !
!                                                                       !
!  Written by Pedro Brandimarte, Apr 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    April 2013                                      !
!  ****************************** INPUT ******************************  !
!  integer n           : Order of the squared matrix A                  !
!  character(1) LorU   : Indicates whether the upper or lower           !
!                        triangular part of the matrix 'A' is stored    !
!  integer ipiv(n)     : Contains details of the interchanges and the   !
!                        block structure of D as returned by 'zhetrf'   !
!  ************************** INPUT/OUTPUT ***************************  !
!  complex*8 A(n,n)    : Hermitian matrix factored as                   !
!                        A = P*L*D*L^H*P^T (returned by 'zhetrf') and   !
!                        overwritten by the inverse of 'A'              !
!  *******************************************************************  !
  subroutine CHECKzhetri (n, LorU, A, ipiv)

!
!   Modules
!
#ifdef MPI
    use parallel, only: MPI_Comm_MyWorld

    include "mpif.h"
#endif

!   Input variables.
    integer, intent(in) :: n
    integer, intent(in) :: ipiv(n)
    complex(8), intent(inout) :: A(n,n)
    character(len=1), intent(in) :: LorU

!   Local variables.
    integer :: info
    complex(8), dimension(:), allocatable :: work
    external :: zhetri
#ifdef MPI
    integer MPIerror
#endif

!   Allocate a workspace array.
    allocate (work(n))

!   Computes the inverse of 'A'.
    call zhetri (LorU, n, A, n, ipiv, work, info)

    if (info .lt. 0) then
       write(0,'(a,/,a,i,a)') "ERROR: In lapack zhetri: ",              &
            " The ", -info, " argument had an illegal value!"
#ifdef MPI
       call MPI_Abort (MPI_Comm_MyWorld, 1, MPIerror)
#endif
       stop
    elseif (info .gt. 0) then
       write(0,'(a,/,a,i,a,/,a)') "ERROR: In lapack zhetri: ",          &
            " The ", info, "-th diagonal element of D is zero, D is",   &
            " singular, and the inversion could not be completed!"
#ifdef MPI
       call MPI_Abort (MPI_Comm_MyWorld, 1, MPIerror)
#endif
       stop
    endif

!   Free memory.
    deallocate (work)


  end subroutine CHECKzhetri


!  *******************************************************************  !


END MODULE idsrdr_check
