!  *******************************************************************  !
!  Copyright (c) Smeagol Authors (2003-2005):  A. R. Rocha,             !
!                                              V. Garcia-Suarez,        !
!                                              S. Bailey,               !
!                                              C. J. Lambert,           !
!                                              J. Ferrer and            !
!                                              S. Sanvito               !
!                                                                       !
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS  !
!  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT    !
!  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS    !
!  FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE      !
!  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,  !
!  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES             !
!  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR   !
!  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)   !
!  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,  !
!  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)        !
!  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED  !
!  OF THE POSSIBILITY OF SUCH DAMAGE.                                   !
!                                                                       !
!  SMEAGOL IS DISTRIBUTED ONLY THROUGH THE OFICIAL WEBSITE              !
!  (www.smeagol.tcd.ie) UPON COMPLETION OF "SMEAGOL ACADEMIC LICENSE".  !
!                                                                       !
!  FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL smeagol@tcd.ie  !
!  *******************************************************************  !

!  *******************************************************************  !
!                              MATRIXMULT                               !
!  *******************************************************************  !
!  Description:                                                         !
!                                                                       !
!  Written by                                                           !
!  ***************************** HISTORY *****************************  !
!  Original version:                                                    !
!  ****************************** INPUT ******************************  !
!  complex*8 A(N1,N2)            :                                      !
!  complex*8 B(N2,N3)            :                                      !
!  integer N1                    :                                      !
!  integer N2                    :                                      !
!  integer N3                    :                                      !
!  ************************** INPUT/OUTPUT ***************************  !
!  complex*8 C(N1,N3)            :                                      !
!  *******************************************************************  !
      subroutine MATRIXMULT (A, B, C, N1, N2, N3)

      implicit none

!     Input variables.
      INTEGER :: N1, N2, N3
      COMPLEX(8), DIMENSION (N1,N2) :: A
      COMPLEX(8), DIMENSION (N2,N3) :: B
      COMPLEX(8), DIMENSION (N1,N3) :: C

!     Local variables.
      INTEGER :: I, J, L

      C = (0.D0,0.D0)
      DO I = 1,N1
         DO J = 1,N3
            DO L = 1,N2
               C(I,J) = C(I,J)+A(I,L)*B(L,J)
            ENDDO
         ENDDO
      ENDDO


      end subroutine


!  *******************************************************************  !
!                              MATRIXMUL_1                              !
!  *******************************************************************  !
!  Description:                                                         !
!                                                                       !
!  Written by                                                           !
!  ***************************** HISTORY *****************************  !
!  Original version:                                                    !
!  ****************************** INPUT ******************************  !
!  complex*8 B(N1,N1)            :                                      !
!  integer N1                    :                                      !
!  ************************** INPUT/OUTPUT ***************************  !
!  complex*8 A(N1,N1)            :                                      !
!  *******************************************************************  !
      subroutine MATRIXMUL_1 (A, B, N1)

      implicit none

!     Input variables.
      INTEGER :: N1
      COMPLEX(8), DIMENSION (N1,N1) :: A, B

!     Local variables.
      INTEGER :: I, J, L
      COMPLEX(8), DIMENSION (N1,N1) :: C

      C = (0.D0,0.D0)
      DO I = 1,N1
         DO J = 1,N1
            DO L = 1,N1
               C(I,J) = C(I,J)+A(I,L)*B(L,J)
            ENDDO
         ENDDO
      ENDDO
      A = C


      end subroutine
