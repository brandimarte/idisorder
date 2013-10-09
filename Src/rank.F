! 
! Copyright (c) Smeagol Authors:
! A. R. Rocha, V. Garcia-Suarez, S. Bailey, C. J. Lambert, J. Ferrer and
! S. Sanvito 2003-2005
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! SMEAGOL IS DISTRIBUTED ONLY THROUGH THE OFICIAL WEBSITE (www.smeagol.tcd.ie)
! UPON COMPLETION OF THE "SMEAGOL ACADEMIC LICENSE" .
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
	SUBROUTINE RANK (n, A, B, side_rank)

C *****************************************************************
C Determines the direction of decimation by performing an
C SVD transformation on matrices A and B
C
C Written by Alexandre Reily Rocha, June 2003
C Computational Spintronics Group
C Trinity College Dublin
C e-mail: rochaa@tcd.ie
C ***************************** HISTORY ***********************************
C Original version:	June 2005
C ******************************************************************
	implicit none
	
	integer :: n
	double complex, dimension (n,n) :: A,B
	double complex, dimension (:,:) , allocatable  :: Mat
	double precision, dimension (:), allocatable :: S
	double complex :: vt, u
	double complex, dimension (:), allocatable :: WORK
	double precision, dimension (:), allocatable :: RWORK
	character(LEN=1) :: side_rank
	
	integer :: I,J,INFO,irankl,irankr

        allocate(Mat(2*n,n),WORK(2*n*64+2*n),RWORK(5*n),S(2*n))
        Mat=0.d0
        WORK=0.d0
        RWORK=0.d0
	
	Mat(1:n,:)=A
	mat(n+1:2*n,:)=B

	Call ZGESVD('N','N',2*n,n,Mat,2*N,S,U,1,VT,1,WORK,
     &   2*n*64+2*n,RWORK,INFO)
        irankl=0
        Do I=1,n
#ifdef SMEAGOLDEBUG	 
	 write(118,*) S(I)
#endif
	 If (S(I)/S(1) .GT. 1.0D-06) Then
	  irankl=irankl+1
	 EndIf
	EndDo
#ifdef SMEAGOLDEBUG
	write(118,*)
#endif
	
	Do I=1,n
	 Do J=1,n
	  Mat(I,J)=DCONJG(A(J,I))
	  Mat(n+I,J)=DCONJG(B(J,I))
	 EndDo
	EndDo
	
	Call ZGESVD('N','N',2*n,n,Mat,2*N,S,U,1,VT,1,WORK,
     &   2*n*64+2*n,RWORK,INFO)
     
        irankr=0
        Do I=1,n
#ifdef SMEAGOLDEBUG
	 write(118,*) S(I)
#endif
	 If (S(I)/S(1) .GT. 1.0D-6) Then
	  irankr=irankr+1
	 EndIf
	EndDo
#ifdef SMEAGOLDEBUG
        write(118,*)
#endif
	
	If (MIN(irankl,irankr) .EQ. n) Then
	 side_rank='0'
	ElseIf (irankl .LE. irankr) Then
	 side_rank='N'
	Else
	 side_rank='C'
	EndIf

        deallocate(Mat,WORK,RWORK,S)
	
	End Subroutine rank
