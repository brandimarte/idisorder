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
! UPON COMPLETION OF THE "SMEAGOL ACADEMIC LICENSE".
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
	SUBROUTINE GENSVD(n,S1Imp,H1Imp,side_rank,number,VH,VS,Q)
	
C *****************************************************************
C Calculates the Generalized Singular Value Decomposiation of 
C a pair of matrices.
C For references see: 
C
C Written by Alexandre Reily Rocha, June 2003 
C Computational Spintronics Group
C Trinity College Dublin
C e-mail: rochaa@tcd.ie
C ***************************** HISTORY ***********************************
C Original version:	June 2003
C ***************************** INPUT ***************************** 	
C integer n			: Dimension of the matrices
C double complex S1Imp(n,n) 	: Overlap matrix S1
C double complex H1Imp(n,n) 	: Coupling matrix H1
C character side_rank		: Direction of GSVD
C ***************************** OUTPUT ****************************
C integer number		: dimension of non-singular matrix
C double complex VH(n,n)	: GSVD transformed H1 (U*H1)
C double complex VS(n,n)	: GSVD transformed S1 (Vdag*S1)
C double complex Q(n,n)		: Singular value matrix
C *********************** AUXILIARY ******************************
C double complex S1E(n,n)	: Auxiliary overlap matrix
C double complex H1E(n,n)	: Auxiliary coupling matrix
C ****************************************************************

	implicit none
	
	double precision, parameter :: TOLA=3.0D-9,TOLB=3.0D-9
	character(LEN=1) :: side_rank
	integer :: n, number
	double complex, dimension (:,:), allocatable :: S1E, H1E
	double complex, dimension (n,n) :: S1Imp, H1Imp
	double complex, dimension (n,n) :: VH,VS,Q
	double complex, dimension (:,:), allocatable :: U, VMat
	double complex, dimension (:), allocatable :: work
	double precision, dimension (:), allocatable :: rwork
	double complex, dimension (:), allocatable   :: tau
	integer, dimension (n) :: ipiv
	
	integer :: info, K, L, I, J
	
        allocate(S1E(n,n),H1E(n,n),U(n,n),VMat(n,n))
        allocate(work(3*n),rwork(2*n),tau(n))

	If (side_rank .EQ. 'N') Then
	 H1E=H1Imp
	 S1E=S1Imp
	Else
	 Do I=1,N
	  Do J=1,N
	   H1E(I,J)=DCONJG(H1Imp(J,I))
	   S1E(I,J)=DCONJG(S1Imp(J,I))
	  EndDo
	 EndDo
	EndIf
	
	Call ZGGSVP('U','V','Q',n,n,n,H1E,n,S1E,n,TOLA,TOLB,K,L,
     &    U,n,Vmat,n,Q,n,Ipiv,RWORK,TAU,WORK,INFO)
	    
	number=n-K-L
	write(6,'(/,a30)')  "   gensvd: Leads decimation   "
	write(6,'(a30,i6)') "   gensvd: Dim of H1 and S1 : ", n
	write(6,'(a30,i6)') "   gensvd: Rank of H1:        ", K
	write(6,'(a30,i6)') "   gensvd: Rank of (H1,S1):   ", L
	write(6,'(a30,i6)') 
     &   "   gensvd: Decimated states:  ", number
	if (side_rank .eq. "N") then
	 write(6,'(a35)')   
     &    "   gensvd: Decimation from the left"
	else
	 write(6,'(a36)') 
     &    "   gensvd: Decimation from the right"
	endif
	 
        If (side_rank .EQ. 'N') Then
	 Call ZGEMM('N','N',N,N,N,(1.0D0,0.D0),U,N,H1E,N,
     &    (0.D0,0.D0),VH,N)
	 Call ZGEMM('N','N',N,N,N,(1.0D0,0.D0),VMat,N,S1E,N,
     &    (0.D0,0.D0),VS,N)	    
        Else
	  Call ZGEMM ('C','C',N,N,N,(1.0D0,0.D0),H1E,N,
     &     U,N,(0.D0,0.D0),VH,N)
	  Call ZGEMM ('C','C',N,N,N,(1.0D0,0.D0),S1E,N,
     &     VMat,N,(0.D0,0.D0),VS,N)
     
     
        EndIf

        deallocate(S1E,H1E,U,VMat)
        deallocate(work,rwork,tau)
	
        End Subroutine GENSVD
	  
        
