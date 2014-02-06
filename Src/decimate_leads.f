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
	SUBROUTINE DECIMATE_LEADS(side_rank,side_lead,n,number,H0E,
     &               H1E,H1E_dag,h0correc)

C *****************************************************************
C Calculates the decimated Hamiltonian for the leads and the correction
C
C Written by Alexandre Reily Rocha, June 2003
C Computational Spintronics Group
C Trinity College Dublin
C e-mail: rochaa@tcd.ie
C ***************************** HISTORY ***********************************
C Original version:	June 2003
C ***************************** INPUT *****************************
C character side_lead             : 'L' (left) or 'R' (right) depending
C                                   on the self-energy calculated
C character side_rank             : 'L' (left) or 'R' (right) depending
C                                   on the direction of decimation
C integer n                       : Dimension of the basis orbitals on
C				    the lead
C integer number                  : number of decimated states
C
C complex*8 H0E(n,n)              : Hamiltonian of the lead
C complex*8 H1E(n,n)              : Coupling Matrix
C complex*8 H1E_dag(n,n)          : Complex conjugate of the Coupling
C **************************** OUTPUT *******************************
C complex*8 H0E(n,n)              : Decimated Hamiltonian of the lead
C complex*8 H1E(n,n)              : Decimated Coupling Matrix
C complex*8 H1E_dag(n,n)          : Decimated Complex conjugate of the Coupling
C complex*8 h0correc(n,n)	  : Decimated correction to the last surface
C				    of the lead (due to finiteness of the lead)
C *******************************************************************

	
	IMPLICIT NONE
	
	CHARACTER(LEN=1) :: side_lead,side_rank

	INTEGER :: number,n

	DOUBLE COMPLEX, DIMENSION (n,n) :: H0E,H1E,H1E_dag,h0correc

	DOUBLE COMPLEX, DIMENSION (:,:), ALLOCATABLE ::aux_mat

        allocate(aux_mat(2*N,2*N))

	aux_mat(1:n,n+1:2*n)=H1E
	aux_mat(n+1:2*n,1:n)=H1E_dag
	h0correc=H0E
	
	IF (side_rank .EQ. 'C') THEN
	 aux_mat(1:n,1:n)=0.D0
	 aux_mat(n+1:2*n,n+1:2*n)=H0E
	 CALL DECIMATION_L(n,number,aux_mat,h0correc)
	 H0E=aux_mat(1:n,1:n)
	 H1E=aux_mat(1:n,n+1:2*n)
	 H1E_dag=aux_mat(n+1:2*n,1:n)
	 If (side_lead .EQ. 'L') Then
	  h0correc=aux_mat(n+1:2*n,n+1:2*n)-aux_mat(1:n,1:n)
	 EndIf
	ELSE
	 aux_mat(1:n,1:n)=H0E
	 aux_mat(n+1:2*n,n+1:2*n)=0.D0
	 CALL DECIMATION_R(n,number,aux_mat,h0correc)
	 H0E=aux_mat(n+1:2*n,n+1:2*n)
	 H1E=aux_mat(1:n,n+1:2*n)
	 H1E_dag=aux_mat(n+1:2*n,1:n)
	 If (side_lead .EQ. 'R') Then
	  h0correc=aux_mat(1:n,1:n)-aux_mat(n+1:2*n,n+1:2*n)
	 EndIf
	ENDIF

        deallocate(aux_mat)
	
	CONTAINS

	 SUBROUTINE DECIMATION_L(n,number,mat,h0corr)

	 IMPLICIT NONE
	 
	 INTEGER :: n,number
	 DOUBLE COMPLEX, DIMENSION (2*n,2*n) :: mat
	 DOUBLE COMPLEX, DIMENSION (n,n) :: h0corr
         DOUBLE COMPLEX, DIMENSION (:,:), ALLOCATABLE  :: aux
 
	 INTEGER :: I,J,L
	
         allocate(aux(n,n))
 
	 DO I=n+1,n+number
	  DO J=1,2*n
	   DO L=1,2*n
	    IF (((J .GT. n) .AND. (J .LE. I)) .OR. ((L .GT. n) 
     &       .AND. (L .LE. I))) THEN
	    ELSE
	     mat(J,L)=mat(J,L)-mat(J,I)*mat(I,L)/mat(I,I)
	    ENDIF
	   ENDDO
	  ENDDO
	 ENDDO

	 aux=mat(1:n,1:n)+h0corr
	 mat(1:n,1:n)=h0corr+mat(1:n,1:n)
	 h0corr=aux
	 	 
	 DO I=1,number
	  DO J=I+1,2*n
	   DO L=I+1,2*n
	    IF (((J .GT. n) .AND. (J .LE. n+number)) .OR. ((L .GT. n) 
     &       .AND. (L .LE. n+number))) THEN
            Else
	     mat(J,L)=mat(J,L)-mat(J,I)*mat(I,L)/mat(I,I)
	    Endif
	   ENDDO
	  ENDDO
	 ENDDO

         DEALLOCATE(aux)

       	 END SUBROUTINE DECIMATION_L
	 
	 SUBROUTINE DECIMATION_R(n,number,mat,h0corr)
	 
	 IMPLICIT NONE
	 
	 INTEGER :: n,number
	 DOUBLE COMPLEX, DIMENSION (2*n,2*n) :: mat
	 DOUBLE COMPLEX, DIMENSION (n,n) :: h0corr
	 DOUBLE COMPLEX, DIMENSION (:,:), ALLOCATABLE :: aux
         	 
	 INTEGER :: I,J,L
		 
         ALLOCATE(aux(n,n))

	 DO I=1,number
	  DO J=I+1,2*n
	   DO L=I+1,2*n
	    mat(J,L)=mat(J,L)-mat(J,I)*mat(I,L)/mat(I,I)
	   ENDDO
	  ENDDO
	 ENDDO

	 aux=mat(n+1:2*n,n+1:2*n)+h0corr
	 mat(n+1:2*n,n+1:2*n)=h0corr+mat(n+1:2*n,n+1:2*n)
	 h0corr=aux
	 
	 DO I=n+1,n+number
	  DO J=number+1,2*n
	   DO L=number+1,2*n
	    IF (((L .GT. n) .AND. (L .LE. I)) .OR. (
     &       (J .GT. n) .AND. (J .LE. I))) THEN
	    ELSE
	     mat(J,L)=mat(J,L)-mat(J,I)*mat(I,L)/mat(I,I)
	    ENDIF
	   ENDDO
	  ENDDO
	 ENDDO

         DEALLOCATE(aux)

       	 END SUBROUTINE DECIMATION_R

	END SUBROUTINE DECIMATE_LEADS

