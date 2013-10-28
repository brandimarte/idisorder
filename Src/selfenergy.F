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
	SUBROUTINE SELFENERGY(SIDE,side_rank,number,N2,Ei,H0,H1,S0,S1,
     &		    Sigma,Q,INFO,nrchan)

C *****************************************************************
C Calculates the self-energies of the right and left hand-side leads
C
C Written by Alexandre Reily Rocha, June 2003 
C Computational Spintronics Group
C Trinity College Dublin
C e-mail: rochaa@tcd.ie
C ***************************** HISTORY ***********************************
C Original version:	June 2003
C ***************************** INPUT ***************************** 
C character SIDE		  : 'L' (left) or 'R' (right) depending 
C				    on the self-energy calculated
C character side_rank		  : 'L' (left) or 'R' (right) depending
C				    on the direction of decimation
C integer number		  : number of decimated states
C integer N2		          : Dimension of the basis orbitals on
C                                   the lead
C complex*16 Ei            	  : Energy

C integer ISPIN			  : value of spin
C complex*8 H0(NL,NL,NSPIN)	  : Hamiltonian of the lead lead
C complex*8 H1(N2,N2,NSPIN)	  : Coupling Matrix
C complex*8 S0(N2,N2)	          : on-site Overlap Matrix
C complex*8 S1(N2,N2)   	  : neighbour Overlap Matrix
C complex*8 Q(N2,N2)		  : decimation matrix
C integer INFO			  : if INFO=0, self-energy calculated 
C				    with success, NOT 0 otherwise
C integer nrchan		  : number of open channels
C **************************** OUTPUT *******************************
C complex*16 Sigma                : Calculated self-energy
C *************************** CONSTANTS *****************************
C real*8 T                        : Electronic temperature in Ry
C **************************** AUXILIARIES **************************
 
	IMPLICIT NONE

	INCLUDE "const.h"

        DOUBLE COMPLEX, PARAMETER :: smear=(0.D0,0.0D0)

        INTEGER :: N2,number,N3,nrchan,nlchan,I,J,INFO

	DOUBLE COMPLEX :: Ei

	CHARACTER(LEN=1) :: SIDE,side_rank

	DOUBLE COMPLEX, DIMENSION (:,:), ALLOCATABLE :: Gr
	
	DOUBLE COMPLEX, DIMENSION (:,:), ALLOCATABLE :: 
     .   T1_aux,T1_dag,T1,H0_aux,
     .   H1_aux,H1_dag_aux,h0corr,Gr_2,Gr_square

        DOUBLE COMPLEX, DIMENSION (N2,N2) :: Sigma,Q,S0,S1,H0,H1 
	
	DOUBLE COMPLEX, DIMENSION (:,:), ALLOCATABLE :: Gr_1
	
	INTEGER, DIMENSION (:), ALLOCATABLE :: IPIV

	DOUBLE COMPLEX, DIMENSION (:), ALLOCATABLE :: WRK

	INTEGER, DIMENSION (:), ALLOCATABLE :: IPIV2

	DOUBLE COMPLEX, DIMENSION (:), ALLOCATABLE :: WRK2	

	external :: DECIMATE_LEADS, LEADS, MATRIXMULT
	
        allocate(T1_aux(N2,N2),T1_dag(N2,N2),T1(N2,N2),H0_aux(N2,N2),
     .   H1_aux(N2,N2),H1_dag_aux(N2,N2),h0corr(N2,N2),Gr_2(N2,N2),
     .   Gr_square(N2,N2))  

        allocate(Gr(2*(N2-number),2*(N2-number)),
     .   Gr_1(N2-number,N2-number))

        allocate(IPIV(n2-number),IPIV2(n2),WRK2(N2**2),
     .   WRK((N2-number)**2))
      
	N3=N2-number
	
	H0_aux=H0(:,:)-Ei*S0
	Call ZGEMM ('C','N',N2,N2,N2,(1.0D0,0.D0),Q,N2,H0_aux,N2,
     &   (0.D0,0.D0),T1_aux,N2)
       	Call ZGEMM ('N','N',N2,N2,N2,(1.0D0,0.D0),T1_aux,N2,Q,N2,
     &   (0.D0,0.D0),H0_aux,N2)
             
	Gr=(0.D0,0.D0)
	nrchan=0
	nlchan=0

	If (side_rank .NE. '0') Then
	 If (side_rank .EQ. 'N') Then
	 
	  T1_aux=H1(:,:)-Ei*S1
	  Call ZGEMM ('C','N',N2,N2,N2,(1.0D0,0.D0),Q,N2,T1_aux,N2,
     &    (0.D0,0.D0),H1_aux,N2)
     
	  DO I=1,N2
	   DO J=1,N2
	    T1_aux(I,J)=DCONJG(H1(J,I))-Ei*DCONJG(S1(J,I))
	   ENDDO
	  ENDDO
       	  Call ZGEMM('N','N',N2,N2,N2,(1.0D0,0.D0),T1_aux,N2,Q,N2,
     &    (0.D0,0.D0),H1_dag_aux,N2)
	 Else
	  T1_aux=H1(:,:)-Ei*S1
	  Call ZGEMM ('N','N',N2,N2,N2,(1.0D0,0.D0),T1_aux,N2,Q,N2,
     &    (0.D0,0.D0),H1_aux,N2)
     
	  DO I=1,N2
	   DO J=1,N2
	    T1_aux(I,J)=DCONJG(H1(J,I))-Ei*DCONJG(S1(J,I))
	   ENDDO
	  ENDDO
       	  Call ZGEMM('C','N',N2,N2,N2,(1.0D0,0.D0),Q,N2,T1_aux,N2,
     &    (0.D0,0.D0),H1_dag_aux,N2)
	 EndIf

         Call DECIMATE_LEADS(side_rank,side,N2,number,H0_aux,H1_aux,
     &    H1_dag_aux,h0corr)

     	Else
	 H0_aux=H0-Ei*S0
	 H1_aux=H1-Ei*S1
	 h0corr=0.D0
	 DO I=1,N2
	  DO J=1,N2
	   H1_dag_aux(I,J)=DCONJG(H1(J,I))-Ei*DCONJG(S1(J,I))
	  ENDDO
	 ENDDO
	EndIf 

	Call LEADS(SIDE,N3,2*N3,4*N3,Ei,H0_aux(number+1:N2,number+1:n2),
     &   H1_aux(number+1:N2,number+1:n2),
     &   H1_dag_aux(number+1:N2,number+1:n2),Gr,nrchan,nlchan,info)
        

	If (info .EQ. 0) Then      
	If (side .eq. 'L') Then
	 If (side_rank .EQ. 'N') Then
	  T1(1:n3,:)=H1_aux(number+1:N2,:)
	  T1_dag(:,1:n3)= H1_dag_aux(:,number+1:n2)
	  
          Gr_1=Gr(1:N3,1:N3)
	  
	  Call MATRIXMULT(Gr_1,T1(1:n3,:),T1_aux(1:n3,:),N3,N3,N2)
C          CALL ZGEMM('N','N',N3,N2,N3,(1.d0,0.d0),Gr_1,N3,T1(1:N3,:),N3,
C     &      (0.d0,0.d0),T1_aux(1:N3,:),N3)
	  Call MATRIXMULT(T1_dag(:,1:n3),T1_aux(1:n3,:),Sigma,N2,N3,N2)
C          CALL ZGEMM('N','N',N2,N2,N3,(1.d0,0.d0),T1_dag(:,1:N3),N2,
C     &       T1_aux(1:N3,:),N3,(0.d0,0.d0),Sigma,N2)


	  Gr_2=-h0corr-Sigma
          Call ZGETRF(n2,n2,Gr_2,n2,IPIV2,INFO)
	  Call ZGETRI(n2,Gr_2,n2,IPIV2,WRK2,n2**2,INFO)	   

          IF (IMAG(smear).GT. 1.D-7) THEN
           Call ZGEMM('N','N',n2,n2,n2,(1.D0,0.D0),Gr_2,n2,Gr_2,n2,
     &      (0.D0,0.D0),Gr_square,n2)
           Gr_2=Gr_2 - smear*Gr_square
          ENDIF

	  T1(:,:)=H1(:,:)-Ei*S1
       	  Call ZGEMM('N','C',N2,N2,N2,(1.0D0,0.D0),T1,N2,Q,N2,
     &    (0.D0,0.D0),T1_aux,N2)
       	  Call ZGEMM('C','N',N2,N2,N2,(1.0D0,0.D0),Q,N2,T1_aux,N2,
     &    (0.D0,0.D0),T1,N2)     
     	  
	  Do i=1,n2
	   Do j=1,n2
	    T1_dag(I,J)= DCONJG(H1(J,I))-Ei*DCONJG(S1(J,I))
           EndDo
	  EndDo
       	  Call ZGEMM('N','N',N2,N2,N2,(1.0D0,0.D0),T1_dag,N2,Q,N2,
     &    (0.D0,0.D0),T1_aux,N2)
       	  Call ZGEMM('N','N',N2,N2,N2,(1.0D0,0.D0),Q,N2,T1_aux,N2,
     &    (0.D0,0.D0),T1_dag,N2)     
	  
  	  Call MATRIXMULT(Gr_2,T1,T1_aux,N2,N2,N2)
C          CALL ZGEMM('N','N',N2,N2,N2,(1.d0,0.d0),Gr_2,N2,T1,N2,
C     &     (0.d0,0.d0),T1_aux,N2)
	  Call MATRIXMULT(T1_dag,T1_aux,Sigma,N2,N2,N2)
C          CALL ZGEMM('N','N',N2,N2,N2,(1.d0,0.d0),T1_dag,N2,T1_aux,N2,
C     &     (0.d0,0.d0),Sigma,N2)
          
	 Else
	 
	  T1(1:n3,:)=H1(number+1:n2,:)-Ei*S1(number+1:n2,:)
	  Do I=1,N2
	   Do J=1,N3
	    T1_dag(I,J)=DCONJG(H1(number+J,I))-Ei*DCONJG(S1(number+J,I))
	   EndDo
	  EndDo
	
     	  Gr_1=Gr(1:N3,1:N3)

	  If (side_rank .NE. '0') Then	  
           Call ZGETRF(n3,n3,Gr_1,n3,IPIV,INFO)
	   Call ZGETRI(n3,Gr_1,n3,IPIV,WRK,n3**2,INFO)
	   Gr_1=Gr_1-h0corr(number+1:n2,number+1:n2)
           Call ZGETRF(n3,n3,Gr_1,n3,IPIV,INFO)
	   Call ZGETRI(n3,Gr_1,n3,IPIV,WRK,n3**2,INFO)
	  EndIf

          IF (IMAG(smear).GT. 1.D-7) THEN
           Call ZGEMM('N','N',n3,n3,n3,(1.D0,0.D0),Gr_1,n3,Gr_1,n3,
     &      (0.D0,0.D0),Gr_square(1:n3,1:n3),n3)
           Gr_1 = Gr_1 - smear*Gr_square(1:n3,1:n3)
          ENDIF
                                               
	  Call MATRIXMULT(Gr_1,T1(1:n3,:),T1_aux(1:n3,:),N3,N3,N2)
C          CALL ZGEMM('N','N',N3,N2,N3,(1.d0,0.d0),Gr_1,N3,T1(1:n3,:),N3,
C     &     (0.d0,0.d0),T1_aux(1:n3,:),N3)

	  Call MATRIXMULT(T1_dag(:,1:n3),T1_aux(1:n3,:),Sigma,N2,N3,N2)
C          CALL ZGEMM('N','N',N2,N2,N3,(1.d0,0.d0),T1_dag(:,1:n3),N2,
C     &     T1_aux(1:n3,:),N3,(0.d0,0.d0),Sigma,N2)	 
 
	 EndIf
	
	Else
	 If (side_rank .EQ. 'C') Then
	  T1(:,1:n3)=H1_aux(:,number+1:n2)
	  T1_dag(1:n3,:)= H1_dag_aux(number+1:N2,:)
	  
          Gr_1=Gr(n3+1:2*N3,n3+1:2*N3)
	  
	  Call MATRIXMULT(Gr_1,T1_dag(1:n3,:),T1_aux(1:n3,:),N3,N3,N2)
	  Call MATRIXMULT(T1(:,1:n3),T1_aux(1:n3,:),Sigma,N2,N3,N2)
	  
	  Gr_2=-h0corr-Sigma
          Call ZGETRF(n2,n2,Gr_2,n2,IPIV2,INFO)
	  Call ZGETRI(n2,Gr_2,n2,IPIV2,WRK2,n2**2,INFO)	   

          IF (IMAG(smear).GT. 1.D-7) THEN
           Call ZGEMM('N','N',n2,n2,n2,(1.D0,0.D0),Gr_2,n2,Gr_2,n2,
     &      (0.D0,0.D0),Gr_square,n2)
           Gr_2=Gr_2 - smear*Gr_square
          ENDIF
                                               
	  T1(:,:)=H1(:,:)-Ei*S1
       	  Call ZGEMM('N','N',N2,N2,N2,(1.0D0,0.D0),T1,N2,Q,N2,
     &    (0.D0,0.D0),T1_aux,N2)
       	  Call ZGEMM('N','N',N2,N2,N2,(1.0D0,0.D0),Q,N2,T1_aux,N2,
     &    (0.D0,0.D0),T1,N2)     
     	  
	  Do i=1,n2
	   Do j=1,n2
	    T1_dag(I,J)= H1(J,I)-Ei*S1(J,I)
           EndDo
	  EndDo
       	  Call ZGEMM('N','C',N2,N2,N2,(1.0D0,0.D0),T1_dag,N2,Q,N2,
     &    (0.D0,0.D0),T1_aux,N2)
       	  Call ZGEMM('C','N',N2,N2,N2,(1.0D0,0.D0),Q,N2,T1_aux,N2,
     &    (0.D0,0.D0),T1_dag,N2)     
	  
  	  Call MATRIXMULT(Gr_2,T1_dag,T1_aux,N2,N2,N2)
	  Call MATRIXMULT(T1,T1_aux,Sigma,N2,N2,N2)
	 Else
	 
	  T1(:,1:N3)=H1(:,number+1:n2)-Ei*S1(:,number+1:n2)
	  Do I=1,N3
	   Do J=1,N2
	    T1_dag(I,J)=DCONJG(H1(J,number+I))-Ei*DCONJG(S1(J,number+I))
	   EndDo
	  EndDo
	
     	  Gr_1=Gr(n3+1:2*N3,n3+1:2*N3)
	  
	  If (side_rank .NE. '0') Then
           Call ZGETRF(n3,n3,Gr_1,n3,IPIV,INFO)
	   Call ZGETRI(n3,Gr_1,n3,IPIV,WRK,n3**2,INFO)
	   Gr_1=Gr_1-h0corr(number+1:n2,number+1:n2)
           Call ZGETRF(n3,n3,Gr_1,n3,IPIV,INFO)
	   Call ZGETRI(n3,Gr_1,n3,IPIV,WRK,n3**2,INFO)
	  EndIf

          IF (IMAG(smear).GT. 1.D-7) THEN
           Call ZGEMM('N','N',n3,n3,n3,(1.D0,0.D0),Gr_1,n3,Gr_1,n3,
     &      (0.D0,0.D0),Gr_square(1:n3,1:n3),n3)
           Gr_1 = Gr_1 - smear*Gr_square(1:n3,1:n3)
          ENDIF
                                               
	  Call MATRIXMULT(Gr_1,T1_dag(1:n3,:),T1_aux(1:n3,:),N3,N3,N2)
	  Call MATRIXMULT(T1(:,1:n3),T1_aux(1:n3,:),Sigma,N2,N3,N2)   
	 EndIf
	EndIf	  
	Else
	 RETURN	
	EndIf

        deallocate(T1_aux,T1_dag,T1,H0_aux,
     .   H1_aux,H1_dag_aux,h0corr,Gr_2,Gr_square)  

        deallocate(Gr,Gr_1)

        deallocate(IPIV,IPIV2,WRK2,WRK)
	
	END SUBROUTINE SELFENERGY

