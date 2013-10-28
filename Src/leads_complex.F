C $Id: leads.f, v 1.1 03/29/2001 $
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

C ********************************************************************
C Subroutines for calculating the Surface Green function
C Copyright by S.Sanvito, 2001
C
C Send comments/suggestions/bug-reports to:
C
C   Stefano Sanvito: sanvitos@tcd.ie
C *************************** HISTORY ********************************
C Original Version		2001
C 
C Altered to include non-orthogonal basis set	January 2003
C Alexandre Reily Rocha, rochaa@tcd.ie
C 
C Altered to include in the Smeagol package
C Alexandre Reily Rocha  rochaa@tcd.ie		June 2003
C ********************************************************************
C The surface Green functions are calculated by using the technique of
C S.Sanvito et al. PRB, 69, 11936 (1999), see also S.Sanvito PhD Thesis
C
C *************************** DESCRIPTION ****************************
C This code calculates the Greens function g_o and group velocites
C v_g of a lead with n degrees of freedom, described by
C Hamiltonian  Matrices H_o and H_1. The corresponding Schrodinger 
C equation is
c E Psi_i = H_o Psi_i + H_1 Psi_1 + H_-1 Psi_-1
C where H_-1(i,j) = H_1^*(j,i). H_o, H_1 are nxn matrices,
C whereas g_0 is a (nn x nn) matrix, where nn=2n.
C
C The user chooses the energy E & the code calculates 
C the phase-factors  z = e^{ika}
C It organises these into open & closed channels
C and right & left going states, calculates the group velocities and
C finally the Green's function in a semi-infinite lead g_o.
C
C The inputs are the array dimensions n, nn=2n, nnn=2nn, energy E &
C the Hamiltonian matrices h0 and h1. These are destroyed on-exit.
C
C On-exit A contains g_o and vg contains the group velocities,
C ho, h1 the left and right going eigenvectors respectively and
C rdual, ldual the right and left-going dual vectors.
C ie rdual is the inverse of h0 and ldual is the inverse of h1
C
C The number of open (closed) right-going channels is nrchan (nrcchan)
C The number of open (closed) left-going channels is nlchan (nlcchan)
C These satisfy nrcchan=nlcchan, nrchan=nlchan, nlchan+nlcchan=n
C
C Elements 1 to nrchan (nrchan+1 to n) of vg(i) contain the group
C velocities of the open (closed) right-going channels.
C Elements n+1 to n+nlchan (n+nlchan+1 to nn) of vg(i) contain the group
C velocities of the open (closed) left-going channels.
C Columns 1 to nrchan (nrchan+1 to n) of h0(i,j) contain the eigenvectors
C of the open (closed) right-going channels.
C Columns 1 to nlchan (nlchan+1 to n) of h1(i,j) contain the eigenvectors
C of the open (closed) left-going channels.
C Rows 1 to nrchan (nrchan+1 to n) of rdual(i,j) contain the eigenvectors
C of the open (closed) right-going channels.
C Rows 1 to nlchan (nlchan+1 to n) of ldual(i,j) contain the eigenvectors
C of the open (closed) left-going channels.
C
C *************************** INPUT **********************************
C integer      n            : Dimension Hamiltonian
C real*8       e            : Energy
C complex*16   h0(n,n)      : On-site Hamiltonian
C complex*16   h1(n,n)      : Coupling-site Hamiltonian
C
C  Alterations made by Alexandre Reily Rocha to Introduce
C  non-orthogonal basis set (10/09/2002)
C complex*16   S_0(n,n)     : On-site Coupling matrix
C complex*16   S_1(n,n)     : Coupling-site Overlap matrix
C
C *************************** OUTPUT **********************************
C integer      nrchan	    : number-left going channel
C integer      nlchan	    : number right-going channel
C integer      nrcchan	    : number-left decaying channel
C integer      nlcchan	    : number right-decaying channel
C complex*16   rdual(n,n)   : right-going dual vectors (k|
C complex*16   ldual(n,n)   : left-going dual vectors (k-bar|
C complex*16   h0(n,n)      : vectors |k> 
C complex*16   h1(n,n)      : vectors |k-bar> 
C complex*16   A(2*n,2*n)   : Green function (see below)
C complex*16   vg(2*n)      : Group velocity
C complex*16   z(2*n)       : e^(ik)
C complex*16   h2(n,n)      : Operator V 
C
C *************************** AUXILIARY *******************************
C complex*16   hh(n,n)
C complex*16   hi(n,n)
C complex*16   hj(n,n)
C complex*16   hl(n,n)
C complex*16   hk(n,n)
C complex*16   aux(n,n)
C *********************************************************************

      subroutine leads(side,n,nn,nnn,Ei,h_2,v02,v02_dag,
     +              A,nrchan,nlchan,info)

      IMPLICIT NONE 

      DOUBLE PRECISION, PARAMETER :: 
CH1CondMin=1.0D-11,H1NormMax=3.0D+3,
     +                          InvTolCond=5.0D-8,InvTolNorm=5.0D+8,
     +				XInvTol=1.0D-6,
     +				SVDTol=1.0D-7,
     +				zero=1.D-6
      CHARACTER(LEN=1) :: side
      INTEGER :: n,i,j,k,l,nn,info,nchan,nrchan,nlchan,nlcchan,nrcchan,
     + nnn

      DOUBLE COMPLEX, DIMENSION (:,:), ALLOCATABLE :: vec,vec1
      DOUBLE COMPLEX, DIMENSION (:), ALLOCATABLE   :: wrk
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: wrk1
      DOUBLE COMPLEX  :: zi,Ei

      INTEGER, DIMENSION (:), ALLOCATABLE   :: rchan,lchan

      DOUBLE COMPLEX, DIMENSION (:,:), ALLOCATABLE :: 
     +   hj,hi,hh,h1,h1_dag,h0,h2,SJ_aux
      DOUBLE COMPLEX, DIMENSION (:), ALLOCATABLE   :: wk1

      DOUBLE COMPLEX, DIMENSION (n,n) :: h_2,v02,v02_dag,rdual,ldual
      DOUBLE COMPLEX, DIMENSION (nn,nn) :: A
      DOUBLE COMPLEX, DIMENSION	(nn)  :: vg,z

      INTEGER :: deg,degbar,s,dimo,start,fin
      INTEGER, DIMENSION (n)	      :: mul,mulbar

      INTEGER, DIMENSION (n)	      :: IPIV,IPIV2
      DOUBLE PRECISION :: ANORM,RCOND,anorm_sav,rcond_sav
      double precision, dimension (n) :: FERR,BERR
      double complex, dimension (:,:), ALLOCATABLE :: Ident
      DOUBLE COMPLEX, DIMENSION (:), ALLOCATABLE :: WORKEIG
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: 
     &  RCONDE,RCONDV,SCALE
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: RWORK
      INTEGER :: ILO,IHI
     
      PARAMETER(zi=(0.0d0,1.0d0))

      ALLOCATE(vec(nn,nn),vec1(nn,nn),wrk(nnn),wrk1(nnn),rchan(n),
     &  lchan(n))
      ALLOCATE(SJ_aux(n,n),wk1(n))
      ALLOCATE(hj(n,n),hi(n,n),hh(n,n),h1(n,n),h1_dag(n,n),h0(n,n),
     & h2(n,n))

      ALLOCATE(WORKEIG(2*nn**2+2*nn),Ident(n,n),
     &  RCONDE(nn),RCONDV(nn),SCALE(nn),RWORK(5*n))

      Ident=0.0
      do i=1,n
       Ident(i,i)=1.0D0
      enddo

     
      info=0

      h0=h_2
      h1=v02
      h1_dag=v02_dag
      
      vg=(0.0d0,0.0d0)
      z=(0.0d0,0.0d0)
      A=(0.0d0,0.0d0)
      rchan=0
      lchan=0

C########################################################################
C SOLVE THE DISPERSION AS K=K(E)
C########################################################################

      h2 = h1
      hj = h1
      h0 = -h0

      DO I=1,n
       WRK1(I)=0.D0
       Do J=1,n
        WRK1(I)=WRK1(I)+CDABS(h2(J,I))
       EndDo
      ENDDO
      ANORM=MAXVAL(WRK1(1:n))
      CALL ZGETRF(n,n,h2,n,IPIV,INFO)
      CALL ZGECON('1',n,h2,n,ANORM,RCOND,WRK(1:2*n),WRK1(1:2*n),INFO)
      anorm_sav = anorm
      rcond_sav = rcond
#ifdef SMEAGOLDEBUG
       WRITE(116,'(a8,2d16.8)') 'Energy :',DREAL(Ei),DIMAG(Ei)
       WRITE(116,'(a38,d11.3,d11.3,i4)') 
     &  'H1^{-1} : 1-Norm, Cond number, Info: ',ANORM,RCOND,INFO
#endif
      If ((RCOND .gt. InvTolCond) .AND. (ANORM .lt. InvTolNorm)) Then
       hi=h2
       IPIV2=IPIV
       CALL ZGETRI(n,h2,n,IPIV,WORKEIG(1:2*n**2+n),2*n**2+n,INFO)

       if (rcond .lt. XInvTol) then
        call ZGERFS('N',n,n,hj,n,hi,n,IPIV2,Ident,n,h2,n,
     &  FERR, BERR, WRK(1:2*n),WRK1(1:n),INFO)      
       endif
       hi = h2
      Else
#ifdef SMEAGOLDEBUG
       WRITE(116,'(a16)') 'Calling SVD H1:'
#endif
       h2=h1
       CALL ZGESVD('A','A',n,n,h2,n,FERR,rdual,n,hi,n,
     .  WORKEIG(1:2*n**2+n),2*n**2+n,RWORK,INFO)
       ldual=hi
#ifdef SMEAGOLDEBUG
       WRITE(116,'(a4,2d16.7)') 'SVD ',FERR(1),FERR(n)
#endif
       Do I=1,n
        Do J=1,n
	 If ((FERR(J) .GT. ABS(FERR(1))*SVDTol).AND.(J.NE.n)) Then
	  h2(I,J)=DCONJG(hi(J,I))/FERR(J)
	  hj(I,J)=rdual(I,J)*FERR(J)
	 Else
	  h2(I,J)=DCONJG(hi(J,I))/(ABS(FERR(1))*SVDTol)
	  hj(I,J)=rdual(I,J)*(ABS(FERR(1))*SVDTol)
	 EndIf
	EndDo
       EndDo
       
       CALL ZGEMM('N','C',n,n,n,(1.D0,0.D0),h2,n,rdual,n,
     &  (0.D0,0.D0),hi,n)
       CALL ZGEMM('N','N',n,n,n,(1.D0,0.D0),hj,n,ldual,n,
     &  (0.D0,0.D0),h1,n)   
       
       if (DIMAG(Ei) .lt. zero) then
        Do I=1,n
         Do J=1,n
	  h1_dag(I,J)=DCONJG(h1(J,I))
	 EndDo
        EndDo        

       else
        CALL ZGESVD('A','A',n,n,h1_dag,n,FERR,rdual,n,ldual,n,
     .  WORKEIG(1:2*n**2+n),2*n**2+n,RWORK,INFO)
        
        Do I=1,n
         Do J=1,n
	  If ((FERR(J) .GT. ABS(FERR(1))*SVDTol).AND.(J.NE.N)) Then
	   h2(I,J)=rdual(I,J)*FERR(J)
	  Else
	   h2(I,J)=rdual(I,J)*ABS(FERR(1))*SVDTol
	  EndIf
	 EndDo
        EndDo
       
        CALL ZGEMM('N','N',n,n,n,(1.D0,0.D0),h2,n,ldual,n,
     &  (0.D0,0.D0),rdual,n)              
        
        h1_dag=rdual
       Endif
      EndIf 

C================================================= 
C After re-definition the matrices contain:
C 
C  Alterations made by Alexandre Reily Rocha to Introduce
C  non-orthogonal basis set (10/09/2002)
C
C  	h0 = E*S_0 - H_0 : H_0^E
C  	h1 = H_1	 : H_1^E
C  	hh = H_1-E*S_1   : H_1^E
C  	hi = (H_1-E*S_1)^{-1} : H_1^E^{-1}
C  	h2 = (H_1-E*S_1)^{-1} : H_1^E^{-1}
C=================================================

C==========================================================
C  set up the matrix A , which is defined as follows:
C
C      A = | (H_1^E^{-1}H_0^E)   -H_1^E^{-1}H_1^E^{-1} |
C          |          1                      0         |
C
C      where  H_{-1}(i,j) = H_{1}^{*}(j,i)
C
C  and solve the eigenvalues problem   
C
C==========================================================

      do j=1,n,1
       l=j+n
       A(l,j) = (1.0d0,0.0d0)
       do i=1,n,1
        do k=1,n,1
         A(i,j)= A(i,j) + (hi(i,k) * h0(k,j))
         A(i,l)= A(i,l) - hi(i,k) * h1_dag(k,j)
        enddo
       enddo
      enddo

      Call ZGEEVX('P','N','V','V',nn,A,nn,z,vec,nn,vec1,nn,ILO,
     & IHI,SCALE,ANORM,rconde,rcondv,workeig,2*nn**2+2*nn,wrk1,INFO)

       do i=1,nn
	if (dreal(sum(vec1(1:n,i)*dconjg(vec1(1:n,i)))) .lt.
     +   dreal(sum(vec1(n+1:nn,i)*dconjg(vec1(n+1:nn,i))))) then
	  vec1(1:n,i)=vec1(n+1:nn,i)
	 endif
	enddo
	
C==========================================================
C z   = eigenvalues e^{ika}
C vec = eigenvectors
C==========================================================

C########################################################################
C
C  Now calculate the group velocities , which are normalised
C  expectation values of the group velocity operator
C  Vg = i ( H_{1} z - H_{-1}/z )
C
C  The group velocities are, vg = <vec| Vg |vec> / <vec|vec>
C
C########################################################################

      do i=1,nn,1
       h2= z(i)*h1 - h1_dag/z(i)

       do j=1,n
	h2(j,1) = sum(h2(j,1:n)*vec1(1:n,i))
       end do

       vg(i) = zi*sum(dconjg(vec1(1:n,i))*h2(1:n,1))/
     *  (sum(dconjg(vec1(1:n,i))*vec1(1:n,i)))

      enddo 

C########################################################################
C  Now sort out the states into right and left going ones: 
C  
C	Left going open    =  vg<0 & |z|=1
C	Left going closed  =  |z|>1 
C	Right going open   =  vg>0 & |z|=1
C	Right going closed =  |z|<1
C
C  Check also symmetry
C
C########################################################################

      nchan=0
      nlchan=0
      nrchan=0
      nlcchan=0
      nrcchan=0

      do i=1,nn,1
       if(dabs(1.0d0-cdabs(z(i))).le.5.0d-5) then 
        nchan=nchan+1 
        if(dreal(vg(i)).gt.0.0d0) then
         nrchan=nrchan+1
	 rchan(nrchan)=i
        elseif(dreal(vg(i)).lt.0.0d0) then
         nlchan=nlchan+1
         lchan(nlchan)=i
        else
         print*,'vg(',i,') = 0'
        endif
       elseif(abs(z(i)).lt.(1.0d0-5.0d-5)) then
        nrcchan=nrcchan+1
	IF ((n+1-nrcchan) .EQ. 0) THEN
	 PRINT*, n,i,nrcchan,'righ'
         print*, z(i),'final' 
         print*, z
	 INFO = 1
	 RETURN 
	ENDIF  
        rchan((n+1-nrcchan))=i
       elseif(abs(z(i)).gt. (1.0d0+5.0d-5)) then 
	 nlcchan=nlcchan+1
	 IF ((n+1-nlcchan) .EQ. 0) THEN
	  PRINT*, n,i,nlcchan
          print*, z(i),'final'
          print*, z
	  INFO = 1
	  RETURN 
	 ENDIF
         lchan((n+1-nlcchan))=i
       else
        print*,'error during sorting in leads'
       endif
      enddo

#ifdef SMEAGOLDEBUG      
      WRITE(116,'(a13,4i5)') 'left, right : ',
     & nlchan,nlcchan,nrchan,nrcchan
#endif	
      if((nlchan+nrchan).ne.nchan) then
       print*,'left + right open channels not equal to total'
       INFO = 1
       RETURN 
      endif

      if((nlchan+nrchan+nlcchan+nrcchan).ne.nn) then
       print*,'total channels not equal to number of group velocities'
       INFO = 1
       RETURN
      endif

      if(nlchan+nlcchan.ne.nrchan+nrcchan) then
       print*,'no. of left channels .ne. number of right channels'
       print*,'Norm :',anorm_sav,' Cond num : ',rcond_sav
       INFO = 1
       RETURN 
      endif

C########################################################################
C  NOW SOLVE THE PROBLEM OF DEGENERATE EIGENCHANNELS
C   
C  The code diagonalizes the current operator over each subspace
C  of degenerate vectors.
C
C  Subroutine identify stores sequencially the degenerate eigenvectors
C  in rchan (lchan).
C  mul (mulbar) = the multeplicity of each eigenvalue
C  
C  Subspace parameters
C  dimo=Subspace dimension
C  start=label of the first vector belonging to the subspace
C
C########################################################################

	IF((nrchan.gt.0).and.(nlchan.gt.0)) then

	 call identify(z,nn,rchan,n,nrchan,mul,deg)

	 IF(deg.eq.1) then

#ifdef SMEAGOLDEBUG	 
	  WRITE(116,*) 'degeneracy right =',deg
#endif
	  start=1
	  dimo=0

	  DO s=1,n

	   IF ((mul(s).ne.1).and.(mul(s).ne.0)) THEN

	    start=start+dimo
	    dimo=mul(s)
	    fin=start+dimo-1
C==========================================================
C  |k> vectors
C==========================================================

	    h2= z(rchan(start))*h1 - h1_dag/z(rchan(start))

	    DO k=start,fin
	     DO j=1,n
	      wk1(j) = sum(h2(j,1:n)*vec1(1:n,rchan(k)))
	     ENDDO
	
	     DO j=start,fin

	      h0(j,k)= zi*sum(dconjg(vec1(1:n,rchan(j)))*wk1(1:n))

	      SJ_aux(j,k)=SUM(DCONJG(vec1(1:n,rchan(j)))*
     .          vec1(1:n,rchan(k)))

	     ENDDO
	    ENDDO

	    call diagonal(n,h0,SJ_aux,dimo,dimo+dimo,start)

C==========================================================
C  new |k> eigenvectors are stored in hi and than copied in vec1
C==========================================================

	    hh=(0.d0,0.d0)

C==========================================================
C  Build the diagonal basis |k-tilda>
C==========================================================

	    do i=start,fin
	     do j=1,n
	      do l=start,fin
	       hh(j,i)=hh(j,i)+vec1(j,rchan(l))*h0(l,i)*
     *                    DSQRT(DREAL(SJ_aux(i,i)))
	      enddo
	     enddo
	    enddo
	   
C==========================================================
C  Store the new basis and the new group velocity
C==========================================================

	    do i=start,fin
	     do j=1,n
	      vec1(j,rchan(i))=hh(j,i)	
	     enddo
	    enddo

	    DO i=start,fin

	     h0(i,i)=(0.d0,0.d0)

	     DO j=1,n
	      wk1(j) = sum(h2(j,1:n)*vec1(1:n,rchan(i)))
	     ENDDO

	     h0(i,i)= zi*sum(dconjg(vec1(1:n,rchan(i)))*wk1(1:n))
	     
	     vg(rchan(i))=h0(i,i)

	    ENDDO

	   elseif(mul(s).eq.1) then

	    start=start+dimo
	    dimo=mul(s)

	   endif

	  enddo

	 elseif(deg.eq.0) then 
	 endif

C==========================================================
C  |k-bar> vectors
C==========================================================

	call identify(z,nn,lchan,n,nlchan,mulbar,degbar)

	if(degbar.eq.1) then
#ifdef SMEAGOLDEBUG	
	 WRITE(116,*) 'degeneracy left =',degbar
#endif
	 start=1
	 dimo=0

	  do s=1,n

	   if((mulbar(s).ne.1).and.(mulbar(s).ne.0)) then

	    start=start+dimo
	    dimo=mulbar(s)
	    fin=start+dimo-1

	    h2 = z(lchan(start))*h1 - h1_dag/z(lchan(start))
	
	    DO k=start,fin
	     DO j=1,n
	      wk1(j) = sum(h2(j,1:n)*vec1(1:n,lchan(k)))
	     ENDDO
	
	     DO j=start,fin

	      h0(j,k)= zi*sum(dconjg(vec1(1:n,lchan(j)))*wk1(1:n))
	      
	      SJ_aux(j,k)=SUM(DCONJG(vec1(1:n,lchan(j)))*
     .          vec1(1:n,lchan(k)))

	     ENDDO
	    ENDDO          

	    call diagonal(n,h0,SJ_aux,dimo,dimo+dimo,start)

C==========================================================
C  New |k-bar> eigenvectors are stored in hi and than copied in vec1
C==========================================================


	    hh = (0.d0,0.d0)

C==========================================================
C  Build the diagonal basis |k-bar-tilda>
C==========================================================

	   do i=start,fin
	    do j=1,n
	     do l=start,fin
	      hh(j,i)=hh(j,i)+vec1(j,lchan(l))*h0(l,i)*
     *                   DSQRT(DREAL(SJ_aux(i,i)))
	     enddo
	    enddo
	   enddo

C==========================================================
C  Store the new basis and the new group velocity
C==========================================================

	   do i=start,fin
	    do j=1,n
 	     vec1(j,lchan(i))=hh(j,i)	
	    enddo
	   enddo


	   DO i=start,fin

	    h0(i,i)=(0.d0,0.d0)

	    DO j=1,n
	     wk1(j) = sum(h2(j,1:n)*vec1(1:n,lchan(i)))
	    ENDDO
	
	    h0(i,i)= zi*sum(dconjg(vec1(1:n,lchan(i)))*wk1(1:n))
	    
	    vg(lchan(i))=h0(i,i)

	   ENDDO
	   	     
	   elseif(mulbar(s).eq.1) then

	    start=start+dimo
	    dimo=mulbar(s)

	   endif

	  enddo

	 elseif(degbar.eq.0) then 
	 endif

	 elseif((nrchan.le.0).and.(nlchan.le.0)) then
	 endif
C==========================================================
C  Store the right-going, right-decaying eigenvectors h0 and hi
C  Store the left-going, left-decaying eigenvectors h1 and h2
C  hi and h2 are simply used as workspace when computing the inverse
C  of h0 and h1
C==========================================================

        do i=1,n
         h0(:,i)=vec1(1:n,rchan(i))
         hi(:,i)=vec1(1:n,rchan(i))
         h1(:,i)=vec1(1:n,lchan(i))
         h2(:,i)=vec1(1:n,lchan(i))
        enddo 

C==========================================================
C  Now construct the dual vector (vec| such that (vec|vec> = I, the unit 
C  matrix therefore if all of |vec> are written in C , D=C^{-1} , DC = I 
C  and so (vec_{k}| is the k'th row of D.
C==========================================================

      hj=hi
      DO I=1,n
       WRK1(I)=0.D0
       Do J=1,n
        WRK1(I)=WRK1(I)+CDABS(hi(J,I))
       EndDo
      ENDDO
      ANORM=MAXVAL(WRK1(1:n))

      CALL ZGETRF(n,n,hi,n,IPIV,INFO)
      CALL ZGECON('1',n,hi,n,ANORM,RCOND,WRK(1:2*n),WRK1(1:2*n),INFO)
#ifdef SMEAGOLDEBUG
      WRITE(116,'(a8,d12.4,d12.4,i4,a2)')
     & 'rdual : ',ANORM,RCOND,INFO,side
#endif
      If ((RCOND .gt. InvTolCond) .AND. (ANORM .lt. InvTolNorm)) Then
       rdual=hi
       IPIV2=IPIV
       CALL ZGETRI(n,rdual,n,IPIV,WORKEIG(1:2*n**2+n),2*n**2+n,INFO)
       if (rcond .lt. XInvTol) then
        call ZGERFS('N',n,n,hj,n,hi,n,IPIV2,Ident,n,rdual,n,
     &	FERR, BERR, WRK(1:2*n),WRK1(1:n),INFO)      
       endif   		     
      Else
#ifdef SMEAGOLDEBUG
       WRITE(116,'(a18)') 'Calling SVD rdual:'
#endif
       CALL ZGESVD('A','A',n,n,hj,n,FERR,hi,n,rdual,n,
     .  WORKEIG(1:2*n**2+n),2*n**2+n,RWORK,INFO)
     
       Do I=1,n
        Do J=1,n
	 If ((FERR(J) .GT. ABS(FERR(1))*SVDTol).AND.(J.NE.N)) Then
	  hj(I,J)=DCONJG(rdual(J,I))/FERR(J)
	 Else
	  hj(I,J)=DCONJG(rdual(J,I))/(ABS(FERR(1))*SVDTol)
	 EndIf
	EndDo
       EndDo
       
       Do I=1,n
        Do J=1,n
	 rdual(I,J)=0.D0
         Do L=1,n
	  rdual(I,J)=rdual(I,J)+hj(I,L)*DCONJG(hi(J,L))
	 EndDo
	EndDo
       EndDo
      EndIf 
      
      hj=h2
      Do I=1,n
       WRK1(I)=0.D0
       Do J=1,n
        WRK1(I)=WRK1(I)+CDABS(h2(J,I))
       EndDo
      ENDDO
      ANORM=MAXVAL(WRK1(1:n))

      CALL ZGETRF(n,n,h2,n,IPIV,INFO)
      CALL ZGECON('1',n,h2,n,ANORM,RCOND,WRK(1:2*n),WRK1(1:2*n),INFO)
#ifdef SMEAGOLDEBUG
      WRITE(116,'(a8,d12.4,d12.4,i4,a2)')
     & 'ldual : ',ANORM,RCOND,INFO,side
#endif
      If ((RCOND .gt. InvTolCond) .AND. (ANORM .lt. InvTolNorm)) Then
       ldual=h2
       IPIV2=IPIV
       CALL ZGETRI(n,ldual,n,IPIV,WORKEIG(1:2*n**2+n),2*n**2+n,INFO)
       if (rcond .lt. XInvTol) then
        Call ZGERFS('N',n,n,hj,n,h2,n,IPIV2,Ident,n,ldual,n,
     &	 FERR, BERR, WRK(1:2*n),WRK1(1:n),INFO)      
       endif
      Else
#ifdef SMEAGOLDEBUG
       WRITE(116,'(a18)') 'Calling SVD ldual:'
#endif
       CALL ZGESVD('A','A',n,n,hj,n,FERR,h2,n,ldual,n,
     .  WORKEIG(1:2*n**2+n),2*n**2+n,RWORK,INFO)
           
       Do I=1,n
        Do J=1,n
	 If ((FERR(J) .GT. (ABS(FERR(1))*SVDTol)).AND.(J.NE.N)) Then
	  hj(I,J)=DCONJG(ldual(J,I))/FERR(J)
	 Else
	  hj(I,J)=DCONJG(ldual(J,I))/(ABS(FERR(1))*SVDTol)
	 EndIf
	EndDo
       EndDo
       
       Do I=1,n
        Do J=1,n
	 ldual(I,J)=0.D0
         Do L=1,n
	  ldual(I,J)=ldual(I,J)+hj(I,L)*DCONJG(h2(J,L))
	 EndDo
	EndDo
       EndDo
      EndIf 	


C########################################################################
C
C  Now calculate the Green's Function , which is defined as follows:
C    g_{o} = | Right lead	g(x,x) rg |
C	     | g(x,x) lg	Left lead |
C  where:
C      Right lead  is the surface Green fn. for the right lead,
C      Left lead   is the surface Green fn. for the left lead,
C      g(x,x) lg   is the infinite left going Green fn,
C      g(x,x) rg   is the infinite right going Green fn.
C
C########################################################################

C==========================================================
C  First store the operator V in the workspace array  h2 and hh
C
C  h2 is used as output
C
C==========================================================

      do i=1,n
       do j=1,n
        hi(i,j)=(0.0d0,0.0d0)
        do l=1,n
         hi(i,j)=hi(i,j)+vec1(i,rchan(l))*rdual(l,j)/z(rchan(l))
     +               -vec1(i,lchan(l))*ldual(l,j)/z(lchan(l))
        enddo
       enddo
      enddo
	
      DO I=1,n
       WRK1(I)=0.D0
       Do J=1,n
        WRK1(I)=WRK1(I)+CDABS(hi(J,I))
       EndDo
      ENDDO
      ANORM=MAXVAL(WRK1(1:n))
      
      CALL ZGEMM('N','N',n,n,n,(1.D0,0.D0),h1_dag,n,hi,n,
     &  (0.D0,0.D0),hh,n)
     
C      do i=1,n
C       do j=1,n
C        h2(i,j)=(0.0d0,0.0d0)
C        do l=1,n
C         h2(i,j)=h2(i,j)+h1_dag(i,l)*hi(l,j)
C        enddo
C       enddo
C      enddo

C      hh=h2

C==========================================================
C  Then place the inverse of V in hi:
C==========================================================

      hj=hh
      DO I=1,n
       WRK1(I)=0.D0
       Do J=1,n
        WRK1(I)=WRK1(I)+CDABS(hh(J,I))
       EndDo
      ENDDO
      ANORM=MAXVAL(WRK1(1:n))
      CALL ZGETRF(n,n,hh,n,IPIV,INFO)
      CALL ZGECON('1',n,hh,n,ANORM,RCOND,WRK(1:2*n),WRK1(1:2*n),INFO)
#ifdef SMEAGOLDEBUG
       WRITE(116,'(a9,d12.4,d12.4,i4)') 'V^{-1} : ',ANORM,RCOND,INFO
#endif
      If ((RCOND .gt. InvTolCond) .AND. (ANORM .lt. InvTolNorm)) Then
       hi=hh
       IPIV2=IPIV
       CALL ZGETRI(n,hh,n,IPIV,WORKEIG(1:2*n**2+n),2*n**2+n,INFO)
       if (rcond .lt. XInvTol) then
        call ZGERFS('N',n,n,hj,n,hi,n,IPIV2,Ident,n,hh,n,
     &	 FERR, BERR, WRK(1:2*n),WRK1(1:n),INFO)      
       endif   		     
       hi=hh
      Else
       CALL ZGESVD('A','A',n,n,hj,n,FERR,hh,n,hi,n,
     .  WORKEIG(1:2*n**2+n),2*n**2+n,RWORK,INFO)
       
       Do I=1,n
        Do J=1,n
	 If ((FERR(J) .GT. ABS(FERR(1))*SVDTol).AND.(J.NE.N)) Then
	  hj(I,J)=DCONJG(hi(J,I))/FERR(J)
	 Else
	  hj(I,J)=DCONJG(hi(J,I))/(ABS(FERR(1))*SVDTol)
	 EndIf
	EndDo
       EndDo
       
       CALL ZGEMM('N','C',n,n,n,(1.D0,0.D0),hj,n,hh,n,
     &  (0.D0,0.D0),hi,n)
     
      EndIf 

C==========================================================
C  Now construct Green functions and place in A
C==========================================================

C==========================================================
C  Matrices common to the construction of the left and right lead:
C	hh is |k>e^{ik}(k|,
C	h1 is |k-bar>e^{-ik-bar}(k-bar|.
C==========================================================

      do i=1,n
       do j=1,n
       hh(i,j)=(0.0d0,0.0d0)
        do l=1,n
         hh(i,j)=hh(i,j)+
     +   vec1(i,rchan(l))*rdual(l,j)*z(rchan(l))
        enddo
       enddo
      enddo


      do i=1,n
       do j=1,n
       h1(i,j)=(0.0d0,0.0d0)
        do l=1,n
         h1(i,j)=h1(i,j)+
     +   vec1(i,lchan(l))*ldual(l,j)/z(lchan(l))
        enddo
       enddo	
      enddo

C==========================================================
C  Left lead
C ==========================================================
      If (side .eq. 'L') Then

       CALL ZGEMM('N','N',n,n,n,(-1.D0,0.D0),h1,n,hh,n,
     &  (0.D0,0.D0),h0,n)
     
C       do i=1,n
C        do j=1,n
C        h0(i,j)=(0.d0,0.d0)
C         do l=1,n
C          h0(i,j)=h0(i,j)-h1(i,l)*hh(l,j)
C         enddo
C        enddo
C       enddo

       do j=1,n
        h0(j,j)=(1.d0,0.d0)+h0(j,j)
       enddo

       CALL ZGEMM('N','N',n,n,n,(1.D0,0.D0),h0,n,hi,n,
     &  (0.D0,0.D0),a(1:n,1:n),n)
      
C       do i=1,n
C        do j=1,n
C	 a(i,j)=(0.d0,0.d0)
C         do l=1,n
C          a(i,j)=a(i,j)+h0(i,l)*hi(l,j)
C         enddo
C        enddo
C       enddo
      EndIf
 
C==========================================================
C  Right lead
C==========================================================
      If (side .eq. 'R') then

       CALL ZGEMM('N','N',n,n,n,(-1.D0,0.D0),hh,n,h1,n,
     &  (0.D0,0.D0),h0,n)
     
C       do i=1,n
C        do j=1,n
C        h0(i,j)=(0.d0,0.d0)
C         do l=1,n
C          h0(i,j)=h0(i,j)-hh(i,l)*h1(l,j)
C         enddo
C        enddo
C       enddo

       do j=1,n
        h0(j,j)=(1.d0,0.d0)+h0(j,j)
       enddo
      
       CALL ZGEMM('N','N',n,n,n,(1.D0,0.D0),h0,n,hi,n,
     &  (0.D0,0.D0),a(n+1:2*n,n+1:2*n),n)
C       do i=1,n
C        do j=1,n
C	 a(i+n,j+n)=(0.d0,0.d0)
C         do l=1,n
C          a(i+n,j+n)=a(i+n,j+n)+h0(i,l)*hi(l,j)
C         enddo
C        enddo
C       enddo
      EndIf
      
       info=0
C      Else
C       info=1
C       print*, ANORM,RCOND
C      EndIf

      DEALLOCATE(vec,vec1,wrk,wrk1,rchan,lchan)
      DEALLOCATE(SJ_aux,wk1)
      DEALLOCATE(hj,hi,hh,h1,h1_dag,h0,
     & h2)

      DEALLOCATE(WORKEIG,Ident,
     &  RCONDE,RCONDV,SCALE,RWORK)
 
      return
      end
