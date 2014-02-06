C $Id: diagonal.f, v 1.1 04/03/2001 $
! 
! Copyright (c) Trinity College Dublin:
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
      subroutine diagonal(n,C,D,m,mm,start)

C ********************************************************************
C Subroutines for diagonalizing the current operator in the degenerate
C space
C Copyright by S.Sanvito, 2001
C
C Send comments/suggestions/bug-reports to:
C
C   Stefano Sanvito: sanvitos@tcd.ie
C ************************** HISTORY *********************************
C Original Version		2001
C
C Altered to include non-orthogonal basis set by
C  Alexandre Reily Rocha, 2003
C email: rochaa@tcd.ie
C *************************** INPUT **********************************
C integer      n           : dimension Hamiltonian 
C complex*16   C(n,n)      : Current operator
C *************************** OUTPUT *********************************
C complex*16   C(n,n)      : Transformation Matrix
C ********************************************************************

      implicit none 

	integer n,i,j,m,in,mm,start

	DOUBLE COMPLEX, DIMENSION (n,n) :: C,D
	COMPLEX*16 B(m,m),E(m,m),lambda(m),dwrk(mm)

         real*8 dwrk1(3*m-2)	

	do i=1,m
	 do j=1,m
	  B(i,j)=(0.d0,0.d0)
	  E(i,j)=(0.d0,0.d0)
	 enddo
	enddo

	do i=1,m
	 do j=1,m
	  B(i,j)=C(start+i-1,start+j-1)
	  E(i,j)=D(start+i-1,start+j-1)
	 enddo
	enddo
	
C       call zgeev('N','V',m,B,m,lambda,dvec,m,dvec1,m,dwrk,mm,dwrk1,in) 
	CALL ZHEGV(1,'V','U',m,B,m,E,m,lambda,dwrk,mm,dwrk1,IN)

c check se la matrice di trasformazione funziona
c check equazione agli autovalori

c	do i=1,m
c	do j=1,m
c	B(i,j)=C(start+i-1,start+j-1)
c	enddo
c	enddo

c	do i=1,m
c	do j=1,m
c	C(i,j)=(0.d0,0.d0)
c	enddo
c	enddo

c check orthogonality eigenvectors
c la subroutine mi da' una base ortonormale

	do i=1,n
	 do j=1,n
	  C(i,j)=(0.d0,0.d0)
	 enddo
	enddo

	do i=1,m
	 do j=1,m
	  C(start+i-1,start+j-1)=B(i,j)
	 enddo
	enddo

      return
      end
