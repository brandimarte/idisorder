C $Id: identify.f, v 1.1 03/29/2001 $
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

C *********************************************************************
C Subroutines for identifing degenerate eigenchannel
C Copyright by S.Sanvito, 2001
C
C Send comments/suggestions/bug-reports to:
C
C   Stefano Sanvito: sanvitos@tcd.ie
C *********************************************************************

C *************************** INPUT ***********************************
C integer      n           : Dimension Hamiltonian
C integer      nn          : nn=2*n
C integer      z(nn)       : exp(ik)
C integer      chan(nn)    : Vector with labels for the channels 
C *************************** OUTPUT **********************************
C integer      chan(nn)    : Vector with labels for the channels 
C integer      mul(n)      : Multeplicity degenerate eigenvalues
C *************************** DESCRIPTION *****************************
C chan is the vector containing the label of all the channels
C from 1 to ochan I've the open channel
C the deg matrix contains the degeneration information
C *********************************************************************

        subroutine identify(z,nn,chan,n,ochan,mul,deg)

        implicit none 

        integer ochan,nn,n,i,stepc,j,deg,count,count2

        complex*16 z(nn)

        integer chan(n),bchan(n),newchan(ochan+1),mul(n)
	real*8 tol

	parameter(tol=1.d-5)

c Initialize parameter

        do i=1,ochan+1
         newchan(i)=1
        enddo

        do i=1,n
         mul(i)=0
         bchan(i)=chan(i)
        enddo

        deg=0
        stepc=1
        count=1
        count2=1

        do j=1,ochan
         if(chan(j).ne.0) then
          do i=1,ochan

           if(chan(i).ne.0) then
            if(cdabs(z(chan(i))-z(bchan(j))).le.tol) then
             newchan(stepc)=chan(i)
             chan(i)=0
             stepc=stepc+1
            endif
           endif

          enddo
         endif
        enddo


        do i=1,ochan
         chan(i)=newchan(i)
        enddo


        do i=1,ochan
         if((cdabs(z(newchan(i+1))-z(newchan(i))).le.tol) .AND. 
     &      (i+1 .LE. ochan)) then
          count=count+1	 
          deg=1 
         else
          mul(count2)=count
          count2=count2+1
          count=1
         endif
        enddo

      return
      end
