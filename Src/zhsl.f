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
      subroutine zhsl(trnspin, nspin, nuo, dnuo, no, maxnh,
     .               numh, listhptr, indxuo, listh, S, H,
     .               S0, S1, H0, H1)
C *********************************************************************
C Subroutine which calculates S0, S1, H0 and H1 
C Written by V. M. Garcia-Suarez
C Departamento de Fisica
C Universidad de Oviedo
C e-mail: victor@condmat.uniovi.es
C ***************************** HISTORY *******************************
C Original version:	June 2003
C **************************** INPUT **********************************
C integer trnspin             : True value of the spin
C integer nspin               : Number of spin components
C integer nuo                 : Number of basis orbitals in unit cell
C integer dnuo                : Number of basis orbitals in unit cell
C                               including spin components
C integer no                  : Number of basis orbitals in supercell
C integer maxnh               : Maximum number of orbitals interacting
C integer numh(nuo)           : Number of nonzero elements of each row
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               Hamiltonian matrix
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C integer listh(maxnh)        : Nonzero Hamiltonian-matrix element
C                               column indexes for each matrix row
C real*8  S(maxnh)            : Overlap in sparse form
C real*8  H(maxnh,trnspin)    : Hamiltonian in sparse form
C **************************** OUTPUT *********************************
C complex*8  S0(dnuo,dnuo)    : Overlap in the unit cell
C complex*8  S1(dnuo,dnuo)    : Overlap that connects unit cells along
C                               the transport direction
C complex*8  H0(dnuo,dnuo,nspin) : Hamiltonian in the unit cell
C complex*8  H1(dnuo,dnuo,nspin) : Hamiltonian that connects unit cells
C                               along the transport direction
C *********************************************************************
C
C Modules
C
C      use sys 

      implicit none

      integer
     .  trnspin, nspin, nuo, dnuo, no, maxnh, numh(nuo),
     .  listhptr(nuo), indxuo(no), listh(maxnh)
      double precision
     .  S(maxnh), H(maxnh,trnspin)
      double complex
     .  S0(dnuo,dnuo), S1(dnuo,dnuo), H0(dnuo,dnuo,nspin),
     .  H1(dnuo,dnuo,nspin)

C  Internal variables .............................................
      integer
     .  iuo, juo, j, jo, ind
      double complex
     .  ii
      
C ....................

      ii = (0.d0,1.d0)

c Find S0, S1, H0 and H1
      S0 = 0.0d0
      S1 = 0.0d0
      H0 = 0.0d0
      H1 = 0.0d0
      if (trnspin.le.2) then
        do iuo = 1, nuo
          do j = 1, numh(iuo)
            ind = listhptr(iuo) + j
            jo = listh(ind)
            juo = indxuo(jo)
            if(jo.le.nuo) then
              S0(juo,iuo) = S(ind)
              H0(juo,iuo,:) = H(ind,:)
            else if(jo.gt.nuo.and.jo.le.2*nuo) then
              S1(iuo,juo) = S(ind)
              H1(iuo,juo,:) = H(ind,:)
            endif
          enddo
        enddo
      else 
        do iuo = 1, nuo
          do j = 1, numh(iuo)
            ind = listhptr(iuo) + j
            jo = listh(ind)
            juo = indxuo(jo)
            if(jo.le.nuo) then
              S0(juo,iuo) = S(ind)
              S0(nuo+juo,nuo+iuo) = S(ind)
              H0(juo,iuo,1) = H(ind,1)
              H0(juo,nuo+iuo,1) = H(ind,3) - ii*H(ind,4)
              H0(nuo+juo,nuo+iuo,1) = H(ind,2)
            else if(jo.gt.nuo.and.jo.le.2*nuo) then
              S1(iuo,juo) = S(ind)
              S1(nuo+iuo,nuo+juo) = S(ind)
              H1(iuo,juo,1) = H(ind,1)
              H1(iuo,nuo+juo,1) = H(ind,3) - ii*H(ind,4)
              H1(nuo+iuo,nuo+juo,1) = H(ind,2)
            endif
          enddo
        enddo
      endif

c Hermiticity of H (non-collinear spin)
      if(trnspin.gt.2) then
        do iuo = 1, nuo
          do juo = 1, nuo
            H0(nuo+juo,iuo,1) = dconjg(H0(iuo,nuo+juo,1))
            H1(nuo+juo,iuo,1) = dconjg(H1(iuo,nuo+juo,1))
          enddo
        enddo
      endif

      return
      end

