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
      subroutine hsunits(trnspin, nspin, dnuoL, nsc, iter,
     .                   istep, iv, gamma, ik, temp,
     .                   H0L, H1L, S0L, S1L, slabeli)
C *******************************************************************
C Reads and calculates the Hamiltonians, overlaps and density matrices
C of the leads.
C Written by V. M. Garcia-Suarez.
C Departamento de Fisica
C Universidad de Oviedo
C e-mail: victor@condmat.uniovi.es
C ***************************** HISTORY *******************************
C Original version:	June 2003
C *****************************INPUT*********************************
C For simplicity, I list only the left lead variables
C integer trnspin              : True value of the spin
C integer nspin                : Number of spin components
C integer dnuoL                : Number of basis orbitals in the left 
C                                lead, including spin components
C integer dnuoR                : Number of basis orbitals in the right 
C                                lead, including spin components
C integer nsc(2)               : Number of unit cells along parallel
C                                directions
C integer iter                 : Scf iteration in SIESTA
C integer istep                : Molecular dynamics iteration
C integer iv                   : Bias potential iteration
C logical gamma                : Calculation with parallel k-points
C integer ik                   : k-point index
C real*8  temp                 : Electronic temperature
C *****************************OUTPUT********************************
C complex*8  H0L(dnuoL,dnuoL,nspin) : Hamiltonian in the unit cell of the
C                                left lead
C complex*8  H1L(dnuoL,dnuoL,nspin) : Conexion between unit cells in the
C                                left lead
C complex*8  S0L(dnuoL,dnuoL)    : Overlaps in the unit cell of the
C                                left lead
C complex*8  S1L(dnuoL,dnuoL)    : Overlaps between unit cells in the
C                                left lead
C complex*8  H0R(dnuoR,dnuoR,nspin) : Hamiltonian in the unit cell of the
C                                right lead
C complex*8  H1R(dnuoR,dnuoR,nspin) : Conexion between unit cells in the
C                                right lead
C complex*8  S0R(dnuoR,dnuoR)    : Overlaps in the unit cell of the
C                                right lead
C complex*8  S1R(dnuoR,dnuoR)    : Overlaps between unit cells in the
C                                right lead
C *******************************************************************

      implicit none
      
      integer
     .  trnspin, nspin, dnuoL, nsc(2), iter, istep, iv, ik
      double precision
     .  temp
      double precision 
     .  H0L(dnuoL,dnuoL,nspin), H1L(dnuoL,dnuoL,nspin),
     .  S0L(dnuoL,dnuoL), S1L(dnuoL,dnuoL)
      logical
     .  gamma
      character
     .  slabel*20, slabeli*30

C Internal variables
      integer
     .  nspinL, nscL(2), 
     .  iu, iu1, io, iuo, ind, j
      integer ::
     .  nuoL, noL, maxnhL
      integer, allocatable ::
     .  numhL(:), listhptrL(:), indxuoL(:),
     .  listhL(:)
      double precision
     .  efL, tempL
      double precision, allocatable ::
     .  xijL(:,:), SL(:), HL(:,:), foo(:)
      character
     .  paste*25
      
      external
     .  io_assign, io_close, hsl

      if (iter.eq.1 .and. istep.eq.0 .and. iv.eq.0 .and. ik.eq.1) then
C Read data
        call io_assign(iu1)
        open(iu1,file=paste(slabeli,'.DAT'),status='old')
        read(iu1,*) slabel, nuoL, nspinL, maxnhL, efL, tempL,
     .              nscL(1), nscL(2), noL
        write(6,*) slabel, nuoL, nspinL, maxnhL, efL, tempL,
     .              nscL(1), nscL(2), noL


c Allocate arrays
        allocate(numhL(nuoL),listhptrL(nuoL))
        allocate(indxuoL(noL))
        allocate(listhL(maxnhL))
        allocate(xijL(3,maxnhL))

        allocate(SL(maxnhL))
        allocate(HL(maxnhL,trnspin))
        allocate(foo(trnspin))

c Read data of the left lead
        do iuo = 1, nuoL
          read(iu1,*) numhL(iuo), listhptrL(iuo)
        enddo
        do io = 1, noL
          read(iu1,*) indxuoL(io)
        enddo
        do iuo = 1, nuoL
          do j = 1, numhL(iuo)
            ind = listhptrL(iuo) + j
            read(iu1,*) listhL(ind)
            if (.not.gamma) then
              read(iu1,*) xijL(1,ind), xijL(2,ind), xijL(3,ind)
            endif
          enddo
        enddo
      
        call io_close(iu1)

c Compare supercell of leads and EM
        if (nsc(1).ne.nscL(1) .or. nsc(2).ne.nscL(2)) then
          write(6,'(a)') 'ERROR: The left supercell along parallel'
          write(6,'(a)') 'directions is different of the supercell'
          write(6,'(a)') 'of the EM. Change the size of it'
          write(6,'(a,2i4)') 'nsc = ', nsc(1), nsc(2)
          write(6,'(a,2i4)') 'nscL = ', nscL(1), nscL(2)
        endif
          
c Verify if the spin at the leads and the EM is the same
        if (nspinL.ne.trnspin) then
          write(6,'(a)') 'ERROR: The spin at the left lead is not'
          write(6,'(a)') 'the same as in the extended molecule'
          stop
        endif

c Check the temperature
        if (tempL-temp.gt.1.d-4) then
          write(6,'(a)')'WARNING: The temperature at the left lead is'
          write(6,'(a)')'not the same as in the extended molecule'
        endif

c Initialize overlaps and Hamiltonians
        HL = 0.d0
        SL = 0.d0

c Read overlaps and Hamiltonians
        call io_assign(iu)
        open(iu,file=paste(slabel,'.HST'),status='old')
        do iuo = 1, nuoL
          do j = 1, numhL(iuo) 
            ind = listhptrL(iuo) + j
            read(iu,*) SL(ind)
            read(iu,*) foo
            HL(ind,:) = foo
          enddo
        enddo
        call io_close(iu)
        

        if(gamma) then 
c Calculate S0, S1, H0 and H1 (left and right)
          call hsl(trnspin, nspin, nuoL, dnuoL, noL, maxnhL,
     .             numhL, listhptrL, indxuoL, listhL, SL, HL,
     .             S0L, S1L, H0L, H1L)
        endif
          
      endif

c Calculate S0, S1, H0 and H1 (left and right) if not gamma point
      if(.not.gamma) then
       write(6,*) "calculation must be performed in the gamma point"
       stop
C        call hslk(trnspin, nspin, nuoL, dnuoL, noL, maxnhL, numhL,
C     .            listhptrL, indxuoL, listhL, xijL, kpoint, nsc,
C     .            SL, HL, S0L, S1L, H0L, H1L)
C        call hslk(trnspin, nspin, nuoR, dnuoR, noR, maxnhR, numhR,
C     .            listhptrR, indxuoR, listhR, xijR, kpoint, nsc,
C     .            SR, HR, S0R, S1R, H0R, H1R)
      endif

      deallocate(numhL,listhptrL,indxuoL,listhL,xijL,SL,HL,foo)

      end

