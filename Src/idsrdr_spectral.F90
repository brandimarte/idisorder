!  *******************************************************************  !
!  I-Disorder Fortran Code                                              !
!                                                                       !
!  Written by Alexandre Reily Rocha and Pedro Brandimarte, 2007-2013    !
!                                                                       !
!  Copyright (c), All Rights Reserved                                   !
!                                                                       !
!  This program is free software. You can redistribute it and/or        !
!  modify it under the terms of the GNU General Public License          !
!  (version 3 or later) as published by the Free Software Foundation    !
!  <http://fsf.org/>.                                                   !
!                                                                       !
!  This program is distributed in the hope that it will be useful, but  !
!  WITHOUT ANY WARRANTY, without even the implied warranty of           !
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU     !
!  General Public License for more details (file 'LICENSE_GPL'          !
!  distributed along with this program or at                            !
!  <http://www.gnu.org/licenses/gpl.html>).                             !
!  *******************************************************************  !
!                        MODULE idsrdr_spectral                         !
!  *******************************************************************  !
!  Description: computes the spectral function and the density of       !
!  states (DOS) of dynamic units (units where electron-phonon           !
!  interaction is considered).                                          !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *******************************************************************  !

MODULE idsrdr_spectral

!
!   Modules
!
  use parallel,        only: 
  use idsrdr_options,  only: 
  use idsrdr_engrid,   only: 
  use idsrdr_ephcoupl, only: 
  use idsrdr_units,    only: 
  use idsrdr_green,    only: 

  implicit none

  real(8), dimension(:,:,:), allocatable :: spctrl ! spectral function
  real(8), dimension(:,:,:), allocatable :: dos ! density of states

  PUBLIC  :: spectralinit, spectral, writespectral, freespectral
  PRIVATE ! default is private


CONTAINS


!  *******************************************************************  !
!                             spectralinit                              !
!  *******************************************************************  !
!  Description: allocate spectral and density of states arrays.         !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer nspin               : Number of spin components              !
!  integer NTenerg_div         : Number of energy grid points per node  !
!  integer neph                : Number of units with e-ph interaction  !
!  *******************************************************************  !
  subroutine spectralinit

!
!   Modules
!
    use idsrdr_options,  only: nspin
    use idsrdr_engrid,   only: NTenerg_div
    use idsrdr_ephcoupl, only: neph

!   Allocate spectral and density of states arrays.
    allocate (spctrl(NTenerg_div,nspin,neph))
    allocate (dos(NTenerg_div,nspin,neph))


  end subroutine spectralinit


!  *******************************************************************  !
!                               spectral                                !
!  *******************************************************************  !
!  Description: computes the spectral function and the density of       !
!  states (DOS).                                                        !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode                       : True if it is the I/O node    !
!  integer nunits                       : Total number of units         !
!  integer neph                         : Number of units with e-ph     !
!                                         interaction                   !
!  integer norbDyn(neph)                : Number of orbitals from       !
!                                         dynamic atoms                 !
!  integer idxF(neph)                   : First dynamic atom orbital    !
!  integer idxL(neph)                   : Last dynamic atom orbital     !
!  integer unit_type(nunits+2)          : Units types                   !
!  integer unitdimensions(ntypeunits+2) : Units number of orbitals      !
!  integer ephIndic(ntypeunits+2)       : E-ph interaction indicator    !
!  TYPE(unitS) Sunits(ntypeunits+2)%S(unitdimensions,unitdimensions) :  !
!                                         [real*8] Units overlap        !
!  TYPE(green) Gr_nn(neph)%G(unitdimensions,unitdimensions) :           !
!                                                  [complex] G^r_{n,n}  !
!  ****************************** INPUT ******************************  !
!  integer ispin                        : Spin component index          !
!  real*8 ienergy                       : Energy grid index             !
!  *******************************************************************  !
  subroutine spectral (ienergy, ispin)

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_options,  only: nunits
    use idsrdr_ephcoupl, only: neph, norbDyn, idxF, idxL
    use idsrdr_units,    only: unit_type, unitdimensions,               &
                               ephIndic, Sunits
    use idsrdr_green,    only: Gr_nn

!   Input variables.
    integer, intent(in) :: ienergy, ispin

!   Local variables.
    integer :: I, J, K
    complex(8), dimension(:,:), allocatable :: Scp, Aux
    external :: zsymm

    if (IOnode) write (6,'(a)', advance='no')                           &
            '      computing spectral function... '

    J = 1 ! e-ph store indexing

    do I = 2,nunits ! over units from left to right

       if (ephIndic(unit_type(I)) == 1) then

!         Allocate auxiliary matrix.
          allocate (Scp(norbDyn(J),norbDyn(J)))
          allocate (Aux(norbDyn(J),norbDyn(J)))

!         Copy dynamic orbitals part of overlap matrix.
          Scp = Sunits(unit_type(I))%S(idxF(J):idxL(J),idxF(J):idxL(J))

!         (Aux = GrMM * Saux)
          call zsymm ('R', 'L', norbDyn(J), norbDyn(J), (1.D0,0.D0),    &
                      Scp, norbDyn(J), Gr_nn(J)%G, norbDyn(J),          &
                      (0.D0,0.D0), Aux, norbDyn)

          spctrl(ienergy,ispin,J) = 0.D0
          dos(ienergy,ispin,J) = 0.D0
          do K = 1,norbDyn(J)
             spctrl(ienergy,ispin,J) = spctrl(ienergy,ispin,J)          &
                                       - DIMAG(Gr_nn(J)%G(K,K))
             dos(ienergy,ispin,J) = dos(ienergy,ispin,J)                &
                                    - DIMAG(Aux(K,K))
          enddo

!         Free memory.
          deallocate (Scp)
          deallocate (Aux)

          J = J + 1

       endif

    enddo


    if (IOnode) write(6,'(a)') " ok!"


  end subroutine spectral


!  *******************************************************************  !
!                             writespectral                             !
!  *******************************************************************  !
!  Description: writes calculated spectral function and density of      !
!  states to output files ('slabel.SPCTR' and 'slabel.DOS').            !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer Node                : Actual node (MPI_Comm_rank)            !
!  integer Nodes               : Total number of nodes (MPI_Comm_size)  !
!  integer nspin               : Number of spin components              !
!  integer label_length        : Length of system label                 !
!  character(label_length) slabel : System Label (for output files)     !
!  character(60) directory     : Working directory                      !
!  integer NTenerg_div         : Number of energy grid points per node  !
!  real*8 Ei(NTenerg_div)      : Energy grid points                     !
!  integer neph                : Number of units with e-ph interaction  !
!  *******************************************************************  !
  subroutine writespectral

!
!   Modules
!
    use parallel,        only: IOnode, Node, Nodes
    use idsrdr_options,  only: nspin, label_length, slabel, directory
    use idsrdr_engrid,   only: NTenerg_div, Ei
    use idsrdr_ephcoupl, only: neph

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    integer :: J, n, e, s, iuSpc, iuDos
    real(8), parameter :: pi = 3.1415926535897932384626433832795028841D0
    real(8), dimension(:), allocatable :: buffEn
    real(8), dimension(:,:), allocatable :: buffSpc, buffDos
    character(len=10) :: suffix
    character(len=label_length+70) :: fnSpc, fnDos
    character(len=label_length+70), external :: paste
    character(len=10), external :: pasbias2
    external :: io_assign, io_close
#ifdef MPI
    integer :: MPIerror
    integer, dimension(MPI_Status_Size) :: MPIstatus
#endif

    do J = 1,neph ! over units with e-ph intereaction

       if (IOnode) then

!         Set file's names and open it.
          write (suffix,'(i3)') J
          suffix = pasbias2 (suffix, '.SPCTR')
          suffix = paste ('_', suffix)
          fnSpc = paste (slabel, suffix)
          fnSpc = paste (directory, fnSpc)
          write (suffix,'(i3)') J
          suffix = pasbias2 (suffix, '.DOS')
          suffix = paste ('_', suffix)
          fnDos = paste (slabel, suffix)
          fnDos = paste (directory, fnDos)
          call io_assign (iuSpc)
          open (iuSpc, file=fnSpc, form='formatted', status='unknown')
          call io_assign (iuDos)
          open (iuDos, file=fnDos, form='formatted', status='unknown')

          write (6, '(/,a,i3,a,a)', advance='no')                       &
               'Writing spectral function ', J, ' to: ', trim(fnSpc)
          write (6, '(/,a,i3,a,a)') 'Writing DOS ', J, ' to: ',         &
               trim(fnDos)

!         Allocate buffers arrays.
          allocate (buffEn(NTenerg_div))
          allocate (buffSpc(NTenerg_div,nspin))
          allocate (buffDos(NTenerg_div,nspin))

       endif

!      Write to the output files (energies in eV (from CODATA - 2012)).
#ifdef MPI
       do n = 0,Nodes-1
          if (Node == n .and. Node == 0) then
#endif
             do e = 1,NTenerg_div
                write (iuSpc, '(/,e17.8e3)', advance='no')              &
                     Ei(e) * 13.60569253D0 !eV
                write (iuDos, '(/,e17.8e3)', advance='no')              &
                     Ei(e) * 13.60569253D0 !eV
                do s = 1,nspin
                   write (iuSpc, '(e17.8e3)', advance='no')             &
                        spctrl(e,s,J) / (13.60569253D0 * pi) !eV
                   write (iuDos, '(e17.8e3)', advance='no')             &
                        dos(e,s,J) / (13.60569253D0 * pi) !eV
                enddo
             enddo
#ifdef MPI
          elseif (Node == n) then
             call MPI_Send (Ei, NTenerg_div, MPI_Double_Precision,      &
                            0, 1, MPI_Comm_world, MPIerror)
             call MPI_Send (spctrl(1,1,J), NTenerg_div*nspin,           &
                            MPI_Double_Precision, 0, 2,                 &
                            MPI_Comm_world, MPIerror)
             call MPI_Send (dos(1,1,J), NTenerg_div*nspin,              &
                            MPI_Double_Precision, 0, 3,                 &
                            MPI_Comm_world, MPIerror)
          elseif (Node == 0) then
             call MPI_Recv (buffEn, NTenerg_div, MPI_Double_Precision,  &
                            n, 1, MPI_Comm_world, MPIstatus, MPIerror)
             call MPI_Recv (buffSpc, NTenerg_div*nspin,                 &
                            MPI_Double_Precision, n, 2,                 &
                            MPI_Comm_world, MPIstatus, MPIerror)
             call MPI_Recv (buffDos, NTenerg_div*nspin,                 &
                            MPI_Double_Precision, n, 3,                 &
                            MPI_Comm_world, MPIstatus, MPIerror)
          endif
          if (n /= 0) then
             call MPI_Barrier (MPI_Comm_world, MPIerror)
             if (Node == 0) then
                do e = 1,NTenerg_div
                   write (iuSpc, '(/,e17.8e3)', advance='no')           &
                        buffEn(e) * 13.60569253D0 !eV
                   write (iuDos, '(/,e17.8e3)', advance='no')           &
                        buffEn(e) * 13.60569253D0 !eV
                   do s = 1,nspin
                      write (iuSpc, '(e17.8e3)', advance='no')          &
                           buffSpc(e,s) / (13.60569253D0 * pi) !eV
                      write (iuDos, '(e17.8e3)', advance='no')          &
                           buffDos(e,s) / (13.60569253D0 * pi) !eV
                   enddo
                enddo
             endif
          endif
       enddo
#endif

!      Close files and free buffers memory.
       if (Node.eq.0) then

          call io_close (iuSpc)
          call io_close (iuDos)

          deallocate (buffEn)
          deallocate (buffSpc)
          deallocate (buffDos)

       endif

    enddo ! do J = 1,neph


  end subroutine writespectral


!  *******************************************************************  !
!                             freespectral                              !
!  *******************************************************************  !
!  Description: free allocated vectors.                                 !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *******************************************************************  !
  subroutine freespectral

!   Free memory.
    deallocate (spctrl)
    deallocate (dos)


  end subroutine freespectral


!  *******************************************************************  !


END MODULE idsrdr_spectral

