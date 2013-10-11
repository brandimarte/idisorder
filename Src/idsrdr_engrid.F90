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
!                         MODULE idsrdr_engrid                          !
!  *******************************************************************  !
!  Description: create and distribute to nodes the energy grid.         !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *******************************************************************  !

MODULE idsrdr_engrid

!
!   Modules
!
  use parallel,        only: 
  use idsrdr_options,  only: 
  use idsrdr_leads,    only: 

  implicit none
  
  PUBLIC ! default is public
  PRIVATE :: energygrid

  integer :: NTenerg_div ! number of energy grid points per node

  real(8), allocatable, dimension (:) :: Ei ! energy grid points
  real(8), allocatable, dimension (:) :: gweight ! energy grid weights


CONTAINS


!  *******************************************************************  !
!                                engrid                                 !
!  *******************************************************************  !
!  Description: subroutine to create and to distribute the energy grid  !
!  to the nodes.                                                        !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer Nodes               : Total number of nodes (MPI_Comm_size)  !
!  integer NTenerg             : Number of transmission energy points   !
!  real*8 TEnergI              : Initial transmission energy            !
!  real*8 TEnergF              : Final transmission energy              !
!  real*8 temp                 : Electronic temperature                 !
!  real*8 EfLead               : Lead Fermi energy                      !
!  ***************************** OUTPUT ******************************  !
!  integer NTenerg_div         : Number of energy grid points per node  !
!  real*8 Ei(NTenerg_div)      : Energy grid points                     !
!  real*8 gweight(NTenerg_div) : Energy grid weights                    !
!  *******************************************************************  !
  subroutine engrid

!
!   Modules
!
    use parallel,        only: IOnode, Nodes
    use idsrdr_options,  only: NTenerg, TEnergI, TEnergF, temp
    use idsrdr_leads,    only: EfLead

    if (IOnode) write (6,'(/,a)') 'engrid: Computing energy grid...'

#ifdef MPI
    NTenerg_div = NTenerg / Nodes
    if (NTenerg_div == 0) NTenerg_div = 1
    NTenerg = Nodes * NTenerg_div ! redefine the total energy grid
#else
    NTenerg_div = NTenerg
#endif

!   Allocate the energy grid points and weights arrays.
    allocate (Ei(NTenerg_div), gweight(NTenerg_div))

!   Compute the energy grid.
    call energygrid (NTenerg, NTenerg_div, TEnergI, TEnergF,            &
                     temp, EfLead, Ei, gweight)
    if (IOnode) write(6,'(/,a)') 'engrid: done!'


  end subroutine engrid


!  *******************************************************************  !
!                              energygrid                               !
!  *******************************************************************  !
!  Description: computes the energy grid.                               !
!                                                                       !
!  Written by Alexandre Reily Rocha, --- 2007.                          !
!  Instituto de Fisica Teorica                                          !
!  Universidade Estadual de Sao Paulo                                   !
!  e-mail: reilya@ift.unesp.br                                          !
!  ***************************** HISTORY *****************************  !
!  Original version:    --- 2007                                        !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer Node                : Actual node (MPI_Comm_rank)            !
!  integer Nodes               : Total number of nodes (MPI_Comm_size)  !
!  character(14) integraltype  : Integration method                     !
!  ****************************** INPUT ******************************  !
!  integer NTenerg             : Number of transmission energy points   !
!  integer NTenerg_div         : Number of energy grid points per node  !
!  real*8 TEnergI              : Initial transmission energy            !
!  real*8 TEnergF              : Final transmission energy              !
!  real*8 temp                 : Electronic temperature                 !
!  real*8 EfLead               : Lead Fermi energy                      !
!  ***************************** OUTPUT ******************************  !
!  real*8 Ei(NTenerg_div)      : Energy grid points                     !
!  real*8 gweight(NTenerg_div) : Energy grid weights                    !
!  *******************************************************************  !
  subroutine energygrid (NTenerg, NTenerg_div, TEnergI, TEnergF,        &
                         temp, EfLead, Ei, gweight)

!
!   Modules
!
    use parallel,        only: IOnode, Node, Nodes
    use idsrdr_options,  only: integraltype

    include "mpif.h"

!   Input variables.
    integer, intent(in) :: NTenerg, NTenerg_div
    real(8), intent(in) :: TenergI, TenergF, temp, EfLead
    real(8), dimension (NTenerg_div), intent(out) :: Ei, gweight


!   Local variables.
    integer :: I
    real(8) :: arctanhinterm
    real(8), parameter :: arctanhinit= -0.999999999999D0
    real(8), parameter :: arctanhfim =  0.999999999999D0
    real(8), dimension (NTenerg) :: Energaux, gweightaux
    external :: gauleg
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

#ifdef IBM
    interface
       real*8 function atanh (%val(x))
         real*8 x
       end function atanh
    end interface
#endif

!   Initialize arrays.  
    Ei = 0.d0
    gweight = 0.d0

    IF (integraltype == 'Sympson') THEN

       do I = 1,NTenerg_div
          if (NTEnerg == 1) then
             Ei(I) = TEnergI
          else
             arctanhinterm = arctanhinit                                &
                  + DBLE(I-1+Node*NTenerg_div) / DBLE(NTenerg-1)        &
                  * DBLE(arctanhfim-arctanhinit)
             Ei(I) = 2.0D0 * temp * ATANH(arctanhinterm) + EfLead
          endif
       enddo
 
       do I = 1,NTenerg_div
          If (MOD(NTenerg,2) /= 0) Then
             if (((Node == 0).and.(I == 1)) .OR.                        &
                  ((Node == Nodes-1).and.(I == NTenerg_div))) then
                gweight(I) = 1.0d0/3.0d0                                &
                     * (arctanhfim-arctanhinit) / DBLE(NTenerg-1)
              elseif (MOD(I+Node*NTenerg_div,2) == 0) then
                gweight(I) = 4.0d0/3.0d0                                &
                     * (arctanhfim-arctanhinit) / DBLE(NTenerg-1)
             else
                gweight(I) = 2.0d0/3.0d0                                &
                     * (arctanhfim-arctanhinit) / DBLE(NTenerg-1) 
             endif
          Else

! Alex: Tomar cuidado, testar se NTenerg_div eh igual a 1.
             if ((Node == 0) .and. (I == 1)) then
                gweight(I) = 1.0d0/3.0d0                                &
                     * (arctanhfim-arctanhinit) / DBLE(NTenerg-1)
             elseif ((Node == Nodes-1) .and. (I == NTenerg_div-1)) then
                gweight(I) = (1.0d0/3.0d0+1.0d0/2.0d0)                  &
                     * (arctanhfim-arctanhinit) / DBLE(NTenerg-1) 
             elseif ((Node == Nodes-1) .and. (I == NTenerg_div)) then
                gweight(I) = 1.0d0/2.0d0                                &
                     * (arctanhfim-arctanhinit) / DBLE(NTenerg-1)
             elseif (MOD(I+Node*NTenerg_div,2) == 0) then
                gweight(I) = 4.0d0/3.0d0                                &
                     * (arctanhfim-arctanhinit) / DBLE(NTenerg-1)
             else
                gweight(I) = 2.0d0/3.0d0                                &
                     * (arctanhfim-arctanhinit) / DBLE(NTenerg-1)
             endif
          EndIf
       enddo

       print*, 0.5*sum(gweight)

    ELSEIF (integraltype .eq. 'Gauss') THEN
      
       if (Node == 0) then 
          call gauleg (arctanhinit, arctanhfim,                         &
                       Energaux, gweightaux, NTEnerg)

          write (6,*) "Precision :", 0.5D0*SUM(gweightaux)

          do I = 1,NTenerg
             if ((Energaux(I) > arctanhinit) .and.                      &
                  (Energaux(I) < arctanhfim)) then
                Energaux(I) = 2.0D0 * temp * ATANH(Energaux(I)) + EfLead
             elseif (Energaux(I) < arctanhinit) then
                Energaux(I) = 2.0D0 * temp * ATANH(arctanhinit) + EfLead
             elseif (Energaux(I) < arctanhfim) then 
                Energaux(I) = 2.0D0 * temp * ATANH(arctanhfim) + EfLead
             else
                print*, "there was some error"
             endif
          enddo
       endif

#ifdef MPI
       call MPI_Scatter (Energaux, NTenerg_div, MPI_Double_Precision,   &
                         Ei, NTenerg_div, MPI_Double_Precision,         &
                         0, MPI_Comm_world, MPIerror)
       call MPI_Scatter (gweightaux, NTenerg_div, MPI_Double_Precision, &
                         gweight, NTenerg_div, MPI_Double_Precision,    &
                         0, MPI_Comm_world, MPIerror)
       write (6,*) "Precision :", 0.5D0*SUM(gweight)
#else
       Ei = Energaux
       gweight = gweightaux
#endif

    ELSEIF (integraltype .eq. 'None') THEN

       do I = 1,NTenerg_div
          if (NTEnerg == 1) then
             Ei(I) = TEnergI 
          else
             Ei(I) = TEnergI + DBLE(I-1+Node*NTenerg_div)               &
                  / DBLE(NTenerg-1) * (TEnergF-TEnergI)
          endif
       enddo

       do I = 1,NTenerg_div
          If (MOD(NTenerg,2) /= 0) Then
             if (((Node == 0).and.(I == 1)) .OR.                        &
                  ((Node == Nodes-1).and.(I == NTenerg_div))) then
                gweight(I) = 1.0d0/3.0d0                                &
                     * (arctanhfim-arctanhinit) / DBLE(NTenerg-1)
             elseif (MOD(I+Node*NTenerg_div,2) == 0) then
                gweight(I) = 4.0d0/3.0d0                                &
                     * (arctanhfim-arctanhinit) / DBLE(NTenerg-1)
             else
                gweight(I) = 2.0d0/3.0d0                                &
                     * (arctanhfim-arctanhinit) / DBLE(NTenerg-1)
             endif
          Else

! Alex: Tomar cuidado, testar se NTenerg_div eh igual a 1.
             if ((Node == 0) .and. (I == 1)) then
                gweight(I) = 1.0d0/3.0d0                                &
                     * (TEnergF-TEnergI) / DBLE(NTenerg-1)
             elseif ((Node == Nodes-1) .and. (I == NTenerg_div-1)) then
                gweight(I) = (1.0d0/3.0d0+1.0d0/2.0d0)                  &
                     * (TEnergF-TEnergI) / DBLE(NTenerg-1)
             elseif ((Node == Nodes-1) .and. (I == NTenerg_div)) then
                gweight(I) = 1.0d0/2.0d0                                &
                     * (TEnergF-TEnergI) / DBLE(NTenerg-1)
             elseif (MOD(I+Node*NTenerg_div,2) == 0) then
                gweight(I) = 4.0d0/3.0d0                                &
                     * (TEnergF-TEnergI) / DBLE(NTenerg-1)
             else
                gweight(I) = 2.0d0/3.0d0                                &
                     * (TEnergF-TEnergI) / DBLE(NTenerg-1)
             endif
          EndIf
       enddo
    ELSE
#ifdef MPI
       call MPI_Abort (MPI_Comm_world, 1, MPIerror)
#else
       stop
#endif
    ENDIF ! IF (integraltype == '...


  end subroutine energygrid


!  *******************************************************************  !
!                               freegrid                                !
!  *******************************************************************  !
!  Description: free allocated vectors.                                 !
!                                                                       !
!  Written by Pedro Brandimarte, Oct 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2013                                    !
!  *******************************************************************  !
  subroutine freegrid


!   Free memory.
    deallocate (Ei, gweight)


  end subroutine freegrid


!  *******************************************************************  !


END MODULE idsrdr_engrid
