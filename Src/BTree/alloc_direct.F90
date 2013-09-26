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
!                          MODULE alloc_direct                          !
!  *******************************************************************  !
!  Description: management of dynamic allocation at a chain file.       !
!  differently from the corresponding modules for internal memory, the  !
!  allocation procedure control here the file ("memory") and create a   !
!  new item by extending the file.                                      !
!                                                                       !
!  Based on P. Lignelet, "Structures de Donnees en Fortran 90/95",      !
!  Masson  (1996).                                                      !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !

MODULE alloc_direct

!
! Modules
!
  use fortran_except
  use fortran_init

  implicit none

  PRIVATE ! default is private

  integer, PUBLIC, parameter :: NULL = 0 ! similar to 'NULL' for pointers

  integer, PUBLIC, parameter :: DCB_UNKNOWN = 16 ! exception

  TYPE DCB; PRIVATE
     integer :: LUN ! file logical unit number
     integer :: EOF ! address of fictitious end of file
     integer :: HEAD ! addres of the top of the free list
     integer, pointer :: HEADER(:) ! user fields
     TYPE(T_INIT) :: INIT
  END TYPE DCB

  PUBLIC :: DCB, IOlength, LOAD_HEADER, ALLOC_DCB, FREE_DCB,            &
            SAVE_HEADER, IN_LIMITS, CLOSE_DCB

  integer, pointer :: LUN, EOF, HEAD, iHEADER(:) ! working variables
  TYPE(DCB), pointer :: p_LUN
  logical :: IS_OK


CONTAINS


!  *******************************************************************  !
!                              equivalence                              !
!  *******************************************************************  !
!  Description: procedure of dynamic redefinition. Also if the          !
!  abstract type DCB was correctly initialized. Returns 'IOSTAT = 0'    !
!  if succeed, otherwise returns 'DCB_UNKNOWN'.                         !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  character(*) file       : File system name                           !
!  integer unit            : Data file logical unit number              !
!  integer LRECL           : Size of items in the file                  !
!  ***************************** OUTPUT ******************************  !
!  integer IOSTAT          : Returns 0 if succeeds, else returns a      !
!                            Fortran error code                         !
!  *******************************************************************  !
  subroutine equivalence (IOSTAT)

!   Input variables.
    integer, intent (out), optional :: IOSTAT

    IS_OK = IS_INIT(p_LUN%INIT)

    if (IS_OK) then
       LUN => p_LUN%LUN
       EOF => p_LUN%EOF
       HEAD => p_LUN%HEAD
       iHEADER => p_LUN%HEADER
       if (present(IOSTAT)) IOSTAT = 0
    else
       call EXCEPTION_MESSAGE (p_LUN%INIT, IOSTAT, DCB_UNKNOWN,         &
            "DCB didn't initialize (by LOAD_HEADER)")
    endif


  end subroutine equivalence


!  *******************************************************************  !
!                               IOlength                                !
!  *******************************************************************  !
!  Description: returns the size of the file header (for the parameter  !
!  'RECL' of 'OPEN' subroutine).                                        !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  integer SIZE_HEADER     : Number of words of the user part from the  !
!                            header                                     !
!  ***************************** OUTPUT ******************************  !
!  integer IOlength        : Size of file header                        !
!  *******************************************************************  !
  integer function IOlength (SIZE_HEADER)

!   Input variables.
    integer, intent(in) :: SIZE_HEADER

!   Local variables.
    integer :: I

    inquire (iolength = IOlength) (I, I=-1,SIZE_HEADER)


  end function IOlength


!  *******************************************************************  !
!                              LOAD_HEADER                              !
!  *******************************************************************  !
!  Description: initializes the 'DCB' of file 'LUN' (previously open    !
!  for direct access). If the file is empty, creates the item 1         !
!  (header). In returns, 'HEADER' is created with the same of the       !
!  header stored in item 1 (or with an undetermided content if the      !
!  file was empty). In case of failure, the optional parameter          !
!  'IOSTAT' returns the Fortran error code.                             !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  integer LUN             : File logical unit number                   !
!  integer SIZE_HEADER     : Number of words of the user part from the  !
!                            header                                     !
!  ***************************** OUTPUT ******************************  !
!  TYPE(DCB) LUN_DCB       : DCB of file associated with 'LUN'          !
!  integer* HEADER(:)      : Header (head of data structure)            !
!  integer IOSTAT          : Returns 0 if succeeds, else returns a      !
!                            Fortran error code                         !
!  *******************************************************************  !
  subroutine LOAD_HEADER (LUN, SIZE_HEADER, LUN_DCB, HEADER, IOSTAT)

!   Input variables.
    integer, intent(in) :: LUN, SIZE_HEADER
    TYPE(DCB), intent(out), target :: LUN_DCB
    integer, pointer, dimension(:) :: HEADER
    integer, intent (out), optional :: IOSTAT

!   Local variables.
    character(11) :: acc, act, frm
    integer :: LRECL, iol, ios
    logical :: opn

    allocate (HEADER(SIZE_HEADER))
    LUN_DCB%HEADER => HEADER
    call INITIALIZE (LUN_DCB%INIT)
    p_LUN => LUN_DCB
    call equivalence
    LUN_DCB%LUN = LUN
    read (LUN, REC=1, IOSTAT=ios) iHEADER, EOF, HEAD

    if (ios /= 0) then ! error reading item 1 => creates it

       inquire (LUN, OPENED=opn, ACCESS=acc, ACTION=act,                &
                FORM=frm, RECL=LRECL)
       inquire (IOLENGTH=iol) iHEADER, EOF, HEAD
       HEAD = NULL
       if (opn .and. acc(:1)//act(:1)//frm(:1)=="DRU"                   &
            .and. iol<=LRECL) then
          EOF = 2
          call write_header
          ios = 0
       else
          print*, "File ", LUN, " does not allow reading item 1"
          EOF = HUGE(0)
       endif

    endif
    if (present(IOSTAT)) IOSTAT = ios


  end subroutine LOAD_HEADER


!  *******************************************************************  !
!                               ALLOC_DCB                               !
!  *******************************************************************  !
!  Description: procedure for allocation of an item, with the priority  !
!  from top of the free list. In case of failure, the previous free     !
!  list is ignored.                                                     !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  TYPE(DCB) LUN_DCB       : DCB of file associated with 'LUN'          !
!  ***************************** OUTPUT ******************************  !
!  integer iADDR           : Header (head of data structure)            !
!  integer IOSTAT          : Returns 0 if succeeds, else returns a      !
!                            Fortran error code                         !
!  *******************************************************************  !
  subroutine ALLOC_DCB (LUN_DCB, iADDR, IOSTAT)

!   Input variables.
    TYPE(DCB), intent(in), target :: LUN_DCB
    integer, intent(out) :: iADDR
    integer, intent (out), optional :: IOSTAT

!   Input variables.
    integer :: Koderr

    p_LUN => LUN_DCB
    call equivalence (IOSTAT) ! returns 'IS_OK'

    if (IS_OK) then
       DO
          if (HEAD == NULL) then ! free list
             iADDR = EOF
             EOF = EOF + 1
             Koderr = 0
          else ! levy on free list
             iADDR = HEAD
             read (LUN, REC=HEAD, IOSTAT=Koderr) HEAD
             if (Koderr /= 0) HEAD = NULL
          endif
          if (Koderr == 0) EXIT
       ENDDO
       call write_header ! bkp
    endif


  end subroutine ALLOC_DCB


!  *******************************************************************  !
!                               FREE_DCB                                !
!  *******************************************************************  !
!  Description: recovers item 1, with key 'iADDR' (stacks on top of     !
!  free list). In case of failure on accessing 'iADDR', it is ignored.  !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  TYPE(DCB) LUN_DCB       : DCB of file associated with 'LUN'          !
!  integer iADDR           : Header (head of data structure)            !
!  ***************************** OUTPUT ******************************  !
!  integer IOSTAT          : Returns 0 if succeeds, else returns a      !
!                            Fortran error code                         !
!  *******************************************************************  !
  subroutine FREE_DCB (LUN_DCB, iADDR, IOSTAT)

!   Input variables.
    TYPE(DCB), intent(in), target :: LUN_DCB
    integer, intent(in) :: iADDR
    integer, intent (out), optional :: IOSTAT

    p_LUN => LUN_DCB
    call equivalence (IOSTAT) ! returns 'IS_OK'

    if (IS_OK) then
       write (LUN, REC=iADDR, ERR=99) HEAD
       HEAD = iADDR
       call write_header ! bkp
    endif


99 end subroutine FREE_DCB


!  *******************************************************************  !
!                              SAVE_HEADER                              !
!  *******************************************************************  !
!  Description: "real time" recording of modifications on master-item.  !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  TYPE(DCB) LUN_DCB       : DCB of file associated with 'LUN'          !
!  ***************************** OUTPUT ******************************  !
!  integer IOSTAT          : Returns 0 if succeeds, else returns a      !
!                            Fortran error code                         !
!  *******************************************************************  !
  subroutine SAVE_HEADER (LUN_DCB, IOSTAT)

!   Input variables.
    TYPE(DCB), intent(in), target :: LUN_DCB
    integer, intent (out), optional :: IOSTAT

    p_LUN => LUN_DCB
    call equivalence (IOSTAT) ! returns 'IS_OK'

    if (IS_OK) call write_header


  end subroutine SAVE_HEADER


!  *******************************************************************  !
!                             write_header                              !
!  *******************************************************************  !
!  Description: "real time" recording of modifications on master-item.  !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !
  subroutine write_header

    write (LUN, REC=1) iHEADER, EOF, HEAD


  end subroutine write_header


!  *******************************************************************  !
!                               IN_LIMITS                               !
!  *******************************************************************  !
!  Description: check if the address 'iADDR' points to the region of    !
!  data of the file (it doesn't veirifies if 'iADDR' belongs or not to  !
!  the free list).                                                      !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  TYPE(DCB) LUN_DCB       : DCB of file associated with 'LUN'          !
!  integer iADDR           : Header (head of data structure)            !
!  ***************************** OUTPUT ******************************  !
!  logical IN_LIMITS       : Returns 0 if succeeds, else returns a      !
!                            Fortran error code                         !
!  *******************************************************************  !
  logical function IN_LIMITS (LUN_DCB, iADDR)

!   Input variables.
    TYPE(DCB), intent(in), target :: LUN_DCB
    integer, intent(in) :: iADDR

    IN_LIMITS = IS_INIT (LUN_DCB%INIT) .and. iADDR>1                    &
         .and. iADDR<LUN_DCB%EOF


  end function IN_LIMITS


!  *******************************************************************  !
!                               CLOSE_DCB                               !
!  *******************************************************************  !
!  Description: close the data structure.                               !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  ****************************** INPUT ******************************  !
!  TYPE(DCB) LUN_DCB       : DCB of file associated with 'LUN'          !
!  *******************************************************************  !
  subroutine CLOSE_DCB (LUN_DCB)

!   Input variables.
    TYPE(DCB), intent(in), target :: LUN_DCB

!   Local variables.
    logical :: F_open

    if (IS_INIT(LUN_DCB%INIT)) then
       inquire (LUN_DCB%LUN, OPENED=F_open)
       if (F_open) close (LUN_DCB%LUN)
    endif


  end subroutine CLOSE_DCB


!  *******************************************************************  !


END MODULE alloc_direct
