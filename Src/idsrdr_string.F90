!  *******************************************************************  !
!  I-Disorder Fortran Code 2007-2014                                    !
!                                                                       !
!  Written by Alexandre Reily Rocha (reilya@ift.unesp.br),              !
!             Pedro Brandimarte (brandimarte@gmail.com) and             !
!             Alberto Torres (alberto.trj@gmail.com).                   !
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
!                         MODULE idsrdr_string                          !
!  *******************************************************************  !
!  Description: some tricks to deal with strings in fortran.            !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *******************************************************************  !

MODULE idsrdr_string


  implicit none
  
  PUBLIC :: STRconcat, STRpaste
  PRIVATE :: 


CONTAINS


!  *******************************************************************  !
!                               STRconcat                               !
!  *******************************************************************  !
!  Description: concatenates the strings 'str1' and 'str2' removing     !
!  the blanc spaces before and after the string 'str1'.                 !
!                                                                       !
!  Written by Pedro Brandimarte, Jul 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    July 2013                                       !
!  ****************************** INPUT ******************************  !
!  character*(*) str1                    : First string                 !
!  character*(*) str2                    : Second string                !
!  ***************************** OUTPUT ******************************  !
!  character*(*) str3                    : Concatenated string          !
!  *******************************************************************  !
  subroutine STRconcat (str1, str2, str3)

!   Input variables.
    character*(*), intent(in) :: str1, str2
    character*(*), intent(out) :: str3

!   Local variables.
    integer :: l, m

    m = len(trim(str1))
    do l = 1,m
       if (str1(l:l) /= ' ') exit
    enddo

    str3 = str1(l:m)//str2


  end subroutine STRconcat


!  *******************************************************************  !
!                               STRpaste                                !
!  *******************************************************************  !
!  Description: concatenates the strings 'str1' and 'str2' removing     !
!  the blanc spaces after the string 'str1' and before 'str2'.          !
!                                                                       !
!  Written by Pedro Brandimarte, Jul 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    July 2013                                       !
!  ****************************** INPUT ******************************  !
!  character*(*) str1                    : First string                 !
!  character*(*) str2                    : Second string                !
!  ***************************** OUTPUT ******************************  !
!  character*(*) str3                    : Concatenated string          !
!  *******************************************************************  !
  subroutine STRpaste (str1, str2, str3)

!   Input variables.
    character*(*), intent(in) :: str1, str2
    character*(*), intent(out) :: str3

!   Local variables.
    integer :: l, m

    m = len(trim(str1))
    do l = m,1,-1
       if (str1(l:l) /= ' ') exit
    enddo

    str3 = str1(1:l)//str2


  end subroutine STRpaste


!  *******************************************************************  !


END MODULE idsrdr_string
