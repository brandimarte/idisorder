!  *******************************************************************  !
!  Copyright (c) Smeagol Authors (2003-2005):  A. R. Rocha,             !
!                                              V. Garcia-Suarez,        !
!                                              S. Bailey,               !
!                                              C. J. Lambert,           !
!                                              J. Ferrer and            !
!                                              S. Sanvito               !
!                                                                       !
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS  !
!  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT    !
!  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS    !
!  FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE      !
!  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,  !
!  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES             !
!  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR   !
!  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)   !
!  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,  !
!  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)        !
!  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED  !
!  OF THE POSSIBILITY OF SUCH DAMAGE.                                   !
!                                                                       !
!  SMEAGOL IS DISTRIBUTED ONLY THROUGH THE OFICIAL WEBSITE              !
!  (www.smeagol.tcd.ie) UPON COMPLETION OF "SMEAGOL ACADEMIC LICENSE".  !
!                                                                       !
!  FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL smeagol@tcd.ie  !
!  *******************************************************************  !

!  *******************************************************************  !
!                                pasbias                                !
!  *******************************************************************  !
!  Description:                                                         !
!                                                                       !
!  Written by V. M. Garcia-Suarez, Oct 2003.                            !
!  Departamento de Fisica                                               !
!  Universidad de Oviedo                                                !
!  e-mail: victor@condmat.uniovi.es                                     !
!  ***************************** HISTORY *****************************  !
!  Original version:    October 2003                                    !
!  ****************************** INPUT ******************************  !
!  character*(*) str1     :                                             !
!  character*(*) str2     :                                             !
!  ***************************** OUTPUT ******************************  !
!  character*(*) pasbias  :                                             !
!  *******************************************************************  !

      character*(*) function pasbias (str1, str2)

!     Input variables.
      character*(*) str1, str2

!     Local variables.
      integer :: l, m, n

      m = len(str1)
      do 10 l = 1, m
         if (str1(l:l) .ne. ' ') then
            n = l
            goto 20
         endif
 10   continue

 20   pasbias = str1(n:m-1)//str2

      end

!  *******************************************************************  !
!                               pasbias2                                !
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
!  character*(*) str1     :                                             !
!  character*(*) str2     :                                             !
!  ***************************** OUTPUT ******************************  !
!  character*(*) pasbias2 :                                             !
!  *******************************************************************  !

      character*(*) function pasbias2 (str1, str2)

!     Input variables.
      character*(*) :: str1, str2

!     Local variables.
      integer :: l, m, n

      m = len(trim(str1))
      do l = 1,m
         if (str1(l:l) .ne. ' ') then
            n = l
            goto 20
         endif
      enddo

 20   pasbias2 = str1(n:m)//str2

      end

