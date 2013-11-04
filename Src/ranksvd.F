C
C   subroutine ranksvd
C
      
      subroutine ranksvd (N, nspin, H1, S1, side_rank,
     &                    number_rank, VH, VS, Q)

!
!   Modules
!
      use parallel,        only: IOnode

      implicit none

#ifdef MPI
      include "mpif.h"
#endif

!     Input variables.
      integer, intent(in) :: N, nspin
      integer, dimension (nspin), intent(out) :: number_rank
      complex(8), dimension (N,N), intent(in) :: S1
      complex(8), dimension (N,N,nspin), intent(in) :: H1
      complex(8), dimension (N,N,nspin), intent(out) :: VH, VS, Q
      character (LEN=1), dimension (nspin), intent(out) :: side_rank

!     Local variables.
      integer :: ISPIN, I
      external :: rank, gensvd
#ifdef MPI
      integer :: MPIerror ! Return error code in MPI routines
#endif

      side_rank = " "

      if (IOnode) then
         do ISPIN = 1,nspin
            call rank (N, S1, H1(:,:,ISPIN), side_rank(ISPIN))
            If (side_rank(ISPIN) .NE. '0') Then

               call gensvd (N, S1, H1(:,:,ISPIN), side_rank(ISPIN),
     &                      number_rank(ISPIN), VH(:,:,ISPIN),
     &                      VS(:,:,ISPIN), Q(:,:,ISPIN))

            Else
               number_rank(ISPIN) = 0
               VH(:,:,ISPIN) = H1(:,:,ISPIN)
               VS(:,:,ISPIN) = S1(:,:)
               Q(:,:,ISPIN) = 0.0
               Do I = 1,N
                  Q(I,I,ISPIN) = 1.D0
               EndDo
            EndIf
         enddo
      endif

#ifdef MPI
      call MPI_Bcast (Q, N*N*nspin, MPI_Double_Complex, 0,
     &                MPI_Comm_world, MPIerror)
      call MPI_Bcast (VH, N*N*nspin, MPI_Double_Complex, 0,
     &                MPI_Comm_world, MPIerror)
      call MPI_Bcast (VS, N*N*nspin, MPI_Double_Complex, 0,
     &                MPI_Comm_world, MPIerror)
      call MPI_Bcast (number_rank, nspin, MPI_Integer, 0,
     &                MPI_Comm_world, MPIerror)
      call MPI_Bcast (side_rank, nspin, MPI_Character, 0,
     &                MPI_Comm_world, MPIerror)
#endif


      end subroutine ranksvd

