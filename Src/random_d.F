C###################################################################
C      Subroutine random.F
C
C
C
C####################################################################
C author: Alexandre Reily Rocha
C e-mail: reilya@if.usp.br
C Universidade de Sao Paulo
C Departamento de Fisica de Materiais e Mecanica
C
C First Version:
C February 2007
C
C####################################################################
C
      subroutine random_d(Ndeffects,avgdist,dist,theta)
       
      use parallel,        only: IOnode
#ifndef IBM
      use ifport
#endif

      implicit none

#ifdef MPI
      include "mpif.h"
#endif

      integer :: Ndeffects,I,MPIError

      integer, dimension (3) :: now

      double precision :: avgdist,Pi,aux

      double precision, dimension (Ndeffects) :: dist,theta


#ifdef IBM
      integer ::  time
      external time
      integer, allocatable, dimension (:) :: seed
#endif


      Pi=4.0d0*ATAN(1.0d0)

      if (IOnode) then

#ifdef IBM
       now=0
       now(1)=time()
       call random_seed(SIZE=now(2))
       allocate(seed(now(2)))
       seed=now(1)+now(2)
       call random_seed(PUT=seed)
C       call random_seed
       call random_number(aux)
       print*, "aux1", aux
       deallocate(seed)
#else
       call itime(now)
       aux=drand(3600*now(1)+60*now(2)+now(3))
#endif
       DO I=1,Ndeffects

#ifdef IBM
        call random_number(aux)
        print*, "aux random",aux
        theta(i) = 2.0d0 * Pi * aux
#else
        theta(i) = 2.0d0 * Pi * drand(0)
#endif

       ENDDO


#ifdef IBM
        now(1)=time()
        call random_seed(GENERATOR=2)
        allocate(seed(now(2)))
        seed=now(1)+now(2)
        call random_seed(PUT=seed)
C        call random_seed
        call random_number(aux)
        print*, "aux", aux
        deallocate(seed)
#else
       call itime(now)
       aux=drand(3600*now(1)+60*now(2)+now(3))
#endif

       DO I=1,NDeffects       
#ifdef IBM
        call random_number(aux)
        dist(i) = 2.d0 * avgdist * aux
#else
        dist(i) = 2.d0 * avgdist * drand(0)
#endif
       ENDDO

      endif

#ifdef MPI
      call MPI_BCAST(theta,NDeffects,MPI_DOUBLE_PRECISION,0,
     & MPI_COMM_WORLD,MPIerror)
      call MPI_BCAST(dist,NDeffects,MPI_DOUBLE_PRECISION,0,
     & MPI_COMM_WORLD,MPIerror)
#endif

      end subroutine random_d

C *************************************************************

      integer function irandomizedeffects(ideffecttypes,
     .                                  unitweight,start)

#ifndef IBM
      use ifport
#endif

      implicit  none

      integer :: ideffecttypes,ii

      integer, dimension (3) :: now

      double precision :: aux,interval

      double precision, dimension (ideffecttypes) :: unitweight

      logical :: start



#ifdef IBM
      integer :: time
      external time
      integer, dimension (:), allocatable :: seed
#endif

      if (start) then

#ifdef IBM
       now(1)=time()
       call random_seed(SIZE=now(2))
       allocate(seed(now(2)))
       seed=now(1)+now(2)
       call random_seed(PUT=seed)
C       call random_seed
       call random_number(aux)
       deallocate(seed)
#else
       call itime(now)
       aux=drand(3600*now(1)+60*now(2)+now(3))
#endif
      else
#ifdef IBM
       call random_number(aux)
#else
       aux=drand(0)
#endif
       print*, "irandomizedeffects: aux 1",aux
      endif

      if (ideffecttypes.eq.0) then 
        irandomizedeffects=1
      else
        ii=1
        interval=unitweight(ii)
        do while ((aux .gt. interval).and.(ii.le.ideffecttypes))
         interval=interval+unitweight(ii+1)
         ii=ii+1
        enddo
         irandomizedeffects=ii  
      endif

      end function irandomizedeffects

C **********************************************************

      integer function irandomize_index(nunits,start)

#ifndef IBM
      use ifport
#endif

      implicit  none

      integer :: nunits

      integer, dimension (3) :: now

      double precision :: aux

      logical :: start


#ifdef IBM
      integer :: time
      external time
      integer, dimension (:), allocatable :: seed
#endif

      if (start) then

#ifdef IBM
       now(1)=time()
       call random_seed(SIZE=now(2))
       allocate(seed(now(2)))
       seed=now(1)+now(2)
       call random_seed(PUT=seed)
C       call random_seed
       call random_number(aux)
       deallocate(seed)
#else
       call itime(now)
       aux=rand(3600*now(1)+60*now(2)+now(3))
#endif
      else
#ifdef IBM
       call random_number(aux)
#else
       aux=drand(0)
#endif
       print*, "irandomize_index: aux 2",aux
      endif

      irandomize_index = NINT(aux*(nunits-1)+2)

      end function irandomize_index 


