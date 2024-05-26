      Program StationaryDiffusion         !(P.Scharoch,R.Szymon 01.2020)
c
      implicit none
c
      real*8 Am(0:10000),Ao(0:10000),Ap(0:10000)
      real*8 b(0:10000),fi(0:10000),h,D(0:10000),Dprim(0:10000)
      real*8 xN
      integer i,N
c
      open(10,file='diffusion.inp',status='unknown')
      open(11,file='diffusion.out',status='unknown')
c
      read(10,*)N,fi(0),fi(N),xN
      h=xN/N
c
      do i=0,N
         D(i)=0.1d0+10.d0*i*h
         Dprim(i)=10.d0
      enddo
c
!-----------------------------------------------------------------------
c
      do i=1,N-1
         Am(i)=D(i)-h*Dprim(i)/2.d0
         Ao(i)=-2.d0*D(i)
         Ap(i)=D(i)+h*Dprim(i)/2.d0
         b(i)=0.d0
      enddo
c
      call GaussEbS(Am,Ao,Ap,b,fi,N)
c
      do i=1,N
	 write(11,'(F15.6,D15.5)')i*h,fi(i)
      enddo
c
      close(10)
      close(11)
c
      stop
      end
c
!=======================================================================
c
      subroutine  GaussEbS(Am,Ao,Ap,b,fi,N)
c 
      implicit none
      real*8 Am(0:10000),Ao(0:10000),Ap(0:10000),b(0:10000),fi(0:10000)
      real*8 alfa(0:10000),beta(0:10000),gammai
      integer N,i
c
      beta(N-1)=fi(N)
      alfa(N-1)=0.0d0
c
      do i=N-1,1,-1
         gammai=-1.d0/(Ao(i)+Ap(i)*alfa(i))
	 alfa(i-1)=Am(i)*gammai
  	 beta(i-1)=(Ap(i)*beta(i)-b(i))*gammai
      enddo
c
      do i=0,N-1
         fi(i+1)=alfa(i)*fi(i)+beta(i)
      enddo
c
      return
      end
c
!=======================================================================
