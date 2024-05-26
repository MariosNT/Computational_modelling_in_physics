      Program CylindricalCapacitor        !(P.Scharoch,R.Szymon 01.2020)
c
      implicit none
c
      real*8 h,r(0:10000),fi(0:10000),finew(0:10000)
      real*8 omega,F,Fold,eps
      integer i,N
c
      open(10,file='capacitor.inp',status='unknown')
      open(11,file='capacitor1.out',status='unknown')
      open(12,file='capacitor2.out',status='unknown')
c
      read(10,*)N,r(0),r(N),omega,eps,fi(0),fi(N)
c
      h=(r(N)-r(0))/N
      Fold=1.0d6
c
      do i=1,N-1
         fi(i)=fi(0)+(fi(N)-fi(0))/(r(N)-r(0))*i*h
c     &         + (i-1)*(N-1-i)*0.001d0
      enddo
c
!-----------------------------------------------------------------------
c
100   do i=1,N-1
      finew(i)=0.5d0*(fi(i+1)+fi(i-1))+
     &         h*(fi(i+1)-fi(i-1))/4.d0/(r(0)+h*i)
      enddo
c
      do i=1,N-1
         fi(i)=finew(i)*omega+fi(i)*(1.d0-omega)
      enddo
c
      F=0
      do i=1,N
         F=F+(fi(i)-fi(i-1))**2*(r(0)+i*h-0.5d0*h)
      enddo
      F=F/2.d0/h
      if(dabs(F-Fold).gt.eps)   then
         Fold=F
         write(11,'(F15.7)')F
         goto 100
      endif
c
      do i=1,N
	 write(12,'(2F14.6)')r(0)+h*i,fi(i)
      enddo
c
      close(10)
      close(11)
      close(12)
c
      stop
      end
c
!=======================================================================
