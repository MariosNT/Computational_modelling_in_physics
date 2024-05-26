      Program QWell                       !(P.Scharoch,R.Szymon 01.2020)
c
      implicit none
      real*8 Vo,a,h,E
      real*8 Feven,Fodd,zero
      real*8 E1,E2,deltE,eps
      integer Np,ip
c
      external Feven,Fodd
c
      common /par/Vo,a
c
      open(10,file='qwell.inp',status='unknown')
      open(11,file='FevenFodd.out',status='unknown')
      open(12,file='Elevels.out',status='unknown')
c
      read(10,*)Vo,a,Np
c
      h=0.99d0*Vo/(Np-1)
c
!-----------------------------------------------------------------------
c
      do ip=1,Np
	 E=-Vo+0.001d0*Vo+(ip-1)*h
	 write(11,'(3F15.9)')E,Feven(E),Fodd(E)
      enddo
c
!-----------------------------------------------------------------------
c
      deltE=Vo/100.d0
      eps=Vo/1000.d0
      E1=-Vo+0.001d0*Vo
      ip=0
c
100   E2=E1+deltE
c
      if(Feven(E1)*Feven(E2).lt.0.d0)   then
         ip=ip+1
         write(12,'(i3,F12.6)')ip,ZERO(E1,E2,eps,Feven)
      endif
c
      if(Fodd(E1)*Fodd(E2).lt.0.d0)   then
         ip=ip+1
         write(12,'(i3,F12.6)')ip,ZERO(E1,E2,eps,Fodd)
      endif
c
      E1=E2
      if(E2.lt.0.99d0*Vo)goto 100
c
      close(10)
      close(11)
      close(12)
c
      stop
      end
c
!=======================================================================
c
      real*8 function ZERO(xll,xrr,eps,fun)
c
      implicit none
      real*8 xl,xr,xm,eps,xll,xrr
      real*8 fun
c
      xl=xll
      xr=xrr
c
100   xm=(xl+xr)/2.d0
      if((xr-xl).lt.eps) goto 200
c 
      if(fun(xl)*fun(xm).lt.0.0d0)then
	 xr=xm
	 goto 100 
      else
    	 xl=xm
	 goto 100
      endif
c
200   ZERO=xm
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION Feven(E)
c
      implicit none
c
      real*8 E,Vo,a,k
c
      common /par/Vo,a
c
      k=dsqrt(2.d0*(Vo+E))
      Feven=dsin(k*a/2.d0)-dcos(k*a/2.d0)*dsqrt(-2.d0*E)/k
c
      return
	end
c
!=======================================================================
c
      real*8 FUNCTION Fodd(E)
c
      implicit none
c
      real*8 E,Vo,a,k
c
      common /par/Vo,a
c
      k=dsqrt(2.d0*(Vo+E))
      Fodd=dsin(k*a/2.d0)+dcos(k*a/2.d0)*k/dsqrt(-2.d0*E)		!
c
      return
      end
c
!====================================================================
