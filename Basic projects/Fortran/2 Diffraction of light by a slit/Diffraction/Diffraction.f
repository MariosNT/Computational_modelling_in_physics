      Program diffraction                 !(P.Scharoch,R.Szymon 01.2020)
c
      implicit none
c
      real*8 Ampl,a,d,yl,yr
      real*8 lambda,k,pi,hy,y,I
      real*8 ReCompA,ImCompA,simpson
      integer Ny,Nint,iy
c
      external ReCompA,ImCompA
c
      parameter (Nint=100, pi=3.141592654d0)
c
      common /dyfrpar/Ampl,k,d,y
c
      open(10,file='Diffraction.inp',status='unknown')
      open(11,file='Diffraction.out',status='unknown')
c
      read(10,*)Ampl,a,d,yl,yr,Ny
      lambda=1.0d0
c
      hy=(yr-yl)/(Ny-1)
      k=2.d0*pi/lambda
c
!-----------------------------------------------------------------------
c
      do iy=1,Ny
         y=yl+(iy-1)*hy
         I=simpson(-a/2.0d0,a/2.0d0,Nint,ReCompA)**2 +
     &	   simpson(-a/2.0d0,a/2.0d0,Nint,ImCompA)**2
c
         write(11,'(2F15.6)')y,I
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
      real*8 FUNCTION simpson(xl,xr,Nh,fun)
c
      implicit none
c
      real*8 x,h,fun,xl,xr
      integer ih, Nh
c
      h=(xr-xl)/Nh/2.d0
c
      simpson=0.d0
	x=xl+h
c
      do ih=1,Nh
         simpson=simpson + h*(fun(x+h)+4.d0*fun(x)+fun(x-h))/3.d0
	 x=x+2.d0*h
      enddo
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION ReCompA(x)
c
      implicit none
c
      real*8 x,k,r,d,y,Ampl
c
      common /dyfrpar/Ampl,k,d,y
c
      r=dsqrt(d**2+(y-x)**2)
      ReCompA=Ampl*dcos(k*r)/dsqrt(r)		
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION ImCompA(x)
c
      implicit none
c
      real*8 x,k,r,d,y,Ampl
c
      common /dyfrpar/Ampl,k,d,y
c
      r=dsqrt(d**2+(y-x)**2)
      ImCompA=Ampl*dsin(k*r)/dsqrt(r)		
c
      return
      end
c
!=======================================================================
