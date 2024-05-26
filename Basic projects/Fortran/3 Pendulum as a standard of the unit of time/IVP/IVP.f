      Program InitialValueProblem         !(P.Scharoch,R.Szymon 01.2020)
c
      implicit none
c
      real*8 h,y0,x0,xE,x(0:1000)
      real*8 yE(0:1000),yAB3(0:1000),yRK5(0:1000)
      integer Np,n
c
      external drfun
c
      open(10,file='IVP.inp',status='unknown')
      open(11,file='IVP.out',status='unknown')
c
      read(10,*)y0,x0,xE,Np
c
      h=(xE-x0)/Np
c
      do n=1,Np
	 x(n)=x0+n*h
      enddo
c
!-----------------------------------------------------------------------
c
      yE(0)=y0
      do n=0,Np-1
         call Euler(x(n),yE(n),h,drfun,yE(n+1))
      enddo
c
!-----------------------------------------------------------------------
c
      yRK5(0)=y0
      do n=0,Np-1
         call RK5(x(n),yRK5(n),h,drfun,yRK5(n+1))
      enddo
c
!-----------------------------------------------------------------------
c
      yAB3(0)=y0
      yAB3(1)=yRK5(1)
      yAB3(2)=yRK5(2)
c
      do n=2,Np-1
         call AB3(x(n),yAB3(n-2),yAB3(n-1),yAB3(n),h,drfun,yAB3(n+1))
      enddo
c
!-----------------------------------------------------------------------
c
      do n=0,Np
         write(11,'(5F14.7)')x(n),yE(n),yAB3(n),yRK5(n),y0*dexp(-x(n))
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
      SUBROUTINE Euler(x,y,h,fun,yp1)
c
      implicit none
c
      real*8 x,y,h,fun,yp1
c
      yp1=y+h*fun(x,y)
c
      return
      end
c
!=======================================================================
c
      SUBROUTINE AB3(x,ym2,ym1,y,h,fun,yp1)
c
      implicit none
c
      real*8 x,ym1,ym2,y,h,fun,yp1
c
      yp1=y+
     &h/12.d0*(5.d0*fun(x-2.d0*h,ym2)-16.d0*fun(x-h,ym1)+23.d0*fun(x,y))
c
      return
      end
c
!=======================================================================
c
      SUBROUTINE RK5(x,y,h,fun,yp1)
c
      implicit none
c
      real*8 x,y,h,fun,yp1
      real*8 k1,k2,k3,k4
c
      k1=h*fun(x,y)
      k2=h*fun(x+h/2.d0,y+k1/2.d0)
      k3=h*fun(x+h/2.d0,y+k2/2.d0)
      k4=h*fun(x+h,y+k3)
      yp1=y+(k1+2.d0*k2+2.d0*k3+k4)/6.d0
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION drfun(x,y)
c
      implicit none
c
      real*8 x,y
c
      drfun=-y
c
      return
      end
c
!=======================================================================
