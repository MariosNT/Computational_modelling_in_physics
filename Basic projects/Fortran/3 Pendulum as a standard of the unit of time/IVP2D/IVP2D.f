      Program InitialValueProblem         !(P.Scharoch,R.Szymon 01.2020)
c
      implicit none
c
      real*8 h,y0,v0,x0,xE,k,x(0:1000000)
      real*8 yE(0:1000000),vE(0:1000000)
      integer Np, n, ih, Nh, multih
c
      common /dyfrpar/k
c
      external drfuny, drfunv
c
      open(10,file='IVP2D.inp',status='unknown')
      open(11,file='IVP2D1.out',status='unknown')
      open(12,file='IVP2D2.out',status='unknown')
c
      read(10,*)y0,v0,x0,xE,Np,multih,Nh,k
c
!-----------------------------------------------------------------------
c
      do ih=1, Nh
         h=(xE-x0)/Np
c
         do n=1,Np
	    x(n)=x0+n*h
         enddo
c
         yE(0)=y0
         vE(0)=v0
c
         do n=0,Np-1
            call Euler(x(n),yE(n),h,drfuny,yE(n+1),vE(n))
            call Euler(x(n),vE(n),h,drfunv,vE(n+1),yE(n))
            if (ih.eq.Nh)    then
               write(11,'(3F14.7)')x(n),yE(n),vE(n)
            endif
         enddo
c
         write(12,'(2F14.7)')-log10(h),yE(Np)**2/2.d0+vE(Np)**2/2.d0
         Np=Np*multih
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
c
      SUBROUTINE Euler(x,y,h,fun,yp1,z)
c
      implicit none
c
      real*8 x,y,h,fun,yp1,z
c
      yp1=y+h*fun(x,y,z)
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION drfuny(x,y,vE)
c
      implicit none
c
      real*8 x,y,vE
c
      drfuny=vE
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION drfunv(x,v,yE)
c
      implicit none
c
      real*8 x,v,yE,k
      common /dyfrpar/k
c
      drfunv=-k*yE
c
      return
      end
c
!=======================================================================
