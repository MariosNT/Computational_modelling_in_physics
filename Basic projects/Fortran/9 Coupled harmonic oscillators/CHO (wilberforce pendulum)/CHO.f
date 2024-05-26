      Program CoupledHarmonicOscillators  !(P.Scharoch,R.Szymon 01.2020)
c
      implicit none
c
      integer Np,n
      real*8 h,xE,k12,k1,k2
      real*8 x(0:1000000)
      real*8 y1(0:1000000),v1(0:1000000),y2(0:1000000),v2(0:1000000)
c
      common /parameters/k12,k1,k2
c
      external drfuny,drfunv
c
      open(10,file='CHO.inp',status='unknown')
      open(11,file='CHO1.out',status='unknown')
      open(12,file='CHO2.out',status='unknown')
c
      read(10,*)Np,x(0),xE,y1(0),v1(0),y2(0),v2(0),k12,k1,k2
c
      h=(xE-x(0))/Np
c
      do n=1,Np
         x(n)=x(0)+n*h
      enddo
c
!-----------------------------------------------------------------------
c
      do n=0,Np-1
         call Euler(x(n),y1(n),y1(n+1),h,drfuny,v1(n),0.d0,.TRUE.)
         call Euler(x(n),v1(n),v1(n+1),h,drfunv,y1(n),y2(n),.TRUE.)
         call Euler(x(n),y2(n),y2(n+1),h,drfuny,v2(n),0.d0,.FALSE.)
         call Euler(x(n),v2(n),v2(n+1),h,drfunv,y2(n),y1(n),.FALSE.)
      enddo
c
      do n=0, Np-1
         write(11,'(4F14.7)') x(n),dsin(y1(n)),dcos(y1(n)),y2(n)
         write(12,'(5F14.7)') x(n),.5*(k1*y1(n)**2+v1(n)**2),
     &                             .5*(k2*y2(n)**2+v2(n)**2),
     &                             .5*(k12*(y2(n)-y1(n))),
     &        .5*(k12*y2(n)*y1(n)+k1*y1(n)**2+k2*y2(n)**2+
     &                             v1(n)**2+v2(n)**2)
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
      SUBROUTINE Euler(x,y,yp1,h,fun,z1,z2,l)
c
      implicit none
c
      real*8 x,y,yp1,h,fun,z1,z2
      logical l
c
      yp1=y+h*fun(x,y,z1,z2,l)
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION drfuny(x,y,v1E,z2,l)
c
      implicit none
c
      real*8 x,y,v1E,z2
      logical l
c
      drfuny=v1E
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION drfunv(x,v,y1E,y2E,l)
c
      implicit none
c
      real*8 x,v,y2E,y1E
      logical l
c
      real*8 k12,k1,k2
      common /parameters/k12,k1,k2
c
      if (l) then
         drfunv=-k2*y1E-k12*y2E
      else
         drfunv=-k1*y1E-k12*y2E
      endif
c
      return
      end
c
!=======================================================================
