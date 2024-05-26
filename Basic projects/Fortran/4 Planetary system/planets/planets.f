      Program Planets                     !(P.Scharoch,R.Szymon 01.2020)
c
      implicit none
c
      real*8 m1,m2,dt,tE,dt2
      real*8 x1(3),x2(3),y1(3),y2(3),t
      real*8 v1x,v2x,v1y,v2y
      real*8 F1x,F2x,F1y,F2y,E,L
c
      integer Np, n
c
      common /masy/m1,m2
c
      open(10,file='planets.inp',status='unknown')
      open(11,file='planets.out',status='unknown')
c
      read(10,*)m1,m2,tE,Np
      read(10,*)x1(1),y1(1),x2(1),y2(1),v1x,v1y,v2x,v2y
c
      dt=(tE-0.d0)/Np
      dt2=dt*dt
c
      call F(x1(1),y1(1),x2(1),y2(1),F1x,F1y,F2x,F2y)
      call EL(x1(1),y1(1),x2(1),y2(1),v1x,v1y,v2x,v2y,E,L)
c
      x1(2)=x1(1)+v1x*dt+0.5*(F1x/m1)*dt2
      y1(2)=y1(1)+v1y*dt+0.5*(F1y/m1)*dt2
      x2(2)=x2(1)+v2x*dt+0.5*(F2x/m2)*dt2
      y2(2)=y2(1)+v2y*dt+0.5*(F2y/m2)*dt2
c
      write(11,'(11E13.5)')
     &  t,x1(1),y1(1),x2(1),y2(1),v1x,v1y,v2x,v2y,E,L
c
!-----------------------------------------------------------------------
c
      do n=1,Np-1
	 t=0.d0+n*dt
c
         call F(x1(2),y1(2),x2(2),y2(2),F1x,F1y,F2x,F2y)
c
         x1(3)=2.d0*x1(2)-x1(1)+dt2*F1x/m1
         y1(3)=2.d0*y1(2)-y1(1)+dt2*F1y/m1
c         x2(3)=2.d0*x2(2)-x2(1)+dt2*F2x/m2
c         y2(3)=2.d0*y2(2)-y2(1)+dt2*F2y/m2
c
         v1x=(x1(3)-x1(1))/2.d0/dt
         v1y=(y1(3)-y1(1))/2.d0/dt
c         v2x=(x2(3)-x2(1))/2.d0/dt
c         v2y=(y2(3)-y2(1))/2.d0/dt
c
         call EL(x1(2),y1(2),x2(2),y2(2),v1x,v1y,v2x,v2y,E,L)
c
         write(11,'(11E13.5)')t,x1(2),y1(2),x2(2),y2(2),
     &                        v1x,v1y,v2x,v2y,E,L
c
         x1(1)=x1(2)
         y1(1)=y1(2)
c         x2(1)=x2(2)
c         y2(1)=y2(2)
c
         x1(2)=x1(3)
         y1(2)=y1(3)
c         x2(2)=x2(3)
c         y2(2)=y2(3)
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
      subroutine F(x1,y1,x2,y2,F1x,F1y,F2x,F2y)
c
      real*8 x1,x2,y1,y2,F1x,F2x,F1y,F2y,m1,m2
      real*8 r1q,r2q,r12q
      real*8 M,G
c
      common /masy/m1,m2
c
      parameter (G=1.d0, M=100.d0)
c
      r1q=(x1*x1+y1*y1)**(1.5)
      r2q=(x2*x2+y2*y2)**(1.5)
      r12q=((x2-x1)**2+(y2-y1)**2)**(1.5)
c
      F1x=-G*M*m1*x1/r1q+G*m1*m2*(x2-x1)/r12q
      F2x=-G*M*m2*x2/r2q+G*m1*m2*(x1-x2)/r12q
      F1y=-G*M*m1*y1/r1q+G*m1*m2*(y2-y1)/r12q
      F2y=-G*M*m2*y2/r2q+G*m1*m2*(y1-y2)/r12q
c
      return
      end
c
!=======================================================================
c
      subroutine EL(x1,y1,x2,y2,v1x,v1y,v2x,v2y,E,L)
c
      real*8 m1,m2,x1,x2,y1,y2,v1x,v2x,v1y,v2y,E,L
      real*8 G,M
c
      common /masy/m1,m2
c
      parameter (G=1.d0, M=100.d0)
c
      r1=dsqrt(x1*x1+y1*y1)
      r2=dsqrt(x2*x2+y2*y2)
      r12=dsqrt((x2-x1)**2+(y2-y1)**2)
c
      v1s=v1x*v1x+v1y*v1y
      v2s=v2x*v2x+v2y*v2y
c
      E=m1*v1s/2.d0+m2*v2s/2.d0-G*m1*M/r1-G*m2*M/r2-G*m1*m2/r12
      L=(x1*v1y*m1-y1*v1x*m1)+(x2*v2y*m2-y2*v2x*m2)
c
      return
      end
c
!====================================================================
