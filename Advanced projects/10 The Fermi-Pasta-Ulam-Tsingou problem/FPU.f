      Program FermiPastaUlam              !(P.Scharoch,R.Szymon 01.2020)
c
      implicit none
c
      integer Np,i,j,k,num,nom
      real*8 h,xi,xE,kk
      real*8 a,mass,alpha,beta
      real*8 norma,pro,pro2
      real*8 y(0:100,1:2),v(0:100,1:2)
      real*8 base(1:100,1:100)
      real*8 vector(1:100,1:100),vector2(1:100,1:100)
      real*8 energy(1:100),energyp
c
      external drfuny,drfunv
c
      common /const/h,a,mass,alpha,beta
c
      open(10,file='FPU.inp',status='unknown')
      open(11,file='FPU1.out',status='unknown')
      open(12,file='FPU2.out',status='unknown')
c
      read(10,*)xE,Np,num,a,mass,alpha,beta,nom
c
      kk=3.14159265d0/((num-1)*a)
      xi=0.d0
      h=(xE-xi)/Np
c
      do j=0,num-1
         y(j,1)=sin(1.d0*kk*j)
         v(j,1)=0.d0
      enddo
c
      do k=1,num-1
         norma=0.
         do j=0, num-1
            base(j,k)=sin(k*kk*j)
            norma=norma+(base(j,k))**2
         enddo
         do j=0, num-1
            base(j,k)=base(j,k)/sqrt(norma)
         enddo
      enddo
c
      do k=1,num-1
         vector(1,k)=0.d0
         vector2(1,k)=0.d0
         vector(num-1,k)=0.d0
         vector2(num-1,k)=0.d0
      enddo
c
      energyp=0.d0
      do k=1, num-1
      pro=0.d0
      do j=1, num-2
         pro=pro+base(j,k)*y(j,1)
      enddo
      do j=1, num-2
         vector(j,k)=pro*base(j,k)
      enddo
      do j=1, num-1
         energyp=energyp+
     &           alpha/2.*(vector(j,k)-vector(j-1,k))**2+
     &           beta/3.*(vector(j,k)-vector(j-1,k))**3
      enddo
      enddo
c
!-----------------------------------------------------------------------
c
      do i=0,Np-1
         if (mod(i,100000).eq.0) write (*,*)xi
         do j=1, num-2
            call RK5(xi,y(j,1),h,drfuny,y(j,2),v(j,1),0.d0,0.d0)
            call RK5(xi,v(j,1),h,drfunv,v(j,2),y(j,1),y(j+1,1),y(j-1,1))
         enddo
c
         if (mod(i,100000).eq.0) then
            energy(1)=0.d0
            do j=1, num-1
               energy(1)=energy(1)+
     &           alpha/2.d0*(y(j,1)-y(j-1,1))**2 +
     &           beta/3.d0*(y(j,1)-y(j-1,1))**3 +
     &           mass/2.d0*(v(j,1))**2
            enddo
c
            do j=0,num-1
               write(11,*)j,y(j,1),v(j,1)
            enddo
c
            do k=1,nom
               energy(k)=0.d0
               pro=0.d0
               pro2=0.d0
               do j=1, num-2
                  pro=pro+base(j,k)*y(j,1)
                  pro2=pro2+base(j,k)*v(j,1)
               enddo
               do j=1, num-2
                  vector(j,k)=pro*base(j,k)
                  vector2(j,k)=pro2*base(j,k)
               enddo
               do j=1, num-1
                  energy(k)=energy(k)+
     &                      alpha/2.d0*(vector(j,k)-vector(j-1,k))**2 +
     &                      beta/3.d0*(vector(j,k)-vector(j-1,k))**3 +
     &                      mass/2.d0*(vector2(j,k))**2
               enddo
            enddo
            write(12,'($f10.2)')xi
            do k=1,nom
               write(12,'($f12.8)')energy(k)/energyp
            enddo
            write(12,*)
         endif
         do j=1, num-2
            y(j,1)=y(j,2)
            v(j,1)=v(j,2)
         enddo
         xi=xi+h
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
      SUBROUTINE RK5(x,y,h,fun,yp1,y0,yjp1,yjm1)
c
      implicit none
c
      real*8 x,y,h,fun,yp1,y0,yjp1,yjm1
      real*8 k1,k2,k3,k4
c
      k1=h*fun(x,y,y0,yjp1,yjm1)
      k2=h*fun(x+h/2.d0,y+k1/2.d0,y0,yjp1,yjm1)
      k3=h*fun(x+h/2.d0,y+k2/2.d0,y0,yjp1,yjm1)
      k4=h*fun(x+h,y+k3,y0,yjp1,yjm1)
      yp1=y+(k1+2.d0*k2+2.d0*k3+k4)/6.d0
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION drfuny(x,y,v,o1,o2)
c
      implicit none
c
      real*8 x,y,v,o1,o2
c
      drfuny=v
c
      return
      end
c
      real*8 FUNCTION drfunv(x,y,y0,yp1,ym1)
c
      implicit none
c
      real*8 x,y,y0,yp1,ym1,h,a,mass,alpha,beta
c
      common /const/h,a,mass,alpha,beta
c
      drfunv=alpha/mass*(yp1+ym1-2.*y0)+
     &       beta/mass*((yp1-y0)**2-(y0-ym1)**2)
c
      return
      end
c
