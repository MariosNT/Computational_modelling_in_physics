      Program Hartree                     !(P.Scharoch,R.Szymon 01.2020)
c
      implicit none
c
      real*8 a,L,v0,h,hh,pi
      real*8 K1,Km,K2,dK,epsK
      real*8 kfun,Sfun,Efun1,Efun2,Efunm
      real*8 y1,ym1,yp1,a1,b1,c,Kzeroold
      integer Np,N,ip,i,j
      real *8 Kzero(1:20)
      real *8 x(1:100000),phi(1:20,1:100000),rho(1:100000),v(1:100000)
      real *8 vnew(1:100000),vzero(1:100000)
      real *8 Am(1:100000),Ap(1:100000),Ao(1:100000),B(1:100000)
c
      external kfun
c
      common /param/pi
c
      open(10,file='Hartree.inp',status='unknown')
      open(11,file='Hartree.out',status='unknown')
      open(12,file='functions.out',status='unknown')
      open(13,file='well.out',status='unknown')
c
      read(10,*)a,Np,v0,N
c
      pi=3.14159265d0
c
      dK=v0/100.d0
      epsK=v0*0.1d0**6
      K1=dK
      L=10.d0*a
c
      h=2.d0*L/Np
      hh=5.d0/12.d0*h**2
c
      vnew(1)=0.d0
      vnew(Np)=0.d0
c
      do i=1,Np
         x(i)=-L+i*h
         if (x(i).LT.-a/2.d0) then
            vzero(i)=N/2.d0+v0
            v(i)=N/2.d0+v0
         elseif (x(i).LE.a/2.d0) then
            vzero(i)=-N*x(i)**2/(2*a)
            v(i)=-N*x(i)**2/(2*a)
         else
            vzero(i)=N/2.d0+v0
            v(i)=N/2.d0+v0
         endif
         write(13,'(2F15.7)')x(i),v(i)
      enddo
c
      Kzeroold=1.d0
      Kzero(1)=2.d0
c
      do while (abs(Kzeroold-Kzero(1)).GE.epsK)
c
         Kzeroold=Kzero(1)
         K1=dK
         ip=0
c
         ym1=0.d0
         y1=exp(-L*sqrt(-kfun(K1,v(1))))
c
         do i=1,Np/2
            call Numerow(x(i+1),h,ym1,y1,yp1,kfun(K1,v(i+1)))
	    ym1=y1
 	    y1=yp1
         enddo
c
         Efun2=2.*(y1-ym1)/y1
c
100      K2=K1+dK
         Efun1=Efun2
c
         ym1=0.d0
         y1=exp(-L*sqrt(-kfun(K1,v(1))))
c
         do i=1,Np/2
            call Numerow(x(i+1),h,ym1,y1,yp1,kfun(K2,v(i+1)))
	    ym1=y1
 	    y1=yp1
         enddo
c
         Efun2=2.*(y1-ym1)/y1
c
         if(Efun1*Efun2.LE.0.d0)   then
            ip=ip+1
c
60          if((K2-K1).LT.epsK)   goto 260
c
            Km=(K1+K2)/2.
c
            ym1=0.
            y1=exp(-L*sqrt(-kfun(K1,v(1))))
c
            do i=1,Np/2
               call Numerow(x(i+1),h,ym1,y1,yp1,kfun(Km,v(i+1)))
	       ym1=y1
 	       y1=yp1
            enddo
c
            Efunm=2.*(y1-ym1)/y1
c
            if((Efun1*Efunm).LT.0.)then
	       K2=Km
	       goto 60
            else
	       K1=Km
	       goto 60
            endif
c
260         if (Km**2/2.d0.GT.v0)goto 360
            Kzero(ip)=Km
            write(11,'(i3,F12.6)')ip,Kzero(ip)**2/2.d0-v0
         endif
c
         K1=K2
         if(K2**2/2.d0.LE.v0)goto 100
c
         if (N.GT.(2*ip))goto 300
c
360      write(11,*)
c
         do i=1,int((N+1)/2.)
            phi(i,1)=0.
            phi(i,2)=exp(-L*sqrt(-kfun(Kzero(i),v(1))))
            phi(i,Np)=0.
            phi(i,Np-1)=(-1.d0)**(mod(i,2)+1)*phi(i,2)
            c=0.
            do j=2,Np/2
               call Numerow(x(j+1),h,phi(i,j-1),phi(i,j),phi(i,j+1),
     &                      kfun(Kzero(i),v(j+1)))
               phi(i,Np-j)=(-1.)**(mod(i,2)+1)*phi(i,j)
            enddo
            do j=2,Np/2
               c=c+2.d0*h*(phi(i,j+1)**2 +
     &         4.d0*phi(i,j)**2+phi(i,j-1)**2)/3.d0
            enddo
            do j=1, Np
               phi(i,j)=phi(i,j)/sqrt(c)
            enddo
         enddo
c
         do i=1, Np
            rho(i)=0.d0
         enddo
c
         do i=1,int((N+1)/2.d0)
            do j=1, Np
               rho(j)=rho(j)-2.d0*phi(i,j)**2
            enddo
         enddo
         if (mod(N,2).EQ.1) then
            do j=1, Np
               rho(j)=rho(j)-phi(int(N/2.d0)+1,j)**2
            enddo
         endif
c
         do i=2,Np
            Am(i)=1.d0
            Ap(i)=1.d0
            Ao(i)=-2.d0
            B(i)=hh/5.d0*(Sfun(x(i+1),rho(i+1)) +
     &      10.d0*Sfun(x(i),rho(i)) + Sfun(x(i-1),rho(i-1)))
         enddo
c
         call GaussEbS(Am,Ao,Ap,B,vnew,Np)
c
         a1=vnew(2)/h
         b1=vnew(Np-1)/h
c
         i=1
         do while (x(i).LE.-a/2.d0)
            i=i+1
            vnew(i)=vnew(i)-a1*(x(i)-x(1))
            vnew(Np-i+1)=vnew(Np-i+1)-b1*(x(Np)-x(Np-i+1))
         enddo
         a1=a1*(x(i)-x(1))
         do while (x(i).LT.0.)
            i=i+1
            vnew(i)=vnew(i)-a1
            vnew(Np-i+1)=vnew(Np-i+1)-a1
         enddo
c
         do i=1,Np
            v(i)=vzero(i)+vnew(i)
         enddo
c
         do i=1,Np
            write(12,'(6F15.7)')x(i),rho(i),vzero(i),vnew(i),v(i)
         enddo
         write(12,*)
      enddo
c
300   close(10)
      close(11)
      close(12)
c
      stop
      end
c
!=======================================================================
c
      SUBROUTINE Numerow(x,h,ym1,y1,yp1,K)
c
      implicit none
c
      real*8 x,h,y1,ym1,yp1
      real*8 hh,Ap1,Am1,A1,K
c
      hh=h**2/12.d0
      Ap1=1.d0+hh*K
      Am1=1.d0+hh*K
      A1= 2.d0*(1.d0-5.d0*hh*K)
c
      yp1=A1/Ap1*y1-Am1/Ap1*ym1
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION kfun(k,v)
c
      implicit none
c
      real*8 k,v
c
      kfun=k**2-2.*v
c
      return
      end
c
!=======================================================================
c
      subroutine  GaussEbS(Am,Ao,Ap,B,vnew,Np)
c
      implicit none
c
      real*8 Am(1:100000),Ao(1:100000),Ap(1:100000),B(1:100000)
      real*8 vnew(1:100000),alfa(1:100000),beta(1:100000),gamma
      integer Np,i
c
      beta(Np-1)=vnew(Np)
      alfa(Np-1)=0.
c
      do i=Np-1,2,-1
         gamma=-1.d0/(Ao(i)+Ap(i)*alfa(i))
         alfa(i-1)=Am(i)*gamma
         beta(i-1)=(Ap(i)*beta(i)-B(i))*gamma
      enddo
c
      do i=1,Np-1
         vnew(i+1)=alfa(i)*vnew(i)+beta(i)
      enddo
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION Sfun(x,rho)
c
      implicit none
c
      real*8 x,rho
      real*8 pi
c
      common /param/pi
c
      Sfun=4.d0*pi*rho
c
      return
      end
c
!=======================================================================
