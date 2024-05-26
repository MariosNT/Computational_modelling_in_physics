      Program PosionsEquation2D           !(P.Scharoch,R.Szymon 01.2020)
c
      implicit none
c
      integer i,j,Nx,Ny
      real*8 hx,x0,xN,hy,y0,yN
      real*8 F,Fold,omega,eps,pi
c
      logical vh
c
      real*8 rho(1:1000,1:1000),phi(1:1000,1:1000),phinew(1:1000,1:1000)
c
      real*8 frho,fphi,fnewphi,functional
c
      common /const/pi,hx,hy
c
      pi=3.14159d0
c
      open(10,file='2DPoisson.inp',status='unknown')
      open(11,file='2DPoisson1.out',status='unknown')
      open(12,file='2DPoisson2.out',status='unknown')
c
      read(10,*)Nx,Ny,x0,xN,y0,yN,omega,eps
c
      hx=(xN-x0)/(Nx-1)
      hy=(yN-y0)/(Ny-1)
      Fold=0.d0
      vh=.TRUE.
c
      do i=1,Nx
         do j=1,Ny
            rho(i,j)=frho(i,j)
            phi(i,j)=fphi(i,j,Nx,Ny)
         enddo
      enddo
c
!-----------------------------------------------------------------------
c
100   if (vh) then
         do i=2,Nx-1
            do j=2,Ny-1
               phinew(i,j)=fnewphi(rho(i,j),phi(i-1,j),phi(i+1,j),
     &                                      phi(i,j-1),phi(i,j+1))
            enddo
         enddo
      else
         do j=2,Ny-1
            do i=2,Nx-1
               phinew(i,j)=fnewphi(rho(i,j),phi(i-1,j),phi(i+1,j),
     &                                      phi(i,j-1),phi(i,j+1))
            enddo
         enddo
      endif
c
      do i=2,Nx-1
         do j=2,Ny-1
            phi(i,j)=phinew(i,j)*omega+phi(i,j)*(1.d0-omega)
         enddo
      enddo
c
      F=0.
c
      do i=2,Nx
         do j=2,Ny
            F= F + functional(rho(i,j),phi(i,j),phi(i-1,j),phi(i,j-1))
         enddo
      enddo
c
      write(12,'(F13.7)')F
      if(dabs(F-Fold).GT.eps) then
         Fold=F
         vh=.NOT.vh
         goto 100
      endif
c
      do i=1,Nx
         do j=1,Ny
           write(11,'(4F13.6)')x0+(i-1)*hx,y0+(j-1)*hy,phi(i,j),rho(i,j)
         enddo
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
      real*8 FUNCTION frho(i,j)
c
      implicit none
c
      integer i,j
c
      frho=0.
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION fphi(i,j,Nx,Ny)
c
      implicit none
c
      integer i,j,Nx,Ny
c
      if (i.EQ.1) then
         fphi=exp((-(j-50)**2)/200.d0)
      elseif (j.EQ.1) then
         fphi=exp((-(i-50)**2)/200.d0)
      elseif (i.EQ.Nx) then
         fphi=exp((-(j-50)**2)/200.d0)
      elseif (j.EQ.Ny) then
         fphi=exp((-(i-50)**2)/200.d0)
      else
         fphi=0.
      endif
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION fnewphi(rho,phi1,phi2,phi3,phi4)
c
      implicit none
c
      real*8 rho,phi1,phi2,phi3,phi4,pi,hx,hy
c
      common /const/pi,hx,hy
c
      fnewphi=(hx**2.d0*hy**2*4.d0*pi*rho+
     &             hy**2*(phi1+phi2)+
     &             hx**2*(phi3+phi4))/
     &            (2.d0*(hx**2+hy**2))
c
      return
      end
c
!============================functional=================================
c
      real*8 FUNCTION functional(rho,phi1,phi2,phi3)
c
      implicit none
c
      real*8 rho,phi1,phi2,phi3,pi,hx,hy
c
      common /const/pi,hx,hy
c
      functional=hx*hy*0.5d0*( (phi1-phi2)**2/hx**2 +
     &           (phi1-phi3)**2/hy**2 )-
     &            hx*hy*4.d0*pi*rho*phi1
c
      return
      end
c
!=======================================================================
