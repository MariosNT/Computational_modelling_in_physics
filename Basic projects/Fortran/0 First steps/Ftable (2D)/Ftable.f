      Program FTABLE                      !(P.Scharoch,R.Szymon 01.2020)
c
      implicit none
      real*8 xl,xr,yl,yr,hx,hy,x,y
      real*8 myfunction
      integer Npx,Npy,ip1,ip2
c
      open(10,file='ftable.inp',status='unknown')
      open(11,file='ftable.out',status='unknown')
c
      read(10,*)xl,xr,Npx,yl,yr,Npy
c
      hx=(xr-xl)/(Npx-1)
      hy=(yr-yl)/(Npy-1)
c
!-----------------------------------------------------------------------
c
      do ip1=1,Npx
         x=xl+(ip1-1)*hx
	 do ip2=1,Npy
            y=yl+(ip2-1)*hy
	    write(11,'(3F15.9)')x,y,myfunction(x,y)
	 enddo
         write(11,*)
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
      real*8 FUNCTION myfunction(x, y)
c
      implicit none
      real*8 x
      real*8 y
c
      myfunction=dsin(x)*dsin(y)
c
      return
      end
c
!=======================================================================
