      Program OneDminimum                 !(P.Scharoch,R.Szymon 01.2020)
c
      implicit none
      real*8 xl,h,xmin
      real*8 fun
      integer n
      logical switch
c
      open(10,file='1Dminimum.inp',status='unknown')
      open(11,file='1Dminimum.out',status='unknown')
c
      read(10,*)xl,h
c
      xmin=xl
c
!-----------------------------------------------------------------------
c
      do n=1,10
         switch=.TRUE.
         do while (switch)
            if (fun(xmin).LE.fun(xmin+h))   then
               write(11,'(3F15.10)')h,xmin,fun(xmin)
               switch=.FALSE.
               xmin=xmin-h
            else
               write(11,'(3F15.10)')h,xmin,fun(xmin)
               xmin=xmin+h
            endif
         enddo
         write(11,*)
         h=h/10.d0
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
      real*8 FUNCTION fun(x)
c
      implicit none
      real*8 x

      fun=dexp(-1.d0/x**2)/x
c
      return
      end
c
!=======================================================================
