      program main
      implicit double precision (a-z)
      integer i,j,k
      dimension mass(301)
      dimension  lo20(301), lo40(301), lo80(301), lo160(301)
      dimension nlo20(301),nlo40(301),nlo80(301),nlo160(301)
      open (11,file="lo20.dat")
      open (12,file="lo40.dat")
      open (13,file="lo80.dat")
      open (14,file="lo160.dat")
      open (21,file="nlo20.dat")
      open (22,file="nlo40.dat")
      open (23,file="nlo80.dat")
      open (24,file="nlo160.dat")
      do i=1,301
         read (11,*) M,lo20(i)
         read (12,*) M,lo40(i)
         read (13,*) M,lo80(i)
         read (14,*) M,lo160(i)
         read (21,*) M,nlo20(i)
         read (22,*) M,nlo40(i)
         read (23,*) M,nlo80(i)
         read (24,*) mass(i),nlo160(i)
      enddo
      close(11)
      close(12)
      close(13)
      close(14)
      close(21)
      close(22)
      close(23)
      close(24)
      klo20 = lo20(301)/lo160(301)
      klo40 = lo40(301)/lo160(301)
      klo80 = lo80(301)/lo160(301)
      knlo20 = nlo20(301)/nlo160(301)
      knlo40 = nlo40(301)/nlo160(301)
      knlo80 = nlo80(301)/nlo160(301)
      do i = 1,301
         write (15,*) mass(i), lo20(i)/klo20, lo40(i)/klo40,
     -        lo80(i)/klo80, lo160(i)
         write (25,*) mass(i), nlo20(i)/knlo20, nlo40(i)/knlo40,
     -        nlo80(i)/knlo80, nlo160(i)
      enddo
      dlo20 = lo20(301)-lo160(301)
      dlo40 = lo40(301)-lo160(301)
      dlo80 = lo80(301)-lo160(301)
      dnlo20 = nlo20(301)-nlo160(301)
      dnlo40 = nlo40(301)-nlo160(301)
      dnlo80 = nlo80(301)-nlo160(301)
      do i = 1,301
         write (16,*) mass(i), lo20(i)-dlo20, lo40(i)-dlo40,
     -        lo80(i)-dlo80, lo160(i)
         write (26,*) mass(i), nlo20(i)-dnlo20, nlo40(i)-dnlo40,
     -        nlo80(i)-dnlo80, nlo160(i)
      enddo
      stop
      end
      
