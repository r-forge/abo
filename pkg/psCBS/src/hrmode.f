c     code to calculate the half-range mode proposed by David Bickel (2002)
c     fortran version of earlier R code; written by Venkat Seshan 4/30/2010
      subroutine hrmode(n, x, hrm)
      integer n
      double precision x(n), hrm

      integer ilo, ihi, n0, ihrng(2)
      double precision xr

      n0 = n
      ilo = 1
      ihi = n0
      xr = x(n) - x(1)

c     iteratively get the densest half range
      do 10 while ((n0 .gt. 2) .and. (xr .gt. 0))
         call hrsel(n0, x(ilo), ihrng)
         ihi = ihrng(2) + ilo - 1
         ilo = ihrng(1) + ilo - 1
         xr = x(ihi) - x(ilo)
         n0 = ihi - ilo + 1
 10   continue

      hrm = (x(ilo)+x(ihi))/2

      return
      end

      subroutine hrsel(n, x, ihrng)
      integer n, ihrng(2)
      double precision x(n)

      integer i, ilo, ihi, j, nov2, n0
      double precision xr, xw

      xr = (x(n) - x(1))/2
      nov2 = n/2

      if (x(nov2) - x(1) .lt. xr) then
c     start from left end if the left half is more dense
         i = 1
         j = nov2
c     first half range window from the smallest value
         do 10 while (x(j) - x(i) .le. xr)
            j = j+1
 10      enddo
         j = j-1
         ilo = i
         ihi = j
         n0 = ihi - ilo + 1
         xw = x(ihi) - x(ilo)
c     now shift the window to the right
         do 30 while (j .lt. n)
            j = j + 1
            do 20 while (x(j) - x(i) .gt. xr)
               i = i + 1
 20         enddo
            if (j-i+1 .gt. n0) then
               ilo = i
               ihi = j
               n0 = ihi - ilo + 1
               xw = x(ihi) - x(ilo)
            endif
            if ((j-i+1 .eq. n0) .and. (x(j) - x(i) .lt. xw)) then
               ilo = i
               ihi = j
               n0 = ihi - ilo + 1
               xw = x(ihi) - x(ilo)
            endif
 30      enddo
      else
c     start from right end if the right half is more dense
         i = nov2
         j = n
c     first half range window from the largest value
         do 60 while (x(j) - x(i) .le. xr)
            i = i-1
 60      enddo
         i = i+1
         ilo = i
         ihi = j
         n0 = ihi - ilo + 1
         xw = x(ihi) - x(ilo)
c     now shift the window to the right
         do 80 while (i .gt. 1)
            i = i - 1
            do 70 while (x(j) - x(i) .gt. xr)
               j = j - 1
 70         enddo
            if (j-i+1 .gt. n0) then
               ilo = i
               ihi = j
               n0 = ihi - ilo + 1
               xw = x(ihi) - x(ilo)
            endif
            if ((j-i+1 .eq. n0) .and. (x(j) - x(i) .lt. xw)) then
               ilo = i
               ihi = j
               n0 = ihi - ilo + 1
               xw = x(ihi) - x(ilo)
            endif
 80      enddo        
      endif
      ihrng(1) = ilo
      ihrng(2) = ihi

      return
      end
