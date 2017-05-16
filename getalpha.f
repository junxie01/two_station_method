****************************************************************
*
*    variable Gaussian Width
*
****************************************************************

      subroutine getalpha(distance,period,alpha)
      real distance,period,alpha
      integer i
      parameter(npts = 6)
      real x(npts), a(npts)
      data x /25., 35., 45., 55., 60., 300./
      data a /0.025,0.011,0.009,0.008,0.007125,0.003/
            
      if(period .le. x(1)) then
        alpha = 0.5 * distance * a(1)
	return
      else if(period .gt. x(npts))then
        alpha = 0.5 * distance * a(npts)
	return
      end if
*      
      i = 1
10    i = i + 1
      if(period .gt. x(i)) go to 10
      
      dx = period - x(i-1)
      dadx = (a(i) - a(i-1)) / (x(i) - x(i-1))
      alpha = a(i-1) + dx * dadx
      alpha = alpha * 0.5 * distance
      
      return
      end
