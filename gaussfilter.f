**********************************************************************
*
*	Gaussian Filter Subroutine
*
*       includes an amplitude factor of sqrt(gwidth/pi) * 1/x0
*         from herrmann 1973 bssa pg 663-671 (eqn 13)
*
**********************************************************************

      subroutine gfilter(y_in,y_out,width,x0,dx,npts)
      real y_in(*),y_out(*), width,x0
      integer npts
      
      real x, wind, arg, pi
      integer k,j
      
      pi = acos(-1.0)
      w0 = 2 * pi * x0
      wc = w0 * sqrt(pi/width)
      dw = 2 * pi * dx
      
      scalefactor = 1.0 / (w0*w0)
      ampcor = pi * sqrt(width/pi) / (w0)
*     took away a factor of two from the paper because i keep
*        the negative frequencies here
*
                  
      do 1 k = 2, npts 
	j = k*2
	w = (k-1) * dw
	arg = -width * ( ((w - w0)*(w - w0)) * scalefactor )
	wind = ampcor * exp(arg)
	if( abs(w-w0) .gt. wc) wind = 0.0

        y_out(j-1) = y_in(j-1) * wind
        y_out(j)   = y_in(j)   * wind
1     continue
*           
*     the DC level and the nyquist frequency
*
      y_out(1) = 0.0
      y_out(2) = 0.0

      return
      end
**********************************************************************
*
*	Pre-whitening Gaussian Filter Subroutine
*
*       includes an amplitude factor of sqrt(gwidth/pi) * 1/x0
*         from herrmann 1973 bssa pg 663-671 (eqn 13)
*
**********************************************************************

      subroutine pw_gfilter(y_in,y_out,width,x0,dx,npts)
      real y_in(*),y_out(*), width,x0
      integer npts
      
      real x, wind, arg, pi
      integer k,j
      
      pi = acos(-1.0)
      w0 = 2 * pi * x0
      wc = w0 * sqrt(pi/width)
      dw = 2 * pi * dx
      
      scalefactor = 1.0 / (w0*w0)
      ampcor = pi * sqrt(width/pi) / (w0)
*     took away a factor of two from the paper because i keep
*        the negative frequencies here
*                  
      do 1 k = 2, npts 
	j = k*2
	w = (k-1) * dw
	arg = -width * ( ((w - w0)*(w - w0)) * scalefactor )
	wind = ampcor * exp(arg)
	if( abs(w-w0) .gt. wc) wind = 0.0
	
	whighten = sqrt(y_in(j-1)*y_in(j-1) + y_in(j)*y_in(j))
	if(whighten .lt. 1e-20) whighten = 1e20
	whighten = 1 / whighten

        y_out(j-1) = y_in(j-1) * wind * whighten
        y_out(j)   = y_in(j)   * wind * whighten
1     continue
*           
*     the DC level and the nyquist frequency
*
      y_out(1) = 0.0
      y_out(2) = 0.0

      return
      end
