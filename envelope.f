***************************************************************************      
*     Compute the envelope of the signal, s. 
*       The formulas used can be found in
*       Claerbout's book "Processing Versus Inversion"
*
*     Charles J. Ammon - Saint Louis University
*
*     April, 1998
***************************************************************************      
*
      subroutine envelope(s,n,dt)
      integer n
      real s(n),dt
      
      parameter(limit = 131104)
      complex spec(limit)
      integer stdout,forward,inverse,halfpts
      
      stdout = 6
      forward = 1
      inverse = -1
      twopi = acos(-1.0)
*
*     watch out for overflow
*      
      if(n .gt. limit) then
         write(stdout,*) ' '
         write(stdout,*)'ERROR: Too many points in envelope routine.'
         write(stdout,*)'decimate the file or increase array bounds.'
         write(stdout,*) ' '
	 stop
      end if
*
*     clean out the matrices
*      
      do 20 i = 1, limit
        spec(i) = cmplx(0.0,0.0)
20    continue
            
***************************************************************************      
*     Compute parameters for Fourier Transform
***************************************************************************      
      nfpts = 1
50    nfpts = nfpts * 2
      if(nfpts .lt. n) go to 50
*
*     watch out for overflow
*      
      if(nfpts .gt. limit) then
         write(stdout,*)' '
         write(stdout,*)'ERROR: Too many points in envelope routine.'
         write(stdout,*)'decimate the file or increase array bounds.'
         write(stdout,*)' '
	 stop
      end if
*      
      halfpts = nfpts / 2
      fnyquist = 1 / (2*dt)
      nyqpt = halfpts
      df = 1.0 / (nfpts * dt)
*
      do 60 i = 1, n
         spec(i) = cmplx(s(i),0.0)
60    continue
*
      call four1(spec,nfpts,forward)
*
***************************************************************************      
*    compute the FT of the envelope using the FT of the signal
*     zero the negative frequencies and scale the positive by 2
*
*      the factor of dt is from the forward transform
*
*    Also compute the FT of the derivative of the analytic signal
*      by scaling the spectrum by i*omega
***************************************************************************      
*      
      do 100 i = 1, halfpts
         w = (i-1) * df * twopi
         spec(i) = 2.0 * spec(i) * dt
	 spec(i+halfpts) = cmplx(0.0,0.0)
100   continue
*
      call four1(spec,nfpts,inverse)
*      
      do 150 i = 1,nfpts
         spec(i) = spec(i) * df
150   continue
*
***************************************************************************      
*     now we are back in the time domain and the analytic signal
*       is stored in the complex array spec
***************************************************************************      
*   
*
***************************************************************************      
*     COMPUTE THE ENVELOPE
***************************************************************************      
*   
      do 200 i = 1, n
         sreal =  real(spec(i))
	 simag = aimag(spec(i))
         s(i) = sqrt(sreal*sreal+simag*simag)
200   continue
*      
      return
      end

***************************************************************************      
***************************************************************************      
***************************************************************************      
***************************************************************************      
*
*     Compute the envelope plus the instantaneous frequency
*       of the signal (in hertz), s. The formulas used can be found in
*       Claerbout's book "Processing Versus Inversion"
*
*     This subroutine is a somewhat memory intensive
*
*     Charles J. Ammon - Saint Louis University
*
*     April, 1998
*
***************************************************************************      
      subroutine envelope_plus(s,n,dt,envelope,ifreq)
      integer n
      real s(n),dt,envelope(n),ifreq(n)
      
      parameter(limit = 131104)
      complex spec(limit),dgdt(limit),num(limit),denom(limit)
      integer stdout,forward,inverse,halfpts
      
      stdout = 6
      forward = 1
      inverse = -1
      twopi = acos(-1.0)
*
*     watch out for overflow
*      
      if(n .gt. limit) then
         write(stdout,*) ' '
         write(stdout,*)'ERROR: Too many points in envelope routine.'
         write(stdout,*)'decimate the file or increase array bounds.'
         write(stdout,*) ' '
	 stop
      end if
*
*     clean out the matrices
*      
      do 20 i = 1, limit
        spec(i)  = cmplx(0.0,0.0)
        dgdt(i)  = cmplx(0.0,0.0)
        num(i)   = cmplx(0.0,0.0)
        denom(i) = cmplx(0.0,0.0)
20    continue
            
***************************************************************************      
*     Compute parameters for Fourier Transform
***************************************************************************      
      nfpts = 1
50    nfpts = nfpts * 2
      if(nfpts .lt. n) go to 50
*
*     watch out for overflow
*      
      if(nfpts .gt. limit) then
         write(stdout,*)' '
         write(stdout,*)'ERROR: Too many points in envelope routine.'
         write(stdout,*)'decimate the file or increase array bounds.'
         write(stdout,*)' '
	 stop
      end if
*      
      halfpts = nfpts / 2
      fnyquist = 1 / (2*dt)
      nyqpt = halfpts
      df = 1.0 / (nfpts * dt)
*
      do 60 i = 1, n
         spec(i) = cmplx(s(i),0.0)
60    continue
*
***************************************************************************      
      call four1(spec,nfpts,forward)
*
*    compute the FT of the envelope using the FT of the signal
*     zero the negative frequencies and scale the positive by 2
*
*      the factor of dt is from the forward transform
*
*    Also compute the FT of the derivative of the analytic signal
*      by scaling the spectrum by i*omega
*      
***************************************************************************      
      do 100 i = 1, halfpts
         w = (i-1) * df * twopi
         spec(i) = 2.0 * spec(i) * dt
	 spec(i+halfpts) = cmplx(0.0,0.0)
	 dgdt(i) = spec(i) * cmplx(0.0,w)
	 dgdt(i+halfpts) = cmplx(0.0,0.0)
100   continue
*
      call four1(spec,nfpts,inverse)
      call four1(dgdt,nfpts,inverse)
*      
      do 150 i = 1,nfpts
         spec(i) = spec(i) * df
	 dgdt(i) = dgdt(i) * df
150   continue
*
***************************************************************************      
*     now we are back in the time domain and the analytic signal
*       is stored in the complex array spec, its derivative in dgdt
***************************************************************************      
*   
*
***************************************************************************      
*     COMPUTE THE ENVELOPE
***************************************************************************      
*   
      do 200 i = 1, n
         sreal =  real(spec(i))
	 simag = aimag(spec(i))
         envelope(i) = sqrt(sreal*sreal+simag*simag)
200   continue
*
*     call wsac1('the_envelope',envelope,n,beg,dt,nerr)
*
*
***************************************************************************      
*     COMPUTE THE INSTANTANEOUS FREQUENCY
***************************************************************************      
*
      do 300 i = 1, nfpts
        num(i)   = conjg(spec(i)) * dgdt(i)
	denom(i) = conjg(spec(i)) * spec(i)
300   continue
*
*     1-2-1 SMOOTHING of the ratio
*
      call csmooth121(num,nfpts)
      call csmooth121(denom,nfpts)
*
*     do the division
*     
      do 310 i = 1,n
        ifreq(i) = aimag( num(i) / denom(i) ) / twopi
310   continue
*      
      return
      end
*
*
***************************************************************************      
      subroutine smooth121(x,n)
***************************************************************************      
      real x(n)
      integer n
      
      s = 1.0 / 4.0
      x(1) = (x(2) + 2*x(1) + x(2)) * s
      do 10 i = 2, n - 1
        x(i) = (x(i-1) + 2*x(i) + x(i+1)) * s
10    continue

      x(n) = (x(n-1) + x(n-1) + 2*x(n)) * s
      
      return
      end
*
***************************************************************************      
      subroutine csmooth121(x,n)
***************************************************************************      
      complex x(n)
      integer n
      
      s = 1.0 / 4.0
      x(1) = (x(2) + 2*x(1) + x(2)) * s
      do 10 i = 2, n - 1
        x(i) = (x(i-1) + 2*x(i) + x(i+1)) * s
10    continue

      x(n) = (x(n-1) + x(n-1) + 2*x(n)) * s
      
      return
      end
***************************************************************************      
***************************************************************************      
***************************************************************************      
***************************************************************************      
c
c     Calculate the envelope of a function
c       based on a code to perform hilbert transforms
c
c     program to calculate the hilbert
c     transform of a ftn of arbitrary length
c     ( < max pts) in the time domain.
c
c     the transform is found by convolving the
c     trace with a 201 point fir filter
c     the filter impulse response is obtained
c     by windowing an ideal hilbert transformer
c     impulse response with a hamming window.
c  
c      x = input sac time series
c      h = FIR hilbert transform filter
c
***************************************************************************      

      subroutine td_envelope(x,npts)
      real x(npts)
      integer npts
      
      integer fpts,maxpts
      parameter (fpts = 401, maxpts = 5000)
      real hf(fpts), envx(maxpts)
      real pi
      
      do 10 i = 1, maxpts
         envx(i) = 0
10    continue

c     compute the filter coefficients
c
      pi = acos(-1.0)
      do 1 i=1,fpts
         k=i-201
         fk=float(k)
         if(k.eq.0)fk=0.00001
         a=0.54+0.46*cos(0.005*pi*fk)
c
c     fk*dx .eq. time
c
         b=1.0/(pi*fk)
         c=float((1-(-1)**k))
         hf(i)=a*b*c
1     continue
c
c     convolve the filter response and the input trace
c
      do 3 i=1,npts
       do 2 j=1,fpts
        kk=i+j-1
        envx(kk)=envx(kk)+x(i)*hf(j)
2     continue
3     continue
c
c     Compute the envelope = sqrt(x*x + y*y)
c     x = time series
c     envx = hilbert transfrom
c
      do 4 i=1,npts
         im = i+200
         envx(i) = x(i)*x(i)+envx(im)*envx(im)
         envx(i) = sqrt(envx(i))
4     continue

       do 5 i=1, npts
         x(i) = envx(i)
5      continue


      return
      end
***************************************************************************      
***************************************************************************      
***************************************************************************      
***************************************************************************      
      
      
