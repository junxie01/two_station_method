      program double_station
      include "/home/junxie/opt/sac/include/sacf.h"
c      include "/Users/junxie/opt/sac/include/sacf.h"
c      implicit none
      parameter (nptsmax=1000000,npmax=500)
      real t1,t2,sd,tb,te
      real tmp_sig(nptsmax),tt
      real stla1,stlo1,stla2,stlo2,vref(npmax)
      real pb,pe,pmin,pmax,dp,tmp_a,trust_amp,trust_ph
      real beg1,beg2,dt,width,period,df,pi,dist1,dist2,dist
      real nperiod(npmax),vb,ve,phv,v_ref,phv_pre,d_phv(npmax)
      real cc(nptsmax,npmax),ytmp(nptsmax),pref(npmax),phv_ref(npmax)
      real en_sig1(nptsmax),en_sig2(nptsmax),summ1,summ2,amax(npmax)
      real sig1(nptsmax),sig2(nptsmax),sig(nptsmax),dum(nptsmax)
      real evla,evlo,daz
      integer halfpts,npts
      integer n1,n2,nn,nno
      integer nwin,wlen,nfft
      integer nb,ne,np,nr,nz,nmax,np_ref
      integer nargs,nsamp1,nsamp2,nerr,i,j,nlen,kk
      character(256) error
      character(80) args,sacfile1,sacfile2,para,phvref
      character(80) names,name1,name2,output,out_put,output_disp
      out_put="test.sac"
      cc=0
      amax=0
      summ1=0
      summ2=0
      pi = acos(-1.0)
      nargs=iargc()
      if(nargs.ne.1)then
          write(*,*)'Usage:'
          write(*,*)'double_station para.dat'
          write(*,*)'para.dat:'
          write(*,*)'file1 file2 pmin pmax vmin vmax 
     &  out_amp out_dsp phv_reference'
          call exit(-1)
      endif
      call getarg(1,para)
      open(99,file=para) ! read in the parameter data
  300 read(99,*,end=200,err=200)sacfile1,sacfile2,pmin,pmax,vb,ve,output,output_disp,phvref
      open(9,file=phvref) ! read in the reference phase velocity
      i=1
  100 read(9,*,err=101,end=101)pref(i),phv_ref(i)
      i=i+1
      goto 100
  101 continue
      close(9)
      np_ref=i-1 ! number of reference phase velocity
!      write(*,*)'np_ref=',np_ref
! read in the first sac file which is the near one
      call rsac1(sacfile1,sig1,nsamp1,beg1,dt,nptsmax,nerr) 
      if(nerr.ne.0) then
         write(*,*)'Error reading in file: ',sacfile1
         call exit(-1)
      endif
      nlen=nsamp1
      write(*,*)"nlen=",nlen
      call getfhv('dist',dist1,nerr)
      call getfhv('stla',stla1,nerr)
      call getfhv('stlo',stlo1,nerr)
      call getfhv('evla',evla,nerr)
      call getfhv('evlo',evlo,nerr)
      call getfhv('b',beg1,nerr)
      do i=1,nsamp1
          summ1=summ1+sig1(i)**2
      enddo
      write(*,*)"sum1=",summ1
! read in the second sac file which is the further one
      call rsac1(sacfile2,sig2,nsamp2,beg2,dt,nptsmax,nerr)
      if(nerr.ne.0) then
         write(*,*)'Error reading in file: ',sacfile2
         call exit(-1)
      endif
      call getfhv('dist',dist2,nerr)
      call getfhv('stla',stla2,nerr)
      call getfhv('stlo',stlo2,nerr)
      do i=1,nsamp2
         summ2=summ2+sig2(i)**2
      enddo

      dist=distance(stla1,stlo1,stla2,stlo2)
      daz=dfaz(evla,evlo,stla1,stlo1,stla2,stlo2)
      !if ( daz .lt. 5)then ! the delta azimuth is smaller than 1
!      dist=1000
          write(*,'(2a30)')sacfile1,sacfile2
          write(*,'("dist="f10.4,"   dt="f7.2)')dist,dt
          npts=1
          do while(npts.lt.nsamp1)
             npts=npts*2
          enddo
          halfpts=npts/2
          df=1.0/(npts*dt)
          write(*,*)'npts=',npts
          write(*,*)'halfpts=',halfpts,'df=',df
      ! do fft transform for both z and r component
          if(dist1>dist2)then
            dum=sig2
            sig2=sig1
            sig1=dum
         endif
         call realft(sig1,halfpts,1)
         call realft(sig2,halfpts,1)
         sig1=sig1 * df
         sig2=sig2 * df
! the travel time respect to 6km/s
         tb=dist/ve 
! the travel time respect to 2km/s
         te=dist/vb
      ! shortest wave length should be one third of distance
      !pe=dist1/ve/3 
      !pb=0.01
      !dp=1.0/dt/2.0
      !dp=1.0/dp
      !pb=max(pb,dp)
      !call getper(nperiod,100,pmin,pmax,np,pb,pe) 
      !      write(*,*)'hello after get period '
      !      write(*,*)'ne=',ne,'nb=',nb
! number of periods
         np=int((pmax-pmin)/dt)+1
!      write(*,*)'np=',np
         do i=1,np
            nperiod(i)=pmin+(i-1)*(pmax-pmin)/(np-1)
         enddo
! the begin time of the cross-correlation
         beg1 = -dt* (nsamp1 - 1.0) + (beg2 - beg1) 
! number to begin the search respect to 6km/s
         nb=(-beg1+tb)/dt
         ne=(-beg1+te)/dt
!         write(*,*)'nb=',nb,"ne=",ne
         open(10,file=output)
! begin search at every peirod
         do i=1,np
            period=nperiod(i)
            do j=2,np_ref
               if(period<pref(j))then
! get the reference velocity
                  vref(i)=phv_ref(j-1)+(period-pref(j-1))*
     &              (phv_ref(j)-phv_ref(j))/(pref(j)-pref(j-1))
                  goto 103
               endif
            enddo
  103       continue
            if(j.eq.np_ref+1.or.vref(i).gt.10)then
               vref(i)=vref(i-1) 
            endif
!         write(*,'("verf("1i5")="f7.4)')i,vref(i)
!         write(*,*)'period=',period
      ! get the alpha parameter
      !   if(period.lt.10)then
      !      write(names,'(f3.1)')period
      !      names='0'//trim(names)
      !   else
      !      write(names,'(f4.1)')period
      !   endif 
      !   names=trim(names)//'.sac'
      !         write(names,*)period
      !   write(*,*)names
            period=1.0/period
! get alpha
            call getalpha(dist1,period,width)
!         width=1
      !   write(*,*)'width=',width
      ! do gaussian filter for both waveforms
            call gfilter(sig1,en_sig1,width,period,df,halfpts)
      !   call getalpha(dist2,period,width)
            call gfilter(sig2,en_sig2,width,period,df,halfpts)
            call realft(en_sig1,halfpts,-1)
            call realft(en_sig2,halfpts,-1)
      !   call wsac0(out_put, dum, en_sig1, nerr)
      !   call exit(-1)
      !   write(*,*)'hello'
      !   name1='gf-z-'//trim(adjustl(names))
      !   call wsac0(name1,dum,en_sig1,nerr)
      !   name2='gf-r-'//trim(adjustl(names))
      !   call wsac0(name2,dum,en_sig2,nerr)
            nwin = 1
            wlen = nlen
      !      wlen = 100
            nfft = 0
! do cross-correlation
         call crscor(en_sig1, en_sig2, nlen, nwin, wlen, SAC_RECTANGLE, 
     &               ytmp, nfft, error)
!         write(*,*)'hello again'
            tmp_sig=0 ! store the cross correlation function
            do j=1,nsamp1-1 
!         do j=1,nsamp1
               tmp_sig(j)=ytmp(nfft-nsamp1+j+1)
            enddo
            do j=1,nsamp2
               tmp_sig(nsamp1+j-1) = ytmp(j)
            enddo
            nfft = nsamp1 + nsamp2 - 1
            !tmp_sig=tmp_sig/sqrt(summ1*summ2)

      !call setnhv('npts',   nfft,    nerr)
      !call setfhv('delta',  dt,   nerr)
      !call setlhv('leven',  .true.,  nerr)
      !call setfhv('b',      beg1,    nerr)
      !call setfhv('o',      0,    nerr)
      !call wsac0(out_put, dum, tmp_sig, nerr)

!         call envelope(tmp_sig,nfft,dt)
      !   call exit(-1)
            nmax=1 ! find the maximam
      !   do j=nb,ne+1
            do j=1,ne-nb+1
               if (tmp_sig(j+nb-1)>amax(i))then
                   amax(i)=tmp_sig(j+nb-1)
                   nmax=j
               endif
      !      phv=dist/(j+beg1/dt) 
      !      write(10,*)nperiod(i),phv,tmp_sig(j)
               cc(j,i)=tmp_sig(j+nb-1) 
            enddo
            cc(:,i)=cc(:,i)/amax(i)
!         trust_amp=0.98*cc(kk,i)
            do j=1,ne-nb+1
      !      cc(j-nb+1,i)=cc(j-nb+1,i)/amax(i)
               kk=j
               tt=(kk+nb-1)*dt+beg1 ! which means the begin one is 6km/s
               phv=dist/tt
      !      phv_pre=dist/((kk-1+beg1)*dt) 
               phv_pre=dist/(tt+dt) 
               if((phv-vref(i))*(phv_pre-vref(i)).le.0)then
                  if(cc(kk+1,i).lt.cc(kk,i))then
                     tmp_a=(cc(kk,i)-cc(kk-1,i))*(cc(kk,i)-cc(kk+1,i))
                     do while(tmp_a.lt.0.and.kk.gt.2)
                        kk=kk-1
                       tmp_a=(cc(kk,i)-cc(kk-1,i))*(cc(kk,i)-cc(kk+1,i))
                     enddo
                  tt=(kk+nb-1)*dt+beg1 ! which means the begin one is 6km/s
                  amax(i)=dist/tt
                  goto 104
                  else
                  tmp_a=(cc(kk,i)-cc(kk-1,i))*(cc(kk,i)-cc(kk+1,i))
                     do while(tmp_a.lt.0.and.kk.lt.ne-nb+1)
                     kk=kk+1
                  tmp_a=(cc(kk,i)-cc(kk-1,i))*(cc(kk,i)-cc(kk+1,i))
                     enddo
                  tt=(kk+nb-1)*dt+beg1 ! which means the begin one is 6km/s
                  amax(i)=dist/tt
                  goto 104
                  endif
               endif
            enddo
  104       continue ! when the maximan founded.
            trust_amp=0.9*cc(kk,i)
            do k=kk,ne-nb+1 ! find the standerd error
               if (cc(k,i).le.trust_amp)then
                     tt=(k+nb-1)*dt+beg1 ! which means the begin one is 6km/s
!j                trust_v=dist/(tt-dt)+(trust_am-cc(k-1,i))*
!     &  (dist/tt-dist/(tt-dt))/(cc(k,i)-cc(k-1,i)) 
                     trust_v=dist/tt
                   goto 105
               endif 
            enddo
  105       continue
            d_phv(i)=abs(trust_v-amax(i))
            do j=1,ne-nb+1
               tt=(j+nb-1)*dt+beg1 ! which means the begin one is 6km/s
               phv=dist/tt
               write(10,*)nperiod(i),phv,cc(j,i)
            enddo
      !   write(*,*)"nmax=",nmax,"amax=",amax(i)
      !   amax(i)=dist/(nmax+beg1/dt)
         enddo
      
         close(10)
         open(11,file=output_disp)
         do i=1,np
            write(11,*)nperiod(i),amax(i),d_phv(i)
         enddo
         close(11)
      !else
      !  write(*,*)'The delta azimuth is larger than 1 degree'
      !kendif
      goto 300
  200 continue 
      close(99)
      end program
