        subroutine getper(wn,NX,pmin,pmax,nper,permin,permax) 
c-----
c       automatically create periods from the [pmin,pmax] limits
c       wn()    R*4 array of periods
c       NX  I   dimension of array
c       pmin    R*4 minimum period
c       pmax    R*4 maximum period
c       nper    I*4 number of periods generated
c       permin  R*4 - minimum period to be used in trace
c       permax  R*4 - maximum period to be used in trace
c----
        integer MAXFRQ
        parameter (MAXFRQ=100)
        common /fval/ wo, noer, wnin, mper
            real*4 wo(MAXFRQ), wnin(MAXFRQ)
            integer*4 noer, mper
        integer NX, nper
        real wn(NX)
        integer key(MAXFRQ) 

        real pmin, pmax
        parameter (NP=40)
        real pfac(NP)
        data pfac/1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
     1      2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
     2      3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,
     3      5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5/

        nper = 0
c-----
c       determine starting power
c-----
        ymxlog = alog10(pmax)   
        ymmin  = alog10(pmin)   
        nocy = ymxlog - ymmin + 1 
        iy = ymmin
        if(ymmin .lt. 0)iy = iy - 1
        do 100 ii=iy,iy+nocy+2
            tenpow = 10.0**ii
            do 200 jj=1,NP
                p = pfac(jj)*tenpow
                if(p .ge. 0.99*pmin .and. p .le. 1.01*pmax
     1          .and. p.ge. 0.99*permin .and. p.le.1.01*permax)then
                    nper = nper + 1
                    mper = nper
                    wn(nper) = p
                    wnin(mper) = 1.0/p
                    key(mper) = mper
                    if(nper.eq.NX)go to 1000
                endif
 200        continue
 100    continue
 1000   continue
C       write(6,*)(wnin(i),i=1,mper)
C              call sort(wnin,key,mper)
C       write(6,*)(wnin(i),i=1,mper)
        return
        end

