      function distance(stla,stlo,evla,evlo)
      real theta,dis
      real theta1,phi1
      real theta2,phi2
      real x1,y1,z1
      real x2,y2,z2
      real mult,m1,m2
      real stla,stlo,evla,evlo,distnce
      integer nargs
      character*40 args
      real pi,rad,R
      R=6371.0
      pi=3.14159265
      rad=pi/180.
      theta1=90-stla
      phi1=stlo+360
      if(stlo>=0)phi1=stlo
      
      theta2=90-evla
      phi2=evlo+360
      if(evlo>=0)phi2=evlo
      
      theta1=theta1*rad
      phi1=phi1*rad
      theta2=theta2*rad
      phi2=phi2*rad
      x1=R * sin(theta1) * cos(phi1)
      y1=R * sin(theta1) * sin(phi1)
      z1=R * cos(theta1)
      
      x2=R * sin(theta2) * cos(phi2)
      y2=R * sin(theta2) * sin(phi2)
      z2=R * cos(theta2)
      
      mult=x1 * x2 + y1 * y2 + z1 * z2
      dis = mult/(R*R)
      if (dis > 1) dis = 1
      theta=acos(dis)
      distance=theta * R
      end function
