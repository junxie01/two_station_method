      function dfaz(evla,evlo,stla1,stlo1,stla2,stlo2)
      real evla,evlo,stla1,stlo1,stla2,stlo2
      real pi,rad,az1,az2,dfaz
      real cosp,delta_o1,delta_o2,sinp,cosaz1,cosaz2
      real delta_o,cos_AB,cos_OB,cos_OA,sin_OB,sin_OA
      pi=3.14159265
      rad=pi/180.
      
      evla=(90-evla)*rad
      stla1=(90-stla1)*rad
      stla2=(90-stla2)*rad
      evlo=evlo*rad
      stlo1=stlo1*rad
      stlo2=stlo2*rad
      delta_o1=stlo1-evlo
      delta_o2=stlo2-evlo
      delta_o=stlo1-stlo2
      
      cos_OA=cos(evla)*cos(stla1)+sin(evla)*sin(stla1)*cos(delta_o1)
      sin_OA=sqrt(1-cos_OA**2)
      cos_OB=cos(evla)*cos(stla2)+sin(evla)*sin(stla2)*cos(delta_o2)
      sin_OB=sqrt(1-cos_OB**2)
      
      cos_AB=cos(stla1)*cos(stla2)+sin(stla1)*sin(stla2)*cos(delta_o)
      
      dfaz=(cos_AB-cos_OA*cos_OB)/(sin_OA*sin_OB)
      dfaz=acos(dfaz)
      end function
