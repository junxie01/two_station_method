#!/bin/bash
nsac=`cat para.dat | wc -l`
for i in `seq $nsac` 
do
   amp=`awk -v n=$i 'NR==n{print $7}' para.dat`
   dsp=`awk -v n=$i 'NR==n{print $8}' para.dat`
   rdsp=`awk -v n=$i 'NR==n{print $9}' para.dat`
   ps=$dsp.ps
   if [ -e $dsp ];then
      n=`cat $dsp | wc -l `
      dat=`awk -v n=$n 'ARGIND==1{p1[NR]=$1;v1[NR]=$2;}
                        ARGIND==2{n2=NR-n;p2[n2]=$1;v2[n2]=$2;}
                        END{a=0;
                            for(i=1;i<=n;i++){
                                for(j=2;j<=n2;j++){
                                    if(p2[j-1]<p1[i]&&p2[j]>p1[i]){
                                       c1true=v2[j-1]+(v2[j]-v1[j-1])*(p1[i]-p2[j-1])/(p2[j]-p2[j-1]);
                                       a=(c1true-v1[i])*(c1true-v1[i])+a;
                                       break;
                                    }
                                }
                           }
                        printf "%d\n",a;
                        }' $dsp $rdsp`
      echo $dat
      if [ $dat -gt 30 ];then
         rm $dsp
      else
      pb=`awk -v n=$i 'NR==n{print $3}' para.dat`
      pe=`awk -v n=$i 'NR==n{print $4}' para.dat`
      vb=`awk -v n=$i 'NR==n{print $5}' para.dat`
      ve=`awk -v n=$i 'NR==n{print $6}' para.dat`
      reg=-R$pb/$pe/$vb/$ve
      psbasemap $reg  -JX6i/4i -Ba20f10:"Period(sec)":/a1f0.5WSne:"Phase velocity(km/sec)"::."$dat": -P -Y2i -K>$ps
      makecpt -Cseis $tscale > a.cpt
      if [ -e $amp ];then
         tscale=`minmax -T0.1/2 $amp`
         xyz2grd $amp -I1/0.1 $reg -Ga.grd
         grdsample a.grd -I0.1/0.01 -Gb.grd
         grdimage b.grd -Ca.cpt $reg -JX6i/4i -B -O -K >>$ps
      fi
      psxy $dsp -R -J -B -W1p,black -Ey -O -K >>$ps
      psxy $dsp -R -J -B -Sa0.1i -Gblack -O -K >>$ps
      psxy $rdsp -R -J -B -Sa0.1i -Gred -O >>$ps 
#awk '{print $3,$5}' $1 | psxy -R -J -B -Sd0.2i -Gred -O -K >>a.ps
#psxy predict.dat -R -J -B -Sa0.1i -Gred -O -K >>a.ps
      fi
   fi
done
