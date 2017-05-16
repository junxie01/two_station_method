#!/bin/bash
tscale=`gmt gmtinfo -T0.1/2 a.out`
echo $scale
pb=`awk '{print $3}' para.dat`
pe=`awk '{print $4}' para.dat`
vb=`awk '{print $5}' para.dat`
ve=`awk '{print $6}' para.dat`
reg=-R50/300/2/6
reg=-R$pb/$pe/$vb/$ve
gmt makecpt -Cseis $tscale > a.cpt
gmt xyz2grd a.out -I1/0.1 $reg -Ga.grd
gmt grdsample a.grd -I0.1/0.01 -Gb.grd
gmt grdimage b.grd -Ca.cpt $reg -JX6i/4i -Ba10f10/a1f0.5WSne -P -Y2i -K>a.ps
gmt psxy maxap.txt -R -J -B -W1p,black -Ex -O -K >>a.ps
gmt psxy maxap.txt -R -J -B -Sa0.1i -Gblack -O -K >>a.ps
#awk '{print $3,$5}' COR_USBOZ_USSCIA.SAC_s_2_DISP.1 | gmt psxy -R -J -B -Sa0.1i -Gred -O -K >>a.ps
#gmt psxy predict.dat -R -J -B -Sa0.1i -Gred -O -K >>a.ps
gmt psxy `awk '{print $8}' para.dat` -R -J -B -Sa0.1i -Gred -O -K >>a.ps 
