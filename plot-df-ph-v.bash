#!/bin/bash
gmtset LABEL_FONT_SIZE 12p
ps=df_ph.ps
sac1=US.BOZ.00.BHZ.M.SACs
sac2=US.SCIA.00.BHZ.M.SACs
stlo1=`saclst stlo f $sac1 | awk '{print $2}'`
stla1=`saclst stla f $sac1 | awk '{print $2}'`
stlo2=`saclst stlo f $sac2 | awk '{print $2}'`
stla2=`saclst stla f $sac2 | awk '{print $2}'`
evlo=`saclst evlo f $sac1 | awk '{print $2}'`
evla=`saclst evla f $sac1 | awk '{print $2}'`
echo $stla1 $stlo1 $stla2 $stlo2
sta1=USBOZ
sta2=USSCIA

pscoast -R0/360/-90/90 -JG-100/39/2.7i -Ba30f15/a15f15 -W0.1p -A100 -P -K -Y8i>$ps
psxy -R -J -B -St0.2i -Gblack -O -K >>$ps<<eof
$stlo1 $stla1
$stlo2 $stla2
eof
psxy -R -J -B -m -W1p,black -O -K >>$ps<<eof
$evlo $evla
$stlo1 $stla1
>
$evlo $evla
$stlo2 $stla2
eof
pstext -R -J -B -O -K>>$ps<<eof
$stlo1 $stla1 10 0 4 BL $sta1
$stlo2 $stla2 10 0 4 BL $sta2
eof
tb1=`saclst dist f $sac1 | awk '{print $2/6}'`
tb2=`saclst dist f $sac2 | awk '{print $2/6}'`
te1=`saclst dist f $sac1 | awk '{print $2/3}'`
te2=`saclst dist f $sac2 | awk '{print $2/3}'`
tb=`echo $tb1 $tb2 | awk '{if($1<$2)print $1;else print $2}'`
te=`echo $te1 $te2 | awk '{if($1<$2)print $2;else print $1}'`
psbasemap -R$tb/$te/0/3 -JX3i/2i -Ba1000f500:"t(sec)":/f1Swne -K -P -O -X3i -Y0.3i>>$ps
pssac $sac2 $sac1 -R -JX -B -O -K -Ent-3 -M1i -W2p,black>>$ps
pstext -R -JX -B -O -K >>$ps<<eof
2000 1.2 10 0 4 ML $sta2
2000 2.2 10 0 4 ML $sta1
eof

sac3=sem.US.SCIA.LHZ.sem.sac
ccf=correlat-boz-scia.sac
sac<<eof
r  $sac1 $sac2
correlate
evaluate to evla &1,stla
evaluate to evlo &1,stlo
evaluate to stla &2,stla
evaluate to stlo &2,stlo
w 1.sac $ccf
cut 0 3000
r $ccf
ch evla %evla
ch evlo %evlo
ch stla %stla
ch stlo %stlo
w over
eof
echo 0 2.5 5.5 50 300 20 1 0.5 0.2 2 $ccf >param.dat
aftani param.dat phv_ref
echo 2 2.5 5.5 50 300 20 1 0.5 0.2 2 $sac3 >param.dat
aftani param.dat phv_ref

sac<<eof
r $sac3
interp delta 1 
w over
bp co 0.00333 0.02 n 4 p 2
ch o 0
w a.sac 
r $ccf
bp co 0.00333 0.02 n 4 p 2
ch o 0
w b.sac 
q
eof
tb=`saclst dist f $sac3 | awk '{print $2/6}'`
te=`saclst dist f $sac3 | awk '{print $2/3}'`
tb=0
te=1500
psbasemap -R$tb/$te/0/2 -JX3i/2i -Ba1000f500:"t(sec)":/f1Swne -K -P -O -X-3.3i -Y-2.7i >>$ps
pssac b.sac -R -JX -B -O -K -Ent-3 -M1i -W2p,black>>$ps
pssac a.sac -R -JX -B -O -K -Ent-3 -M1i -W2p,red>>$ps
pstext -R -JX -O -K -Gred>>$ps<<eof
1000 0.8 15 0 4 BL SEM
eof
pstext -R -JX -O -K -Gblack>>$ps<<eof
1000 1.4 15 0 4 BL CCF
eof

tscale=`minmax -T0.1/2 a.out`
echo $scale
reg=-R50/300/2/6
makecpt -Cseis $tscale > a.cpt
xyz2grd a.out -I1/0.1 $reg -Ga.grd
grdsample a.grd -I0.1/0.01 -Gb.grd
grdimage b.grd -Ca.cpt $reg -JX3i/2i -Ba50f25:"Period(sec)":/a1f0.5wSnE:"Phase velocity(km/sec)": -P -O -K -X3.3i >>$ps
psxy b.out -R -J -B -W1p,black -Ey -O -K >>$ps
psxy b.out -R -J -B -Sa0.1i -Gblack -O -K >>$ps
awk '{print $3,$5}' COR_USBOZ_USSCIA.SAC_s_2_DISP.1 | psxy -R -J -B -Sa0.1i -Gred -O -K >>$ps

psbasemap $reg -JX6.3i/3.8i -Ba50f25:"Period(sec)":/a1f0.5wSnE:"Phase velocity(km/sec)": -P -O -K -X-3.3i -Y-4.7i >>$ps
psxy sem_ph_v.txt -R -JX -B -O -K -W1p,red >>df_ph.ps
psxy eq_ph_v.txt -R -JX -B -O -K -W1p,black >>df_ph.ps
awk '{print $3,$5}' COR_USBOZ_USSCIA.SAC_s_2_DISP.1 | psxy -R -JX -B -O -K -W1p,black >>df_ph.ps
awk '{print $3,$5}' ${ccf}_2_DISP.1 | psxy -R -JX -B -O -K -W1pta/blue >>$ps
awk '{print $3,$5}' ${sac3}_2_DISP.1 | psxy -R -JX -B -O -K -W1pta/red >>$ps
psxy -R -JX -B -O -K -Sv -Gblack>>$ps<<eof
100 5 -60 1.5
200 5.5 -45 2 
250 4.2 90 1.7 
80 3.5 37 1.8
180 3.7 130 2
eof
pstext -R -JX -B -O -K -Gblack>>$ps<<eof
100 5 15 0 4 BR SF synthetic
200 5.5 15 0 4 BR CCFs 
250 4.2 15 0 4 TL two-station 
80 3.5 15 0 4 BR ccf_eq
180 3.7 15 0 4 TL sem_eq_ccf
eof
#awk '{print $3,$5}' 2.sac_2_DISP.1 | psxy -R -JX -B -O -K -W1p,green >>df_ph.ps
#awk '{print $3,$5}' 4.sac_2_DISP.1 | psxy -R -JX -B -O -K -W1p,purple >>df_ph.ps
#awk '{print $3,$5}' 2.sac_2_DISP.1 | psxy -R -JX -B -O -K -W1p,purple >>df_ph.ps
#psxy maxap.txt -R -JX -B -O -K -W1p,red >>df_ph.ps
#psxy predict.dat -R -JX -B -O -K -W1p,black >>df_ph.ps
