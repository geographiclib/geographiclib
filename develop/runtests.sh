#! /bin/sh
SUF=$1
shift
FLAGS=$*
TESTDIR=$HOME/geographiclib/scratch/triaxial
make Geod3Test
while read n l; do
    develop/Geod3Test -e $l $FLAGS < $TESTDIR/test$n.txt > test$n.out$SUF &
    sleep 1
done<<EOF
obl  1 3/4 3 0
seta 1 1   2 1
set  1 3/2 1 2
pro  1 3   0 3
spha 1 0 3 0
sphb 1 0 2 1
sphc 1 0 1 2
sphd 1 0 0 3
hu   1 7577279780/1130005142289 1942065235 6378137
phu  1 30308952520/4489860107041 3178376 971037955
oblx 1 150/197 49 1
prox 1 150/53 1 49
wgs84 1 124/18523 1 0
t321 1 2 3 5
orig 1 3801/10000 3600 201
bill 1 10101/10000 9900 201
EOF
exit

wgs84 e2 = 124/18523
    -> f = 124 / (sqrt(real(340804677)) + 18523)
         = 1/298.2572249059396272...

wgs84 1 595514447126/88957371407509.362414969 1 0

Problem case for phu:
echo 38 83 -38 -98 | ./Geod3Solve -e  1 30308952520/4489860107041 3178376 971037955 -p 20 -i --debug
3.1293641189285856421 -> 3.12936411894543332
err = 2e-11

set  err      bad-line inverse-problem                    alt
phu  151757   357892   38   83 -38  -98 -> -38 82  38 -97 -38 82 142 97
phu  105322   94146    51 -118  51 -119 -> -51 61 -51  62 -51 61 -51 62
pro  29083732 18182    89    1 -90   89 -> -90  1  89  89 -90  1  89 89
sphd 12622381 400668   41 -179 -42  104 -> -90  1  89  76 -90  1  89 76

This is fixed now with --hybridalt.

New problem case:

echo -89.991193792345 0.149841611536 179.999982614218144949 -3.13617456669843866158 |
    ./Geod3Solve $WGS84  -p 10
89.999972874420124 3.871490350067362 179.994355777621166 double
89.999998978116000 -0.000000000000028 179.850175774217501 mpfr

testwgs84a.txt
89.999998978116     0 179.850175774217528486

Simpler case:
echo -89.99119379 0 0.00001701 3.136174567039457 |
    ./Geod3Solve $WGS84  -p 10
89.999974331431957 -3.664060211355958 0.005835683217888 double
89.999998999999772 0.149810776922738 0.149793768059446  long double
89.999998999999988 0.149810810289633 0.149793800317347  quad
89.999998999999988 0.149810810289633 0.149793800317347  mpfr

echo -89.99119379 0 0.00001701 3.13617456 | ./Geod3Solve $WGS84  -p 10
89.999705776074919 -5.873495353757872 0.000509114381075 double
89.999998596669293 0.106758554270750 0.106741566607147  long double
89.999998596669794 0.106758614692179 0.106741604719894  quad
89.999998596669794 0.106758614692179 0.106741604719894  mpfr

s=3.136174506;for d in ../BUILD ../BUILD4; do
  make -C $d > /dev/null
  echo -89.99119379 0 0.00001701 $s | $d/Geod3Solve $WGS84  -p 10
done
-89.665897394122766 -22.822133956468143 0.000000448348734
89.999995502699376 0.033324465367030 0.033307455394745

Corresponding biaxial case:
echo -89.991223317891846397 .14984161153566253 179.9999826142150109 -20002951.0423183 |
    GeodSolve -p 10
89.999998981542007 -0.000000000939220 179.850175773281791 double
89.999998981542 -0.000000000000033 179.850175774220967    mpfr

GeodTest.dat
89.999998981542 0 179.850175774221
0.114 meters from pole

These ellipsoids form a series with a/c = 2
obl  1 3/4 3 0
seta 1 1   2 1
set  1 3/2 1 2
pro  1 3   0 3

These ellipsoids form a series with a/c = 1
spha 1 0 3 0
sphb 1 0 2 1
sphc 1 0 1 2
sphd 1 0 0 3

WGS84
A = 6378137;
f = 1/298.257223563;
B approx round(6378137*(1-1/298.257223563)) = 6356752

Earth approximation Panou et al., Hu et al. a-b = 70 = hu
-t A+35 A-35 B
-t 6378172 6378102 6356752
-e 6378102 0.00670552681260463 0.9967265474113045 0.0032734525886955056
-e 6378102 0.006706 0.996727 0.003273
-e 6378102 7577279780/1130005142289 1942065235 6378137

Prolate "equivalent" = phu
-t A B+35 B-35
-t 6378137 6356787 6356717
-e 6356787 0.006750533824532638 0.003262495093607705 0.9967375049063923
-e 6356787 0.006751 0.003262 0.996738
-e 6356787 30308952520/4489860107041 3178376 971037955

p:[e2=(a^2-c^2)/b^2, k2=(b^2-c^2)/(a^2-c^2), kp2 = (a^2-b^2)/(a^2-c^2)];
p,a=6378172,b=6378102,c=6356752;
p,a=6378137,b=6356787,c=6356717;

a^2/c^2=(1+e2*(1-k2))/(1-e2*k2);
e2=(a^2-c^2)/(a^2*k2+c^2*(1-k2));

obl  1 3/4 3 0
pro  1 3   0 3
spha 1 0 3 0
sphb 1 0 2 1
sphc 1 0 1 2
sphd 1 0 0 3

WGS84 1 595514447126/88957371407509.362414969 1 0
WGS84 1 595514447126/88957371407509.362414969 1 0

direct comparisons:

outb = baseline            hybridalt
testhu.outb 6.782325 1079
testobl.outb 7.776014 3524
testoblx.outb 6.114765 910
testphu.outb 7.491690 41788
testpro.outb 8.348299 2958
testprox.outb 8.718750 28263
testseta.outb 5.303099 896
testset.outb 5.660871 2924
testspha.outb 7.230974 2991
testsphb.outb 4.972728 575
testsphc.outb 5.046986 1419
testsphd.outb 7.132038 3343
testwgs84.outb 346710743.003688 585846147994897

outc = switch solve2 order
testhu.outc 6.825663 1595
testobl.outc 7.807088 3528
testoblx.outc 6.155997 918
testphu.outc 7.529627 41788
testpro.outc 8.393191 2963
testprox.outc 8.747921 28263
testseta.outc 5.305941 896
testset.outc 5.675739 1517
testspha.outc 7.247369 2992
testsphb.outc 4.982759 780
testsphc.outc 5.062428 1851
testsphd.outc 7.158010 3351
testwgs84.outc 12559947783.835535 17953999640092972

biaxial solution for direct
testobl.outd 7.896426 3520
testpro.outd 8.462096 2958
testspha.outd 7.345990 2992
testsphd.outd 7.206284 3343
testwgs84.outd 25990.371086 4115152238

gsolve + hybridalt -- new baseline
testhu.outf 6.815619 1364
testobl.outf 7.778114 3524
testoblx.outf 6.122106 563
testphu.outf 7.490330 41788
testpro.outf 8.365403 2958
testprox.outf 8.724349 28263
testseta.outf 5.295276 896
testset.outf 5.658956 2924
testspha.outf 7.225282 2992
testsphb.outf 4.972187 699
testsphc.outf 5.046515 988
testsphd.outf 7.128320 3343
testwgs84.outf 25990.273243 4115152240

tht fic.Ex fixes result in minor changes for oblx prox set seta
g is new baseline

make -j10 > /dev/null && head -42748 ../testobl.txt | tail -1 | ./Geod3Test $OBL

newt2d working for general case (not umbilical)
h is new baseline

improvements in newt2d (enforce montonicity of f and g)
i is new baseline
umbilical still needs testing

umbilical all bisection
head -3314 ../testset.txt | tail -1
70 0 180 0 0 180 0.46138515584816834054 0.42147199387998978727 0.7978133840499095064 0.68720816474403102059
echo 70 0 180 0.46138515584816834054 | ./Geod3Solve $SET

convergence failure
head -385 ../testprox.txt | tail -1
echo 63 -173 -61.97921997838416712 -4.64409746197940890408 | ./Geod3Solve $PROX

umblilical now working OK with newt2d
j is new baseline

meridional now using newt2d (renamed newt2)
k is new baseline

biaxial special fix
This is accurate but sometimes newt2 takes many iterations to converge.

To fix errors in testwgs84, restore oblpro logic removed between commits

  2025-05-08 3ae2c97add02569df6fd64b7372c06dd75930af3
  2025-05-09 3df35384b9f3a416cd31b34b9186be1c2ec10639

l is new baseline

Fix convergence issues with newt2 & testwgs84.

m is new baseline

mean ncoeffs
r thresh = 9/8 13.1 261 XX
s thresh = 7/8 13.2 261 XX
p thresh = 1/4 17.9 261
n thresh = 1/8 19.7 261
o thresh = 1/16 22.0 261
q thresh = -1/8 150.0 131073 XX

Peculiarity with line 2276 testset.txt inverse problem for -90 180 -16 152
compute with mpfr and -p 13 and -p 20 ->

-90 180 -31.678276203968430803 -16 152 -34.632417851313538856 0.54646403823611170442 0.50489789511517110353 0.80569806066626949844 0.73963052819047320946
-90 180 -31.6782762039684308026562168 -16 152 -34.6324178513135388559304218 0.546464038236111704419875801 0.50489789511517110353 0.80569806066626949844 0.73963052819047320946

Geod3Test $SET gives
double same for both
5 4 7 0 117 5 4 0 113
long double for both
2 1 2 0 97 4 5 0 97
quad
1 1 1 0 1 1 1 1 1719815
1 1 1 0 1 1 1 1 1
mpfr
1 1 1 0 1 1 1 1 2
1 1 1 0 1 1 1 1 1

The reverse direct problem
-16 152 -34.632417851313538856 -0.54646403823611170442
-16 152 -34.6324178513135388559304218 -0.546464038236111704419875801
double
-90.000000000000000 180.000000000000000 -31.678276203968121
-90.000000000000000 180.000000000000000 -31.678276203968121
long double
-90.000000000000000000 180.000000000000000000 -31.678276203968430977
-90.000000000000000000 180.000000000000000000 -31.678276203968430977
quad
-89.999999997721836852932 -179.999999999443832566959 -82.403869305235184494185
-89.999999999998168040374 -179.999999999999136199715 -97.052934618980198366583
mpfr
-89.999999997721836852932 -179.999999999443832566957 -82.403869305288181140469
-89.999999999998168040828 -179.999999999999136202544 -97.052854571677052532343
quad-mpfr
  0.000000000000000000000   -0.000000000000000000002   0.000000000052996646284
  0.000000000000000000454    0.000000000000000002829  -0.000080047303145834240

Large error=655 head -3280 ../testset.txt | tail -1 | ./Geod3Test $SET
Large CNT=147 head -1848 ../testset.txt | tail -1 | ./Geod3Test $SET

switch omg[12] in inverse calc

Redo baseline (minor differences)
                direct R      V        inverse     invdirect R     V
testbill.outa   4.6  309 15.9  3311    2.0   40    4.6  137  13.1   1941
testhu.outa     5.9  666  8.0  1133    2.7   84    5.6  233   8.1   1656
testobl.outa    3.3   27  5.3   246    1.4   12    3.2   25   5.1    243
testoblx.outa   4.8  310  7.7  1031    2.3   36    4.8  123   7.7    784
testorig.outa   4.9  212  6.6   618    2.4   39    5.1  221   6.9   1238
testphu.outa    5.9  672  9.0 27857    4.0 2084    7.0 2115   9.8   4479
testpro.outa    4.7   41  5.0    47    2.4   20    4.6   42   4.9     59
testprox.outa   7.2  318  9.9 19263    4.1  288    7.8 1157   8.9   1821
testset.outa    5.1  160  6.4  1492    2.8   88    6.0 8723   7.2   5132
testseta.outa   4.4  153  6.3  1275    2.4   82    5.0 2733   6.9   2197
testspha.outa   3.1   25  3.6    27    1.1   16    3.1   29   3.7     61
testsphb.outa   4.5  181  5.7   828    2.5   90    5.2 1498   6.4   1664
testsphc.outa   4.6  134  5.7  1419    2.6   88    5.4 2293   6.6   2431
testsphd.outa   2.8   20  3.2    20    1.1   12    2.8   23   3.3     54
testt321.outa   5.2  142  7.4   526    2.8  100    6.2 3985   8.2   1920
testwgs84.outa  4.5   32  5.1    32    2.0   20    4.3   31 126.8 210080

swapomg treatment reduces max inverse error for phu and prox
but max invdirect error becomes unacceptable (= 387838 for prox)

MATLAB results for testset
B     4.5   38  6.0    83   94.1 8424   111.8 8424 116.9   8428
longdoube results for testset
testset.outa    5.1  116  6.5  2490    2.8   86    5.9 4171 121.5 114274565
fix outlier by increase gamma=0 capture multiplier from 2 to 3
testset.outa3   5.1  116  6.6  2490    2.8   86    6.0 4171   7.4   2432

ODE direct only
                 direct R        V
testhu.outb      52.4    995   53.1    991
testobl.outb    121.0   2904  136.8   3731
testoblx.outb   118.0   2648  134.0   4233
testphu.outb     50.6    937   51.2    941
testpro.outb    122.7   2112  133.8   4191
testprox.outb   120.5   2794  129.8   3819
testseta.outb   120.0   3131  135.2   3822
testset.outb    121.1   3388  131.8   3692
testspha.outb    52.0   1198   52.6   1200
testsphb.outb    50.9    882   51.1    880
testsphc.outb    50.8    815   51.0    815
testsphd.outb    51.6    982   52.2    982
testwgs84.outb  102.4   1165  102.4   1165

testset.outb0   121.1   3388  131.8   3692
testset.outb0s with steps
  mean steps = 10.2 mean accel 513
testset.outb01   87.4   1756   92.6   1817 --normp
testset.outb1    35.9    936   39.3   1201 --dense
testset.outb1s with steps
  mean steps = 14.5 mean accel 1176
testset.outb11   41.2    423   45.9    867 --dense --normp
testset.outb11s with steps
  mean steps = 14.4 mean accel = 1118

ellipthresh = 1/8
diagnostics      ncoeff       2dctn     2dcntb     invcntn   invcntb
testhu.outd     20.18  99   4.91  25   0.30  22   8.12  19   1.26  11
testobl.outd    14.50  63   4.14   5   0.00   0   6.61  22   0.48  14
testoblx.outd   35.84 259   4.83  23   0.19  20   7.62  26   0.99  16
testphu.outd    20.67 106   4.93  37   0.31  22   9.57  26   2.32  14
testpro.outd    13.57  63   4.06   5   0.01   1   9.33  25   2.35  15
testprox.outd   22.72 125   5.20  37   0.31  21   8.74  23   1.54  14
testseta.outd   34.08 228   4.93  21   0.13  19   7.42  24   0.80  15
testset.outd    29.84 180   5.09  22   0.16  19   7.62  25   0.83  17
testspha.outd    3.51  63   1.94   2   0.00   0   6.85  18   0.60   9
testsphb.outd   26.21 106   4.91  21   0.15  19   7.84  19   1.00  10
testsphc.outd   26.26 106   4.91  21   0.15  19   8.09  20   1.10  12
testsphd.outd    3.50  63   1.94   2   0.00   0  10.08  26   3.00  17
testwgs84.outd   4.93  63   2.63   5   0.00   0  11.14  58   4.46  52

skip gamma = 0
testset.outc  29.84 180   5.09  22   0.16  19   7.62  25   0.83  17 all
testset.outc  29.94 180   5.04  16   0.11   8   7.73  25   0.84  17 gamma = 0
testset.outc  36.34 224   5.32  16   0.11   8   8.06  25   0.84  16 LD
testset.outc  64.83 414   6.13  17   0.11   8   9.06  27   0.85  16 Quad
testset.outc 148.73 971   7.31  18   0.11   8  10.50  32   0.85  18 MPFR

d=[53,64,113,256];
nc=[ 29.94, 36.34, 64.83,148.73];
nd=[5.04,5.32,6.13,7.31];
ni=[ 7.73, 8.06, 9.06,10.50];
p=polyfit(log(d),log(nc),1)
log(nc)=1.0175*log(d)-0.6392
log2(nc)=1.0175*log2(d)-0.9222
nd=1.4389*log(d)-0.6694
  =0.9973*log2(d)-0.6694
ni=1.7592*log(d)+0.7443
  =1.2194*log2(d)+0.7443
Increase nd by 1, d *= 2^(1/0.9973) = 2.0038
Increase ni by 1, d *= 2^(1/1.2194) = 1.7655

ellipthresh = 1/16 case e
testhu.oute     22.62  99   4.81  25   0.27  22   8.06  19   1.22  11
testobl.oute    15.41  75   4.14   5   0.00   0   6.68  22   0.51  14
testoblx.oute   37.14 259   4.86  23   0.19  20   7.56  26   0.95  17
testphu.oute    23.07 106   4.84  37   0.28  22   9.54  25   2.30  14
testpro.oute    14.48  75   4.06   5   0.01   1   9.42  25   2.40  15
testprox.oute   25.18 125   5.14  37   0.28  21   8.71  23   1.52  14
testseta.oute   37.00 228   4.95  21   0.12  19   7.34  24   0.76  15
testset.oute    33.49 180   5.08  22   0.14  19   7.56  25   0.79  17
testspha.oute    5.07  75   1.94   2   0.00   0   6.95  18   0.64   9
testsphb.oute   30.59 106   4.90  21   0.14  19   7.76  20   0.96  10
testsphc.oute   30.61 106   4.90  21   0.14  19   8.01  20   1.06  12
testsphd.oute    5.06  74   1.94   2   0.00   0  10.19  26   3.05  17
testwgs84.oute   5.68  75   2.63   5   0.00   0  11.19  58   4.48  52

ellipthresh = 1/4 case f
testhu.outf     18.45  99   4.99  25   0.32  22   8.18  19   1.29  11
testobl.outf    14.03  37   4.14   5   0.00   0   6.55  19   0.46  10
testoblx.outf   35.10 259   4.78  23   0.19  20   7.70  26   1.03  16
testphu.outf    18.96 106   5.01  37   0.32  22   9.61  26   2.34  14
testpro.outf    12.99  37   4.06   5   0.01   1   9.26  25   2.32  15
testprox.outf   20.83 125   5.27  37   0.34  21   8.79  23   1.56  14
testseta.outf   32.06 228   4.92  21   0.15  19   7.53  24   0.85  15
testset.outf    27.02 180   5.11  22   0.18  19   7.72  26   0.87  17
testspha.outf    2.32  37   1.94   2   0.00   0   6.75  18   0.56   9
testsphb.outf   22.25 106   4.93  21   0.16  19   7.96  19   1.06  10
testsphc.outf   22.28 106   4.93  21   0.16  19   8.20  20   1.15  12
testsphd.outf    2.31  37   1.94   2   0.00   0   9.97  26   2.95  17
testwgs84.outf   4.24  37   2.63   5   0.00   0  11.08  59   4.43  52
