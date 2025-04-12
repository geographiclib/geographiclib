#! /bin/sh
SUF=$1
while read n l; do
    ./Geod3Test -e $l < ../test$n.txt > test$n.out$SUF &
    sleep 1
done <<EOF
set 1 3/2 1 2
obl 1 3/4 3 0
pro 1 3   0 3
spha 1 0 3 0
sphb 1 0 2 1
sphc 1 0 1 2
sphd 1 0 0 3
EOF

exit
These ellipsoid form a series with a/c = 2
1 3/4 3 0
1 1   2 1
1 3/2 1 2
1 3   0 3
