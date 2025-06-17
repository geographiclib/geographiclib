#! /bin/sh
while read name params; do
# params are b, e2, k2, kp2, bet, omg, alp
    echo $params | ./InvDump > $name.m
done <<EOF
oblmer  1 0.75 1 0 90 0 90
oblmerx 1 0.75 1 0 89.9 0 90
EOF
exit
oblmaj 1 0.75 1 0 0  0 90
oblgen 1 0.75 1 0 45 0 90
oblmer 1 0.75 1 0 90 0 90
promin 1 3    0 1 0 90 0
progen 1 3    0 1 0 45 0
promer 1 3    0 1 0  0 0
EOF
exit
trimaj    1 1.5 0.3333333 0.6666667 0          0         90
tricirca  1 1.5 0.3333333 0.6666667 42.70330   0         90
tricircb  1 1.5 0.3333333 0.6666667 87.52250   0         90
triumb    1 1.5 0.3333333 0.6666667 90         0        135
tritransa 1 1.5 0.3333333 0.6666667 90        10.15216  180
tritransb 1 1.5 0.3333333 0.6666667 90        39.25531  180
trimin    1 1.5 0.3333333 0.6666667 90        90        180
