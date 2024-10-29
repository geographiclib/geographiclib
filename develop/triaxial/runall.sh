#! /bin/sh

n=8
for d in 5 6 4 7 1 2 3; do
    for ((k = 1; k <= n; ++k)); do
        matlab-cli -batch "runtest($d,$n,$k)" > errs-$d-$k.txt& sleep 1
    done
    wait
    cat errs-$d-?.txt | sort > errs-matlab-$d.txt
done
