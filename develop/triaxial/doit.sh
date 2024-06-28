#! /bin/sh
set -e
for ((q=4; q<=32; ++q)); do
    test -f errs-250-$q.txt && continue
    for ((k=q; k<=250; k+=32)); do
        octave-cli --eval "runtest(250,$k)" >& errs-250-$k.txt < /dev/null &
        sleep 1
    done
    wait
done
