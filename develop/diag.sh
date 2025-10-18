#! /bin/sh
# FOR NORMAL RUNS
while true; do
    FILE=$1
    test -z "$FILE" && break
    DIRERRR=`cat $FILE | cut -f6,8 -d' ' | tr ' ' '\n' |
      mean | cut -f4,12 -d' '`
    DIRERRV=`cat $FILE | cut -f7,9 -d' ' | tr ' ' '\n' |
      mean | cut -f4,12 -d' '`
    INVERR=`cat $FILE | cut -f1 -d' ' | mean | cut -f4,12 -d' '`
    INVDIRERRR=`cat $FILE | cut -f2,4 -d' ' | tr ' ' '\n' |
      mean | cut -f4,12 -d' '`
    INVDIRERRV=`cat $FILE | cut -f3,5 -d' ' | tr ' ' '\n' |
      mean | cut -f4,12 -d' '`
    printf "%-14s %4.1f %4.0f %4.1f %5.0f   %4.1f %4.0f   %4.1f %4.0f %5.1f %6.0f\n" \
           $FILE $DIRERRR $DIRERRV $INVERR $INVDIRERRR $INVDIRERRV
    shift
done
exit

# FOR DIAG RUNS (gamma != 0)
while true; do
    FILE=$1
    test -z "$FILE" && break
    printf "%-14s" $FILE
    grep DIAG $FILE | paste ~/git/triaxial/src/testset-gamma.txt - |
        grep '^.1' | cut -f2 | grep DIAG > /tmp/xx.txt
    FILE=/tmp/xx.txt
    NCOEF=`cat $FILE | cut -f2-5 -d' ' | tr ' ' '\n' |
        mean | cut -f4,12 -d' '`
    printf " %6.2f %3.0f" $NCOEF
    NCNTN=`cat $FILE | cut -f6 -d' ' | mean | cut -f4,12 -d' '`
    printf " %6.2f %3.0f" $NCNTN
    NCNTB=`cat $FILE | cut -f7 -d' ' | mean | cut -f4,12 -d' '`
    printf " %6.2f %3.0f" $NCNTB
    ICNTN=`cat $FILE | cut -f8 -d' ' | mean | cut -f4,12 -d' '`
    printf " %6.2f %3.0f" $ICNTN
    ICNTB=`cat $FILE | cut -f9 -d' ' | mean | cut -f4,12 -d' '`
    printf " %6.2f %3.0f" $ICNTB
    echo
    shift
done
exit

# FOR DIAG RUNS
while true; do
    FILE=$1
    test -z "$FILE" && break
    printf "%-14s" $FILE
    NCOEF=`grep DIAG $FILE | cut -f2-5 -d' ' | tr ' ' '\n' |
        mean | cut -f4,12 -d' '`
    printf " %6.2f %3.0f" $NCOEF
    NCNTN=`grep DIAG $FILE | cut -f6 -d' ' | mean | cut -f4,12 -d' '`
    printf " %6.2f %3.0f" $NCNTN
    NCNTB=`grep DIAG $FILE | cut -f7 -d' ' | mean | cut -f4,12 -d' '`
    printf " %6.2f %3.0f" $NCNTB
    ICNTN=`grep DIAG $FILE | cut -f8 -d' ' | mean | cut -f4,12 -d' '`
    printf " %6.2f %3.0f" $ICNTN
    ICNTB=`grep DIAG $FILE | cut -f9 -d' ' | mean | cut -f4,12 -d' '`
    printf " %6.2f %3.0f" $ICNTB
    echo
    shift
done
exit

# FOR ODE RUNS
while true; do
    FILE=$1
    test -z "$FILE" && break
    printf "%-14s" $FILE
    ERRR=`cat $FILE | cut -f1,3 -d' ' | tr ' ' '\n' | mean | cut -f4,12 -d' '`
    ERRV=`cat $FILE | cut -f2,4 -d' ' | tr ' ' '\n' | mean | cut -f4,12 -d' '`
    printf " %6.1f %6.0f %6.1f %6.0f" $ERRR $ERRV
    echo
    shift
done
exit
