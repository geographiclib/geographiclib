#! /bin/sh
while true; do
    FILE=$1
    test -z "$FILE" && break
    DIRERR=`cat $FILE | cut -f6- -d' ' | mean | cut -f4,12 -d' '`
    INVERR=`cat $FILE | cut -f1 -d' ' | mean | cut -f4,12 -d' '`
    INVDIRERR=`cat $FILE | cut -f2-5 -d' ' | mean | cut -f4,12 -d' '`
    printf "%-14s %6.2f %6.0f   %6.2f %6.0f   %6.2f %6.0f\n" \
	   $FILE $DIRERR $INVERR $INVDIRERR
    shift
done
