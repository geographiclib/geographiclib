# $Id$

# Look up a key in a QUERY_STRING and return raw value
lookuprawkey () {
    QUERY="$1"
    KEY="$2"
    echo "$QUERY" | tr '&' '\n' | grep "^$KEY=" | tail -1 | cut -f2- -d=
}

# Decode raw value translating %XX and converting + and , to space
decodevalue () {
    echo "$1" | sed \
	-e 's/\\/%5C/g' -e 's/%\([0-9a-fA-f][0-9a-fA-F]\)/\\x\1/g' -e s/%/%%/g |
    xargs -d '\n' printf | tr -s '+,' ' '
}

# Convert degree (&B0) minute (&#8242;) second (&#8243;) symbols into d ' ".
# and also the Windows variants(?) %81%8B, %81%8C, %81%8D
translate () {
    echo "$1" | sed \
	-e 's/%B0/d/g' -e 's/%26%238242%3B/%27/g' -e 's/%26%238243%3B/%22/g' \
	-e 's/%81%8B/d/g' -e 's/%81%8C/%27/g' -e 's/%81%8D/%22/g'
}

# Look up and decode a key
lookupkey () {
    decodevalue `lookuprawkey "$1" "$2"`
}

# Look up, translate, and decode a key.  If result has unprintable
# characters, log the raw value.
lookupcheckkey () {
    RAWVAL=`lookuprawkey "$1" "$2"`
    VALUE=`translate "$RAWVAL"`
    VALUE=`decodevalue "$VALUE"`
    test `echo "$VALUE" | tr -d '[ -~\n]' | wc -c` -ne 0 &&
    echo `date +"%F %T"` Unprintable "$RAWVAL" >> ../persistent/utilities.log
    echo "$VALUE"
}    

# Encode a string for inclusion into HTML.
encodevalue () {
    echo "$1" | sed -e 's/&/\&amp;/g' -e 's/"/\&quot;/g' \
	-e 's/>/\&gt;/g' -e 's/</\&lt;/g' -e "s/'/\&#39;/g" -e 's/`/\&#96;/g'
}

