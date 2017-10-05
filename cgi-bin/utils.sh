# Look up a key in a QUERY_STRING and return raw value
lookuprawkey () {
    QUERY="$1"
    KEY="$2"
    echo "$QUERY" | tr '&' '\n' | grep "^$KEY=" | tail -1 | cut -f2- -d=
}

# Decode raw value translating %XX, CR-LF to LF, and converting + and ,
# to space
decodevalue () {
    echo "$1" | sed \
        -e 's/\\/%5C/g' \
        -e 's/%0[dD]%0[aA]/%0A/g' \
        -e 's/%\([0-9a-fA-F][0-9a-fA-F]\)/\\x\1/g' -e s/%/%%/g |
    xargs -d '\n' printf | tr -s '+,\t' ' ' | sed -e 's/^ //' -e 's/ $//'
}

# Apply conversions for the various degree, minute, and second symbols
# %A0 -> nothing (non-breaking space)
# [%C2]%B0 [%C2]%BA %81%8B %E2%81%B0 %26%238304%3B (&#8304;) -> d
# %26%23730%3B (&#730;) -> d
# %91 %92 [%C2]%B4 %E2%80%B2 %26%238242%3B (&#8242;) %81%8C -> ' (%27)
# %93 %94 %E2%80%B3 %26%238243%3B (&#8243;) %81%8D -> " (%22)
# %26%238722%3B (&#8722;) -> -
# Remove left/right guillemot symbols %AB %BB used to quote examples.
# Then convert ' ' -> "
translate () {
    echo "$1" | sed \
        -e 's/%C2%A0//g' -e 's/%A0//g' \
        -e 's/%C2%[AB]B//g' -e 's/%[AB]B//g' \
        -e 's/%C2%B0/d/g' -e 's/%C2%BA/d/g' -e 's/%C2%B4/%27/g' \
        -e 's/%26%238304%3B/d/g' -e 's/%E2%81%B0/d/g' \
        -e 's/%26%23730%3B/d/g' \
        -e 's/%B0/d/g' -e 's/%BA/d/g' -e 's/%B4/%27/g' \
        -e 's/%9[12]/%27/g' -e 's/%9[34]/%22/g' \
        -e 's/%E2%80%B2/%27/g' -e 's/%E2%80%B3/%22/g' \
        -e 's/%26%238242%3B/%27/g' -e 's/%26%238243%3B/%22/g' \
        -e 's/%81%8B/d/g' -e 's/%81%8C/%27/g' -e 's/%81%8D/%22/g' \
        -e 's/%B4/%27/g' -e 's/%27%27/%22/g' -e 's/%26%238722%3B/-/g'
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
    test `echo "$VALUE" | tr -d '[ -~\n\t]' | wc -c` -ne 0 &&
    echo `date +"%F %T"` Unprintable "$RAWVAL" >> ../persistent/utilities.log
    echo "$VALUE"
}

# Look up ellipsoid parameter leaving only allowed characters (--/ -> -, ., /)
lookupellipsoid () {
    VALUE=`lookuprawkey "$1" "$2"`
    VALUE=`echo "$VALUE" | sed -e 's/%26%238722%3B/-/g'`
    VALUE=`decodevalue "$VALUE"`
    VALUE=`echo "$VALUE" | tr -cd '[0-9--/Ee]'`
    echo "$VALUE"
}

# Encode a string for inclusion into HTML.
encodevalue () {
    echo "$1" | sed -e 's/&/\&amp;/g' -e 's/"/\&quot;/g' \
        -e 's/>/\&gt;/g' -e 's/</\&lt;/g' -e "s/'/\&#39;/g" -e 's/`/\&#96;/g'
}

# Encode and convert d to &deg;
convertdeg () {
    encodevalue "$1" | sed -e 's/d/\&deg;/g'
}

# Generate GeoHack URL.  $1 $2 are real position; $3 $4 is displayed
# postion; $5 is link color
geohack () {
    echo "<a href=\"http://tools.wmflabs.org/geohack/geohack.php?params=$1;$2\" style=\"color:$5\">$(convertdeg "$3 $4")</a>"
}
