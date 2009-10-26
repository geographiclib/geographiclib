# $Id$
lookupkey () {
    QUERY="$1"
    KEY="$2"
    echo "$QUERY" | tr '&' '\n' | grep "^$KEY=" | tail -1 | cut -f2- -d= |
    sed -e 's/\\/%5C/g' -e 's/%\([0-9a-fA-f][0-9a-fA-F]\)/\\x\1/g' -e s/%/%%/g |
    xargs -d '\n' printf | tr -s '+' ' '
}

encodevalue () {
    VALUE=$1
    echo "$VALUE" | sed -e 's/&/\&amp;/g' -e 's/"/\&quot;/g' \
	-e 's/>/\&gt;/g' -e 's/</\&lt;/g' -e "s/'/\&#39;/g" -e 's/`/\&#96;/g'
}

