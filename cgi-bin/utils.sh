# $Id$
lookupkey () {
    QUERY="$1"
    KEY="$2"
    # Convert degree (&B0) minute (&#8242;) second (&#8243;) symbols into d ' ".
    # and the windows variants(?) %81%8B, %81%8C, %81%8D
    echo "$QUERY" | tr '&' '\n' | grep "^$KEY=" | tail -1 | cut -f2- -d= |
    sed -e 's/%B0/d/g' -e 's/%26%238242%3B/%27/g' -e 's/%26%238243%3B/%22/g' \
	-e 's/%81%8B/d/g' -e 's/%81%8C/%27/g' -e 's/%81%8D/%22/g' \
	-e 's/\\/%5C/g' -e 's/%\([0-9a-fA-f][0-9a-fA-F]\)/\\x\1/g' -e s/%/%%/g |
    xargs -d '\n' printf | tr -s '+' ' '
}

encodevalue () {
    VALUE=$1
    echo "$VALUE" | sed -e 's/&/\&amp;/g' -e 's/"/\&quot;/g' \
	-e 's/>/\&gt;/g' -e 's/</\&lt;/g' -e "s/'/\&#39;/g" -e 's/`/\&#96;/g'
}

