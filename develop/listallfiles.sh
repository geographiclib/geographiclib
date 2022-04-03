#! /bin/sh
cd $HOME
for d in geographiclib git/geographiclib-{c,fortran,java,js,octave,python}; do
    (
        cd $d
        git ls-tree -r HEAD --name-only |
            sed -e "s%^%$HOME/$d/%"
    )
done
