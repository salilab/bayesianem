#!/bin/bash

for f in $@
do
	n=$(basename $f .rmf3)
	dir=$(dirname $f)
	if [[ -f $dir/../stat.$n.out ]]
	then
		echo "$dir/../stat.$n.out"
	else
		echo $f
	fi
done
