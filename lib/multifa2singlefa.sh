#!/bin/bash
# $1: in file name (abs path)
# $2: out dir name (abs path)

# zcat -f: zungip (if needed) & cat
mkdir -p $2

if [ ${1##*.} = "gz" ]; then
  comm="zcat -f"
else
  comm="cat"
fi

eval $comm $1 | awk -v dir="$2" '{
        if (substr($0, 1, 1)==">") {filename=(dir "/" substr($0,2) ".fa")}
        print $0 > filename
}'
