#!/bin/bash
for i in *.tre
do
    md5sum $i | awk '{print $1}' > $i.md5
done