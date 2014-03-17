#! /bin/bash
AUX=`ls *.aux`
for file in $AUX
do
  bibtex $file
done
