#!/usr/bin/bash
# get numbers from *.csv
# turn them into 00XXXX
# copy folders with names 00XXXX to the destination


INPUT=sim.csv
IFS=','

var=`cat sim.csv`
for val in $var; do
	echo $val
	val1=`printf "%06g\n" $val`
	echo $val1
done

#val1 6 numbers



