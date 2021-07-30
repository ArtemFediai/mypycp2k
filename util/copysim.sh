#!/usr/bin/bash
# get numbers from *.csv
# turn them into 00XXXX (6 digits)
# copy folders with names 00XXXX to the destination


IFS=','
path_to_sim='/shared/user_data/New_GW_FromArtemToPatrick/data/sim'
path_to_csv='/home/ws/bh5670/work/b3lyp/100_random_numbers.csv'
path_to_sim_destination='sim'

var=`cat $path_to_csv`
for val in $var; do
	val1=`printf "%06g\n" $val`
	val2=$path_to_sim/$val1
	#echo $val2
	cp -r $val2 $path_to_sim_destination
        echo "I have copied $val2 to $path_to_sim_destination"
done

#val1 6 numbers



