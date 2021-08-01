#!/bin/bash
input="scf_problems.txt"
output="numbers_scf_problems.csv"
while IFS= read -r line
do
  #STRING=$line
  echo ${line:0:6}>>$output
done < "$input"

cat $output | tr '\n' ',' < $output > new_file.txt
