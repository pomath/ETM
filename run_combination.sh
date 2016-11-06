#!/bin/bash

rm *.ps *.srt *.status *.sol *.warning *.prt *.org *.t logfile *.gdl psbase_ensum.* pshist_nrms.SUM.ensum pshist_wrms.SUM.ensum

ls ../glbf/*.glx > temp.gdl

list=`cat temp.gdl | awk '{print substr($1,9,11)}' | sort | uniq`

for item in $list
do

	grep ^"../glbf/$item" temp.gdl | awk '{print $1,"1.0"}' > weights.txt


	lines=`cat weights.txt | wc -l`
	cat weights.txt | awk -v nlines=$lines '{printf " %s %s ", $1, $2}; { if (NR != nlines) printf "+\n"; else printf "\n";}' >> glred.gdl
done

rm temp.gdl 

glred 6 combination.prt logfile glred.gdl globk_comb.cmd > glred.log
#sh_plot_pos -f globk_rep.org 

# save the important information from the log file
# columns: filename chi**2 weight
#echo "Creating previous run file..."
#paste <(grep "\<h.*Chi\>" glred.log | awk '{print $4,$9}') <(grep "\<Scaled by.*h.*.glx\>" glred.log | awk '{print substr($0,82,6)}') > prevrun.log

# extract the data from the prt file with awk
if [ ! -d "ts" ]
then
	mkdir ts
fi
echo "Creating time series file (for Matlab)..."
# grep ^pbo. combination.prt | grep -v "*" | awk '{print $2, $9, $27, $28, $29, $30, $31, $32, $10, $11, $12, $13, $14, $15}' > ts/time_series.txt
grep ^pbo. globk_rep.org | grep -v "*" | awk '{print $2, substr($0,50,10), substr($0,229,14), substr($0,243,14), substr($0,259,10), substr($0,269,8), substr($0,277,8), substr($0,285,8), substr($0,62,14), substr($0,77,14), substr($0,92,14), substr($0,106,8), substr($0,114,8), substr($0,122,8)}' > ts/time_series.txt

