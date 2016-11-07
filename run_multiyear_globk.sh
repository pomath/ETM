#!/bin/bash
year=$1
EXPT="scot"
rm -rf globk_run
mkdir globk_run
hfiles=''
for yr in $year
  do
  hfiles+="`find ${yr}/${yr}_??? -type f -name "h${EXPT}*"` "
done

for file in $hfiles
  do
  if [ ! -d globk_run/${file:(-3)} ]; then
    mkdir globk_run/${file:(-3)}
  fi
  cp $file globk_run/${file:(-3)}/
done

cd globk_run/

if [ "$?" == "1" ]; then
  echo "Error Setting Up Directories"
fi

mkdir gsoln
cp ../*comb.cmd gsoln/


start_day=${hfiles: 23:3}
end_day=${hfiles: -4:3}
echo ${year: 0:4} $start_day ${year: -4:4} $end_day
sh_glred -expt $EXPT -s ${year: 0:4} $start_day ${year: -4:4} $end_day -local -opt H G T >& sh_glred.log

if [ "$?" == "1" ]; then
  echo "Error Running sh_glred"
fi
gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=comb_.pdf gsoln/plot*/*.ps

if [ "$?" == "1" ]; then
  echo "Error Combining Files"
fi

##MAKE MATLAB FILE FOR ETM CODE
if [ ! -d "ts" ]
then
  mkdir ts
  touch ts/time_series.txt
fi
echo "Creating time series file (for Matlab)..."

array=`find . -name '*.org' | sort`

for day in ${array}
  do
grep ^pbo. ${array} | grep -v '*' | awk '{ print $2, $9, $10, $11, $12, $13, $14, $15, $27, $28, $29, $30, $31, $32}' >> ts/time_series.txt
done
awk '!seen[$0]++' ts/time_series.txt > ts/time2.txt
rm ts/time_series.txt
#okular comb_${year}.pdf &
