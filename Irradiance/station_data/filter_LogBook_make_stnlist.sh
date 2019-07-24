#!/bin/bash

# CODE TO FILTER CSV FILE DERIVED FROM GREENEDGE CRUISE XLS LOGBOOK
# USED TO MAKE STATION LIST CONTAINING COORDINATES, SAPLING TIME, STATION CODE,
# ICE COVER AND RELATIVE HUMIDITY (A PROXY FOR FOG)
# Martí Galí Tàpias 25 January 2017

# 0. Define input file name
fxls=LogBookScience_GreenEdge_leg1b.csv

# 1. Filter out rows that do not contain initial data for rosette cast
grep "cast" ${fxls} > log_rosette.tmp1

# ------------------------------------------------------------------------------------------
# Formatting tab delimited file with numeric values only and one header row for easy reading
# Has to be done in this order!
# ------------------------------------------------------------------------------------------ 

# 2. Replace string "Jun" by "06"
sed 's/Jun/06/g' log_rosette.tmp1 > log_rosette.tmp21

# 21. Replace string "Jul" by "07"
sed 's/Jul/07/g' log_rosette.tmp21 > log_rosette.tmp2

# 3. Replace string "2016 " by "2016,"
sed 's/2016 /2016,/g' log_rosette.tmp2 > log_rosette.tmp3

# 4. Replace string ",o." (where . replaces one character) by "," to separate degrees from minutes
sed 's/,o.//g' log_rosette.tmp3 > log_rosette.tmp4

# 5. Replace " " by ","
sed 's/ /,/g' log_rosette.tmp4 > log_rosette.tmp5

# 61. Replace "-06-" by ",6,"
sed 's/-06-/,6,/g' log_rosette.tmp5 > log_rosette.tmp61

# 6. Replace "-07-" by ",7,"
sed 's/-07-/,7,/g' log_rosette.tmp61 > log_rosette.tmp6

# 7. Replace "o" by ","to separate degrees from minutes
sed 's/o/,/g' log_rosette.tmp6 > log_rosette.tmp7

# 8. Replace ":" by ","
sed 's/:/,/g' log_rosette.tmp7 > log_rosette.tmp8

# 9. Replace "9+" by "9.5" (last field, ice cover)
sed 's/9+/9.5/g' log_rosette.tmp8 > log_rosette.tmp9

# 10. Replace "berg" or "Berg" by "0.5" (last field, ice cover)
sed 's/[bB]erg/0.5/g' log_rosette.tmp9 > log_rosette.tmp10

# 11. Remove "cast"
sed 's/cast//g' log_rosette.tmp10 > log_rosette.tmp11

# ------------------------------------------------------------------------------------------
# Now we are ready to format the header line with awk
# ------------------------------------------------------------------------------------------

cp log_rosette.tmp11 ${fxls}.nohead.txt

awk 'BEGIN{printf "day,month,year,H,M,dayUTC,monthUTC,yearUTC,HUTC,MUTC,SUTC,latNdegr,latNmin,lonWdegr,lonWmin,cap,stn,cast,zbot,Wdir,Wspeed,Tair,Twater,Patm,RH,Ice\n"} {print}' log_rosette.tmp11 | cat > ${fxls}.txt

# Tell us that things went well and print file to check

echo "File ${fxls}.txt has been written"
echo ""
cat ${fxls}.txt
