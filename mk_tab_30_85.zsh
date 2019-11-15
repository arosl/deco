#!/usr/bin/env bash

rm -f tab_*.ps tab_*.pdf

outfile=table_30_85.pdf
GFlow=0.3
GFhigh=0.85

time12b="25,30,35,40,50,60,70,80,90,100,120,140"
time12c="40,50,60,70,80,90,100,120,140,160,180,210"
time12d="60,70,80,90,100,120,140,160,180,210,240,270"
time12e="90,100,120,140,160,180,210,240,270,300,330,360"

dgas6=o2
dgas21=ean50
dgas36=ean32
dgas57=21/45@57
dgas72=18/55@72

bgas30=ean32
bgas36=30/30
bgas39=25/35
bgas57=18/45
bgas72=15/55
bgas90=10/65

opt="-f $GFlow,$GFhigh -r 9.0"
num=0

num=$(($num+1))
file="tab_$num.ps"

#ps header
./deco --hdr > $file

#tables
echo "($bgas30 1/2: 15m - 21m) heading" >> $file
./deco $opt -d 15 -g "$bgas30,[$dgas6]" -t $time12e -p >> $file
./deco $opt -d 16 -g "$bgas30,[$dgas6]" -t $time12e -p >> $file
./deco $opt -d 17 -g "$bgas30,[$dgas6]" -t $time12e -p >> $file
./deco $opt -d 18 -g "$bgas30,[$dgas6]" -t $time12e -p >> $file
./deco $opt -d 21 -g "$bgas30,[$dgas6]" -t $time12d -p >> $file

#ps join
./deco --jnr >> $file

#tables
echo "($bgas30 2/2: 24m - 33m) heading" >> $file
./deco $opt -d 24 -g "$bgas30,[$dgas6]" -t $time12c -p >> $file
./deco $opt -d 27 -g "$bgas30,[$dgas6]" -t $time12c -p >> $file
./deco $opt -d 30 -g "$bgas30,[$dgas6]" -t $time12c -p >> $file
./deco $opt -d 33 -g "$bgas30,[$dgas6]" -t $time12c -p >> $file

#ps footer
./deco --ftr >> $file

num=$(($num + 1))
file="tab_$num.ps"

#ps header
./deco --hdr > $file

#tables
echo "($bgas36 1/2: 15m - 21m) heading" >> $file
./deco $opt -d 15 -g "$bgas36,[$dgas6]" -t $time12e -p >> $file
./deco $opt -d 16 -g "$bgas36,[$dgas6]" -t $time12e -p >> $file
./deco $opt -d 17 -g "$bgas36,[$dgas6]" -t $time12e -p >> $file
./deco $opt -d 18 -g "$bgas36,[$dgas6]" -t $time12e -p >> $file
./deco $opt -d 21 -g "$bgas36,[$dgas6]" -t $time12d -p >> $file

#ps join
./deco --jnr >> $file

#tables
echo "($bgas36 2/2: 24m - 33m) heading" >> $file
./deco $opt -d 24 -g "$bgas36,[$dgas6]" -t $time12c -p >> $file
./deco $opt -d 27 -g "$bgas36,[$dgas6]" -t $time12c -p >> $file
./deco $opt -d 30 -g "$bgas36,[$dgas6]" -t $time12c -p >> $file
./deco $opt -d 33 -g "$bgas36,[$dgas6]" -t $time12c -p >> $file

#ps footer
./deco --ftr >> $file

case `uname` in
    Darwin)
        for f in ./tab_*.ps; do pstopdf $f; done
        /System/Library/Automator/Combine\ PDF\ Pages.action/Contents/Resources/join.py --output $outfile tab_*.pdf
        ;;
    Linux)
        for f in ./tab_*.ps; do ps2pdf -sPAPERSIZE=a4 $f; done
        pdfunite tab_*.pdf $outfile
        ;;
    *)
        echo "Only MacOS and Linux will work ¯\_(ツ)_/¯"
        exit 1
        ;;
esac

rm -f tab_*.ps tab_*.pdf