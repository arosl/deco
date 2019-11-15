#

rm -f tab_*.ps tab_*.pdf

set time12a = "10,15,20,25,30,35,40,50,60,70,80,90"
set time12b = "25,30,35,40,50,60,70,80,90,100,120,140"
set time12c = "40,50,60,70,80,90,100,120,140,160,180,210"
set time12d = "60,70,80,90,100,120,140,160,180,210,240,270"
set time12e = "90,100,120,140,160,180,210,240,270,300,330,360"

set dgas6  = o2
set dgas21 = ean50
set dgas36 = ean32
set dgas57 = 21/45@57
set dgas72 = 18/55@72

set bgas30 = ean32
set bgas39 = 25/35
set bgas57 = 18/45
set bgas72 = 15/55
set bgas90 = 10/65

set opt = "-f 0.5,0.9"

set num = 0;

set num = `expr $num + 1`; set file = `printf "tab_%2.2i.ps" $num`

./deco --hdr > $file
echo "(Air 1/2: 15m - 27m) heading" >> $file
./deco $opt -d 15 -g air,[$dgas6] -t $time12e -p >> $file
./deco $opt -d 18 -g air,[$dgas6] -t $time12d -p >> $file
./deco $opt -d 21 -g air,[$dgas6] -t $time12c -p >> $file
./deco $opt -d 24 -g air,[$dgas6] -t $time12c -p >> $file
./deco $opt -d 27 -g air,[$dgas6] -t $time12b -p >> $file

./deco --jnr >> $file

echo "(Air 2/2: 30m - 39m) heading" >> $file
./deco $opt -d 30 -g air,[$dgas6] -t $time12b -p >> $file
./deco $opt -d 33 -g air,[$dgas6] -t $time12b -p >> $file
./deco $opt -d 36 -g air,[$dgas6] -t $time12a -p >> $file
./deco $opt -d 39 -g air,[$dgas6] -t $time12a -p >> $file
./deco --ftr >> $file

set num = `expr $num + 1`; set file = `printf "tab_%2.2i.ps" $num`

./deco --hdr > $file
echo "(EAN32 1/2: 15m - 24m) heading" >> $file
./deco $opt -d 15 -g $bgas30,[$dgas6] -t $time12e -p >> $file
./deco $opt -d 18 -g $bgas30,[$dgas6] -t $time12e -p >> $file
./deco $opt -d 21 -g $bgas30,[$dgas6] -t $time12e -p >> $file
./deco $opt -d 24 -g $bgas30,[$dgas6] -t $time12d -p >> $file
./deco $opt -d 24 -g $bgas30,$dgas21,$dgas6 -t $time12e -p >> $file

./deco --jnr >> $file

echo "(EAN32 2/2: 27m - 33m) heading" >> $file
./deco $opt -d 27 -g $bgas30,[$dgas6] -t $time12d -p >> $file
./deco $opt -d 27 -g $bgas30,$dgas21,$dgas6 -t $time12e -p >> $file
./deco $opt -d 30 -g $bgas30,[$dgas6] -t $time12c -p >> $file
./deco $opt -d 30 -g $bgas30,$dgas21,$dgas6 -t $time12e -p >> $file
./deco $opt -d 33 -g $bgas30,[$dgas6] -t $time12c -p >> $file
./deco --ftr >> $file

set num = `expr $num + 1`; set file = `printf "tab_%2.2i.ps" $num`

./deco --hdr > $file
echo "(Short trimix 1/10: 30m - 36m) heading" >> $file
./deco $opt -d 30 -g $bgas39,[$dgas6] -t $time12a -p >> $file
./deco $opt -d 33 -g $bgas39,[$dgas6] -t $time12a -p >> $file
./deco $opt -d 33 -g $bgas39,$dgas21,[$dgas6] -t $time12a -p >> $file
./deco $opt -d 36 -g $bgas39,$dgas21,[$dgas6] -t $time12a -p >> $file

./deco --jnr >> $file

echo "(Short trimix 2/10: 39m - 48m) heading" >> $file
./deco $opt -d 39 -g $bgas57,$dgas21,[$dgas6] -t $time12a -p >> $file
./deco $opt -d 42 -g $bgas57,$dgas21,[$dgas6] -t $time12a -p >> $file
./deco $opt -d 48 -g $bgas57,$dgas21,[$dgas6] -t $time12a -p >> $file
./deco --ftr >> $file

set num = `expr $num + 1`; set file = `printf "tab_%2.2i.ps" $num`

./deco --hdr > $file
echo "(Short trimix 3/10: 54m) heading" >> $file
./deco $opt -d 54 -g $bgas57,[$dgas36],[$dgas21],[$dgas6] -t $time12a -p >> $file
./deco $opt -d 54 -g $bgas57,$dgas21 -t $time12a -p >> $file

./deco --jnr >> $file

echo "(Short trimix 4/10: 60m) heading" >> $file
./deco $opt -d 60 -g $bgas72,[$dgas36],[$dgas21],[$dgas6] -t $time12a -p >> $file
./deco --ftr >> $file

set num = `expr $num + 1`; set file = `printf "tab_%2.2i.ps" $num`

./deco --hdr > $file
echo "(Short trimix 5/10: 66m) heading" >> $file
./deco $opt -d 66 -g $bgas72,[$dgas36],[$dgas21],[$dgas6] -t $time12a -p >> $file

./deco --jnr >> $file

echo "(Short trimix 6/10: 72m) heading" >> $file
./deco $opt -d 72 -g $bgas90,[$dgas36],[$dgas21],[$dgas6] -t $time12a -p >> $file
./deco --ftr >> $file

set num = `expr $num + 1`; set file = `printf "tab_%2.2i.ps" $num`

./deco --hdr > $file
echo "(Short trimix 7/10: 78m) heading" >> $file
./deco $opt -d 78 -g $bgas90,[$dgas36],[$dgas21],[$dgas6] -t $time12a -p >> $file

./deco --jnr >> $file

echo "(Short trimix 8/10: 84m) heading" >> $file
./deco $opt -d 84 -g $bgas90,$dgas57,[$dgas36],[$dgas21],$dgas6 -t $time12a -p >> $file
./deco --ftr >> $file

set num = `expr $num + 1`; set file = `printf "tab_%2.2i.ps" $num`

./deco --hdr > $file
echo "(Short trimix 9/10: 90m) heading" >> $file
./deco $opt -d 90 -g $bgas90,$dgas57,[$dgas36],[$dgas21],$dgas6 -t $time12a -p >> $file

./deco --jnr >> $file

echo "(Short trimix 10/10: 96m) heading" >> $file
./deco $opt -d 96 -g $bgas90,$dgas57,[$dgas36],[$dgas21],$dgas6 -t $time12a -p >> $file
./deco --ftr >> $file

set num = `expr $num + 1`; set file = `printf "tab_%2.2i.ps" $num`

./deco --hdr > $file
echo "(Long trimix 1/8: 30m - 39m) heading" >> $file
./deco $opt -d 30 -g $bgas39,$dgas21,$dgas6 -t $time12e -p >> $file
./deco $opt -d 33 -g $bgas39,$dgas21,$dgas6 -t $time12e -p >> $file
./deco $opt -d 36 -g $bgas57,$dgas21,$dgas6 -t $time12e -p >> $file
./deco $opt -d 39 -g $bgas57,[$dgas36],$dgas21,$dgas6 -t $time12e -p >> $file

./deco --jnr >> $file

echo "(Long trimix 2/8: 42m - 51m) heading" >> $file
./deco $opt -d 42 -g $bgas57,$dgas36,$dgas21,$dgas6 -t $time12e -p >> $file
./deco $opt -d 45 -g $bgas57,$dgas36,$dgas21,$dgas6 -t $time12e -p >> $file
./deco $opt -d 48 -g $bgas57,$dgas36,$dgas21,$dgas6 -t $time12e -p >> $file
./deco $opt -d 51 -g $bgas57,$dgas36,$dgas21,$dgas6 -t $time12e -p >> $file
./deco --ftr >> $file

set num = `expr $num + 1`; set file = `printf "tab_%2.2i.ps" $num`

./deco --hdr > $file
echo "(Long trimix 3/8: 54m - 60m) heading" >> $file
./deco $opt -d 54 -g $bgas57,$dgas36,$dgas21,$dgas6 -t $time12e -p >> $file
./deco $opt -d 57 -g $bgas72,$dgas36,$dgas21,$dgas6 -t $time12e -p >> $file
./deco $opt -d 60 -g $bgas72,$dgas36,$dgas21,$dgas6 -t $time12e -p >> $file

./deco --jnr >> $file

echo "(Long trimix 4/8: 63m - 66m) heading" >> $file
./deco $opt -d 63 -g $bgas72,$dgas36,$dgas21,$dgas6 -t $time12e -p >> $file
./deco $opt -d 66 -g $bgas72,[$dgas57],$dgas36,$dgas21,$dgas6 -t $time12e -p >> $file
./deco --ftr >> $file

set num = `expr $num + 1`; set file = `printf "tab_%2.2i.ps" $num`

./deco --hdr > $file
echo "(Long trimix 5/8: 69m - 72m) heading" >> $file
./deco $opt -d 69 -g $bgas72,$dgas57,$dgas36,$dgas21,$dgas6 -t $time12e -p >> $file
./deco $opt -d 72 -g $bgas90,$dgas57,$dgas36,$dgas21,$dgas6 -t $time12e -p >> $file

./deco --jnr >> $file

echo "(Long trimix 6/8: 75m - 77m) heading" >> $file
./deco $opt -d 75 -g $bgas90,$dgas57,$dgas36,$dgas21,$dgas6 -t $time12e -p >> $file
./deco $opt -d 78 -g $bgas90,$dgas57,$dgas36,$dgas21,$dgas6 -t $time12e -p >> $file
./deco --ftr >> $file

set num = `expr $num + 1`; set file = `printf "tab_%2.2i.ps" $num`

./deco --hdr > $file
echo "(Long trimix 7/8: 81m - 84m) heading" >> $file
./deco $opt -d 81 -g $bgas90,$dgas57,$dgas36,$dgas21,$dgas6 -t $time12e -p >> $file
./deco $opt -d 84 -g $bgas90,$dgas57,$dgas36,$dgas21,$dgas6 -t $time12e -p >> $file

./deco --jnr >> $file

echo "(Long trimix 8/8: 87m - 90m) heading" >> $file
./deco $opt -d 87 -g $bgas90,$dgas72,$dgas57,$dgas36,$dgas21,$dgas6 -t $time12e -p >> $file
./deco $opt -d 90 -g $bgas90,$dgas72,$dgas57,$dgas36,$dgas21,$dgas6 -t $time12e -p >> $file
./deco --ftr >> $file

foreach file ( tab_*.ps )
  ps2pdf -sPAPERSIZE=a4 $file
end

pdfunite tab_*.pdf table_50_90.pdf

rm -f tab_*.ps tab_*.pdf
