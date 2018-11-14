


var1=i2_toaii.

for ((y=1986;y<=2009;y++))
do
    for ((m=1;m<=12;m++))
        do
            for ((d=1;d<=31;d++))
                do
                    for ((h=0;h<=21;h=h+3))
                        do
i=${y:1}
yy=${i:1}
a=$((100+$m))
mm=${a:1}
b=$((100+$d))
dd=${b:1}
c=$((100+$h))
hh=${c:1}


mv /aos/home/syin/Desktop/TOA/${var1}${yy}${mm}${dd}${hh} /aos/home/syin/TOA1984to2009
done
done
done
done    
