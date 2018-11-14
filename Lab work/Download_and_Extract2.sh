

#https://isccp.giss.nasa.gov/outgoing/FLUX/TOA/i2_toaii.1983mmdd00.tar

var1=https://isccp.giss.nasa.gov/outgoing/FLUX/TOA/i2_toaii.
var2=mmdd
var3=.tar
var4=i2_toaii.

for ((i=1984;i<=2009;i++))
do
for  ((j=0;j<=21;j=j+3))
do
a=$((100+$j))
j=${a:1}
wget ${var1}${i}${var2}${j}${var3}
# wget https://isccp.giss.nasa.gov/outgoing/FLUX/SRF/i2_srfii.1991mmdd00.tar
done
done

for ((m=1984;m<=2009;m++))
do
for  ((n=0;n<=21;n=n+3))
do
b=$((100+$n))
n=${b:1}
tar -xvf ${var4}${m}${var2}${n}${var3}
# tar –xvf i2_srfii.1991mmdd00.tar
done
done


'
for ((i=2008;i<=2008;i++))
do
for  ((j=12;j<=21;j=j+3))
do
a=$((100+$j))
j=${a:1}
wget ${var1}${i}${var2}${j}${var3}
# wget https://isccp.giss.nasa.gov/outgoing/FLUX/SRF/i2_srfii.1991mmdd00.tar
done
done
'

'
for ((m=1983;m<=1986;m++))
do
for  ((n=0;n<=21;n=n+3))
do
b=$((100+$n))
n=${b:1}
tar -xvf ${var4}${m}${var2}${n}${var3}
# tar –xvf i2_srfii.1991mmdd00.tar
done
done


for ((m=2008;m<=2008;m++))
do
for  ((n=12;n<=21;n=n+3))
do
b=$((100+$n))
n=${b:1}
tar -xvf ${var4}${m}${var2}${n}${var3}
# tar –xvf i2_srfii.1991mmdd00.tar
done
done
'



#:<<!
#[root@localhost sh]# var1=/etc/
#[root@localhost sh]# var2=yum.repos.d/
#[root@localhost sh]# var3=${var1}${var2}
#[root@localhost sh]# echo $var3
#/etc/yum.repos.d/

#var=$(printf "%03d" "$no") 

#a=$((1000+$i))
#echo ${a:1}

#!
