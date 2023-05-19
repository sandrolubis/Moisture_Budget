#!/bin/bash

mkdir climo
mkdir eddy
mkdir output

cdo delete,month=2,day=29 /glade/p/univ/uchi0007/ERA5_Data/U/daily/u700/combined.nc u700.daily.nc
cdo delete,month=2,day=29 /glade/p/univ/uchi0007/ERA5_Data/V/daily/v700/combined.nc v700.daily.nc
cdo delete,month=2,day=29 /glade/p/univ/uchi0007/ERA5_Data/SH/daily/q700/combined.nc q700.daily.nc  

cdo ydaymean u700.daily.nc climo.u700.nc
cdo ydaymean v700.daily.nc climo.v700.nc
cdo ydaymean q700.daily.nc climo.q700.nc

cdo merge climo.*.nc climo/climo.nc
rm climo.u700.nc
rm climo.v700.nc
rm climo.q700.nc

cdo ydaysub u700.daily.nc -ydaymean u700.daily.nc u700.anom.nc
cdo ydaysub v700.daily.nc -ydaymean v700.daily.nc v700.anom.nc
cdo ydaysub q700.daily.nc -ydaymean q700.daily.nc q700.anom.nc

rm u700.daily.nc
rm v700.daily.nc
rm q700.daily.nc



