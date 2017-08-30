function [ V_pm ] = TrainData_linetrip( data )
if data==1
load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2015-10-17 Line 328 Tripped and Reclosed 3 times and outage\2015_10_17_Line_328_Outage.mat', 'V_pm');    
t=329+deltat;% this time it is outaged, current from 768 to 0
end 



if data==2% no close
  load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2015-10-17 Line 354 Trip\2015_10_17_Line354_Trip.mat','I_pm','V_pm'); 
  t=1065+deltat;
end

if data==3 % no close
    load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2015-11-27 Line 368 Trip and Close\2015_11_27_Line_368_Trip.mat', 'V_pm','I_pm');
    t=414+deltat;
end

if data==4
load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping Events\2016-10-11 240-601 and 336 Trip\Line_240-601_and_336_Trip.mat', 'V_pm');
t=553+deltat; % the second line trip
end

if data==5
load('H:\Labfile_realdata\Files\Line Tripping Events\2017_03_17_Line_314_Trip_Reclose\2017_03_17_Line_314_Trip_Reclose.mat', 'V_pm');
t=151+deltat;
end

if data==6
load('H:\Labfile_realdata\Files\Line Tripping Events\2017_07_10_Line_3001_Trip_Close\2017_07_10_Line_3001_Trip.mat', 'VoltageSignal');
t=91+deltat; % it is reclosed 
end

end

