function [ V_pm ] = TrainData_faults( data )

if data==1
load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2016-08-12 Line 359, 3271 Trip and Reclose\2016_08_12_359_3271_Trip_Reclose.mat', 'V_pm','I_pm');
t=292+delay;    
    
end


if data==2
load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping Events\2016-10-06 Line 397 Trip\Line_397_Trip.mat', 'V_pm');
t=102+delay; %% it is directly measurements and there is a fault there only one side is meausred
end   

if data==3
load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping Events\2016-10-14 349 Fault\Line_349_Fault.mat', 'V_pm');
t=252+delay; % there is a large fault
end
 
if data==4
load('H:\Labfile_realdata\Files\Line Tripping Events\2017_01_23_Line_357_Trip_Reclose\2017_01_23_Line_357_Trip_Reclose.mat', 'V_pm');
t=117+delay;
end

if data==5
    load('H:\Labfile_realdata\Files\Line Tripping Events\2017_04_19_Line_339_Trip_Reclose\2017_04_19_Line_339_Trip_Reclose.mat', 'V_pm');
    t=110+delay;
end

if data==6
    load('H:\Labfile_realdata\Files\Line Tripping Events\2017_01_23_Line_357_Trip_Reclose\2017_01_23_Line_357_Trip_Reclose.mat', 'V_pm');
    t=951+delay;
end

if data==7
    load('H:\Labfile_realdata\Files\Line Tripping Events\2017_06_24_Line_329_Trip_Reclose\2017_06_24_Line_329_Trip_Reclose.mat', 'V_pm');
    t=280+delay;
end

    
 
 
end

