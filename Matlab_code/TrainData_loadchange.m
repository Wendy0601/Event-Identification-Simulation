function [ V_pm ] = TrainData_loadchange( data )

if data==1
load('H:\Labfile_realdata\ISO_New_England_Files\Pump Storage (Load) Tripping Events\NE_Northfield_Mountain_Trip_2_200MW.mat', 'V_pm');
t=703+deltat;% THIS IS THE 200MW load
end 
if data==2
    load('H:\Labfile_realdata\ISO_New_England_Files\Pump Storage (Load) Tripping Events\NE_Northfield_Mountain_Trip_6_200MW.mat', 'V_pm');
    t=327+deltat;% THIS IS THE 200MW load
end 

if data==3
    load('H:\Labfile_realdata\Bear_Swamp_Pump_Trip\Bear_Swamp_Pump_Trip\2017_06_07_Bear_Swamp_280_MW.mat', 'V_pm');
    t=106+deltat;
end

if data==4
    load('H:\Labfile_realdata\Bear_Swamp_Pump_Trip\Bear_Swamp_Pump_Trip\2017_06_08_Bear_Swamp_280_MW.mat', 'V_pm');
    t=130+deltat;
end

if data==5
    load('H:\Labfile_realdata\Bear_Swamp_Pump_Trip\Bear_Swamp_Pump_Trip\2017_06_09_Bear_Swamp_280_MW.mat', 'V_pm');
    t=96+deltat;
end

if data==6
    load('H:\Labfile_realdata\Bear_Swamp_Pump_Trip\Bear_Swamp_Pump_Trip\2017_06_10_Bear_Swamp_270_MW.mat', 'V_pm');
    t=84+deltat;
end

if data==7
    load('H:\Labfile_realdata\Bear_Swamp_Pump_Trip\Bear_Swamp_Pump_Trip\2017_06_11_Bear_Swamp_270_MW.mat', 'V_pm');
    t=77+deltat;
end

if data==8
load('H:\Labfile_realdata\ISO_New_England_Files\Pump Storage (Load) Tripping Events\NE_Northfield_Mountain_Trip_7_200MW.mat', 'V_pm')
t=515+deltat;% THIS IS THE 200MW load
end 

if data==9
load('H:\Labfile_realdata\ISO_New_England_Files\Pump Storage (Load) Tripping Events\NE_Northfield_Mountain_Trip_8_200MW.mat', 'V_pm')
t=259+deltat;% THIS IS THE 200MW load
end 
    

end

