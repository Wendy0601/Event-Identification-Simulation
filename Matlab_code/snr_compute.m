clc; clear all; close all;
data=1;
deltat=0;
delay=0; 
for data=1:44
% load data
if data==1
load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2015-10-17 Line 328 Tripped and Reclosed 3 times and outage\2015_10_17_Line_328_1st_Trip_and_Reclose.mat', 'V_pm');
t=400 ;
end
if data==2 % after closed
load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2016-06-15 NY Line Tripping Event\2016_06_15_NY_Tripping_Event.mat', 'V_pm','I_pm')
t=356+deltat;%; %328+deltat this is a fault at 464 
end
if data==3
load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2016-07-07 Line T172 Fault and S171 Tripped\2016_07_07_S171_Trip.mat', 'V_pm','I_pm')
t=438+delay;
end
if data==4 % no close
  load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2016-08-06 Line 323,325 Trip and Reclose\2016_08_06_323_325_Trip_Reclose.mat', 'V_pm','I_pm')
t=463+delay;
end
if data==5
   load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2016-08-11 Line 364 Trip\2016_08_11_364_Trip.mat','I_pm','V_pm','V_pa');
%   V_pm=V_pa;
    t=508+delay;% the fault time
end
if data==6
load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2016-08-12 Line 359, 3271 Trip and Reclose\2016_08_12_359_3271_Trip_Reclose.mat', 'V_pm','I_pm')
t=292+delay;        
end
if data==7
 load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2016-07-22 Line T172S Tripped\2016_07_22_T172S_Trip.mat', 'V_pm','I_pm')
t=462+delay;
  
end
if data==8 
load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2015-10-14 Line 3557 Trip\2015_10_14_Line_3557_Trip.mat','I_pm','V_pm');
 t=241+deltat; %%% I do not know what happens there ?
  
end 
if data==9
load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping Events\2016-10-06 Line 397 Trip\Line_397_Trip.mat', 'V_pm')
t=197;%103+deltat;  
end 
if data==10
 load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping Events\2016-10-05 3216 Trip\Line_3216_Trip.mat', 'V_pm')
t=118+deltat; %% it is not directly measurements
end
if data==11
load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping Events\2016-10-11 240-601 and 336 Trip\Line_240-601_and_336_Trip.mat', 'V_pm')
t=300+deltat;    % the first line trip
end
if data==12
load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2015-10-17 Line 328 Tripped and Reclosed 3 times and outage\2015_10_17_Line_328_3nd_Trip_and_Reclose.mat', 'V_pm')
 t=406+deltat; % it is reclosed 
end
if data==13
load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2015-10-17 Line 328 Tripped and Reclosed 3 times and outage\2015_10_17_Line_328_2nd_Trip_and_Reclose.mat', 'V_pm')
t=450+deltat;
end
if data==14
load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2015-10-17 Line 328 Tripped and Reclosed 3 times and outage\2015_10_17_Line_328_Outage.mat', 'V_pm');    
t=329+deltat;% this time it is outaged, current from 768 to 0
end 
if data==15% no close
  load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2015-10-17 Line 354 Trip\2015_10_17_Line354_Trip.mat','I_pm','V_pm'); 
  t=1065+deltat;
end
if data==16 % no close
    load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2015-11-27 Line 368 Trip and Close\2015_11_27_Line_368_Trip.mat', 'V_pm','I_pm')
    t=414+deltat;
end
if data==17
load('H:\Labfile_realdata\Files\Line Tripping Events\2017_03_17_Line_314_Trip_Reclose\2017_03_17_Line_314_Trip_Reclose.mat', 'V_pm');
t=149+deltat;
end
if data==18
load('H:\Labfile_realdata\Files\Line Tripping Events\2017_07_10_Line_3001_Trip_Close\2017_07_10_Line_3001_Trip.mat', 'V_pm');
t=87+deltat; % it is reclosed 
end
if data==19
load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping Events\2016-10-11 240-601 and 336 Trip\Line_240-601_and_336_Trip.mat', 'V_pm')
t=553+deltat; % the second line trip
end

if data==20
    load('H:\Labfile_realdata\Bear_Swamp_Pump_Trip\Bear_Swamp_Pump_Trip\2017_06_10_Bear_Swamp_270_MW.mat', 'V_pm');
    t=84+deltat;
end
if data==21
load('H:\Labfile_realdata\ISO_New_England_Files\Pump Storage (Load) Tripping Events\NE_Northfield_Mountain_Trip_1_200MW.mat', 'V_pm')
t=376+deltat;% THIS IS THE 200MW load
end

if data==22
load('H:\Labfile_realdata\ISO_New_England_Files\Pump Storage (Load) Tripping Events\NE_Northfield_Mountain_Trip_2_200MW.mat', 'V_pm')
t=703+deltat;% THIS IS THE 200MW load
end 

if data==23
load('H:\Labfile_realdata\ISO_New_England_Files\Pump Storage (Load) Tripping Events\NE_Northfield_Mountain_Trip_3_200MW.mat', 'V_pm')
t=384+deltat;% THIS IS THE 200MW load
end 

if data==24
load('H:\Labfile_realdata\ISO_New_England_Files\Pump Storage (Load) Tripping Events\NE_Northfield_Mountain_Trip_4_200MW.mat', 'V_pm')
t=254+deltat;% THIS IS THE 200MW load
end 

if data==25
load('H:\Labfile_realdata\ISO_New_England_Files\Pump Storage (Load) Tripping Events\NE_Northfield_Mountain_Trip_5_200MW.mat', 'V_pm')
t=313+deltat;% THIS IS THE 200MW load
end 

if data==26
load('H:\Labfile_realdata\ISO_New_England_Files\Pump Storage (Load) Tripping Events\NE_Northfield_Mountain_Trip_6_200MW.mat', 'V_pm')
t=327+deltat;% THIS IS THE 200MW load
end 

if data==27
load('H:\Labfile_realdata\ISO_New_England_Files\Pump Storage (Load) Tripping Events\NE_Northfield_Mountain_Trip_7_200MW.mat', 'V_pm')
t=515+deltat;% THIS IS THE 200MW load
end 

if data==28
load('H:\Labfile_realdata\ISO_New_England_Files\Pump Storage (Load) Tripping Events\NE_Northfield_Mountain_Trip_8_200MW.mat', 'V_pm')
t=259+deltat;% THIS IS THE 200MW load
end 

if data==29
    load('H:\Labfile_realdata\Bear_Swamp_Pump_Trip\Bear_Swamp_Pump_Trip\2017_06_09_Bear_Swamp_280_MW.mat', 'V_pm');
    t=96+deltat;
end


if data==30
    load('H:\Labfile_realdata\Bear_Swamp_Pump_Trip\Bear_Swamp_Pump_Trip\2017_06_11_Bear_Swamp_270_MW.mat', 'V_pm');
    t=77+deltat;
end
if data==31
    load('H:\Labfile_realdata\Bear_Swamp_Pump_Trip\Bear_Swamp_Pump_Trip\2017_06_07_Bear_Swamp_280_MW.mat', 'V_pm');
    t=106+deltat;
end

if data==32
    load('H:\Labfile_realdata\Bear_Swamp_Pump_Trip\Bear_Swamp_Pump_Trip\2017_06_08_Bear_Swamp_280_MW.mat', 'V_pm');
    t=130+deltat;
end  


if data==33
load('H:\Labfile_realdata\ISO_New_England_Files\Generator Tripping Events\2016-03-15 NE Phase II HVDC Pole 2 Trip - 650 MW\NE_Phase_II_Pole_2_Trip_650MW.mat', 'V_pm')
t=68+deltat;% THIS IS THE HVDC pole 2 trip 650MW
end 

if data==34
load('H:\Labfile_realdata\ISO_New_England_Files\Generator Tripping Events\2016-05-15 NE Milstone 3 Tripped - 850 MW\NE_Millstone_3_Trip_850MW.mat', 'V_pm')
t=1056+deltat;% THIS IS THE 850MW trip this is the external generator trip
end

if data==35
load('H:\Labfile_realdata\ISO_New_England_Files\Generator Tripping Events\2016-07-06 NY Indian Point 2 Trip - 1000 MW\NY_Indian_Point_2_Trip_1000MW.mat', 'V_pm')
t=1133+deltat;% THIS IS THE 1000MW trip outside
end

if data==36
load('H:\Labfile_realdata\ISO_New_England_Files\Generator Tripping Events\2016-07-08 NE East Shore NHHB Trip - 370 MW\NE_BPT_Energy_Trip_370MW.mat', 'V_pm')
t=354+deltat;% THIS IS THE 370 MWtrip
end

if data==37
load('H:\Labfile_realdata\ISO_New_England_Files\Generator Tripping Events\2016-07-16 NE BPT Energy Singer T1X Trip - 130, 250 MW\NE_BPT_Energy_Trip_130MW_250MW.mat', 'V_pm')
t=408+deltat; % THIS IS THE 130 MWtrip
end

if data==38
load('H:\Labfile_realdata\ISO_New_England_Files\Generator Tripping Events\2016-07-16 NE BPT Energy Singer T1X Trip - 130, 250 MW\NE_BPT_Energy_Trip_130MW_250MW.mat', 'V_pm')
t=824+deltat;% THIS IS THE 250 MWtrip
end

if data==39
load('H:\Labfile_realdata\ISO_New_England_Files\Generator Tripping Events\2016-07-17 NE BPT Energy Singer T1X Trip - 180 MW\NE_BPT_Energy_Trip_180MW.mat', 'V_pm')
t=121+deltat; % This is the 180 MW trip
end

if data==40
load('H:\Labfile_realdata\ISO_New_England_Files\Generator Tripping Events\2016-08-04 NE BPT Energy Singer T1X Trip - 450 MW\NE_BPT_Energy_Trip_450MW.mat', 'V_pm')
t=308+deltat; % This is the 450MW trip
end

if data==41
load('H:\Labfile_realdata\ISO_New_England_Files\Generator Tripping Events\2016-08-31 PJM Salem 2 Trip - 1000 MW\PJM_Salem_2_Trip_1000MW.mat', 'V_pm')
t=2110+deltat; % this is the 1000MW trip
end

if data==42
load('H:\Labfile_realdata\ISO_New_England_Files\Generator Tripping Events\2016-10-05 NE BPT Energy Singer T1X Trip - 535 MW\NE_BPT_Energy_Trip_535MW.mat', 'V_pm')
t=186+deltat;% this is the 535 MW
end
if data==43
load('H:\Labfile_realdata\ISO_New_England_Files\Generator Tripping Events\2016-03-02 NE Seabrook Trip - 1260 MW\NE_Seabrook_Trip_1260MW.mat', 'V_pm')
t=328+deltat;% THIS IS THE seabrook trip 1260MW
end 

if data==44
load('H:\Labfile_realdata\ISO_New_England_Files\Generator Tripping Events\2016-03-10 NE Phase II HVDC Pole 1 Trip - 650 MW\NE_Phase_II_Pole_1_Trip_650MW.mat', 'V_pm')
t=317+deltat;% After the second osillation, it is tripped , thus we choose the second one as the starting pointTHIS IS THE HVDC pole I trip 650MW
end  

X=[];X1=[];snr=[];
% select the steady datasets
X=V_pm(1:t-10,:);
[row,col]=size(X);
% noise filter
M=25;  
for i=1:(row-M)
    X1(i,:)=mean(X(i:(i+M),:));
end
noise=X(M+1:row,:)-X1;
for i=1:col
    snr(i)=10*log10(norm(X1(:,i))^2/norm(noise(:,i))^2);
end
snr(find(isnan(snr)))=[];
SNR_total(data)=mean(snr);% the mean of 136 channels
end
mean(SNR_total)

% SNR of TVA PMU datasets
% clc; clear all; close all;
% load('H:\Labfile_realdata\All_data.mat', 'Mag_data');
% t=1337;
% X=Mag_data(1:t-10,1:5);
% [row,col]=size(X);
% % noise filter
% M=25;  
% for i=1:(row-M)
%     X1(i,:)=mean(X(i:(i+M),:));
% end
% noise=X(M+1:row,:)-X1;
% for i=1:col
%     snr(i)=10*log10(norm(X1(:,i))^2/norm(noise(:,i))^2);
% end
% snr(find(isnan(snr)))=[];

