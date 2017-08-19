clear all; clc; close all;
%%%%%%%%%%%%%%% this file  is used to generate data for different types of
%%%%%%%%%%%%%%% events at different locations
data16m_detail;  
event=2;%4:line trip; 0: three phase 2:line to line faults
sw_con(2,6)=event;% 4--line trip 0--short circuit, 6-- load change
path='H:\identification_data\Line_to_Line_86';
current_path='H:\';
for inde=36:86 
    
%%%%%%%%%%%%%%%%%%% it is used to change event type and locations    for

%%%%%%%%%%%%%%%%%%% Three phase short circuit
if event==0
    sw_con(2,2:3)=line(inde,1:2);
    bus1=line(inde,1);
    bus2=line(inde,2);   
%%%%%%%%%%%%%%%%%% for line outage 
elseif event==2
    sw_con(2,2:3)=line(inde,1:2); 
    bus1=line(inde,1);
    bus2=line(inde,2);
%%%%%%%%%%%%%%%%%% for line outage 
elseif event==4
    sw_con(2,2:3)=line(inde,1:2); 
    bus1=line(inde,1);
    bus2=line(inde,2);
%%%%%%%%%%%%%%%%%%% load change
elseif event==6
    index=[find(line(:,1)==inde);find(line(:,2)==inde)];
    sw_con(2,2:3)=line(index(1),1:2);
    load_con(1,1)=inde;
	lmod_con(1,2)=inde; 
end
[bus_v,bus_freq,theta,mac_ang,mac_spd]=s_simulwt(sw_con,load_con,lmod_con,bus,line);
cd(path);
% fn = sprintf('linetrip_%s-%s.mat',num2str(bus1),num2str(bus2));
% fn = sprintf('Three_Phase%s-%s.mat',num2str(bus1),num2str(bus2));
fn = sprintf('Line_to_Line_%s-%s.mat',num2str(bus1),num2str(bus2));
save(fn) 
cd(current_path)
end
% cd(current_path)
 