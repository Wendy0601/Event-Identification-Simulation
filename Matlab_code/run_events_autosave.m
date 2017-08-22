clear ; clc; close all;
%%%%%%%%%%%%%%% this file  is used to generate data for different types of
%%%%%%%%%%%%%%% events at different locations
data16m_detail;  
event=5; event_num=33;
% event=7;event_num=33;
%4:line trip event_num=86; 0: three phase event_num=68 2:line to line faults event_num=68 6: load change(notice only load_con the changed bus); event_num=68 7: capacitor bank switch on
sw_con(2,6)=event;% 4--line trip 0--short circuit, 6-- load change
path='H:\identification_data\Motor_Off_86';
current_path='H:\';

for inde=1:event_num
    
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
    %%%%%%%%%%%%% Switch off motor load
    elseif event==5 
        load_bus=[1 3 4  7 8  9  15 16 18 20 21 23 24 25 26 27 28 29 33 37 39 40 41 42 44 45 46 47 48 49 50 51 52];
        inde_c=load_bus(inde);
        line_id =find(line(:,1)==inde_c);% notice only the bus needed to change is set to be non-conforming bus
        if isempty(line_id)
            line_id=find(line(:,2)==inde_c);
            sw_con(2,3)=line(line_id(1),1);
            sw_con(2,2)=line(line_id(1),2);
        else
            sw_con(2,2:3)=line(line_id(1),1:2);
        end 
    %%%%%%%%%%%%%%%%%%% load change
    elseif event==6
        line_id =find(line(:,1)==inde);% notice only the bus needed to change is set to be non-conforming bus
        if isempty(line_id)
            line_id=find(line(:,2)==inde);
            sw_con(2,3)=line(line_id(1),1);
            sw_con(2,2)=line(line_id(1),2);
        else
            sw_con(2,2:3)=line(line_id(1),1:2);
        end
        load_con(1,1)=inde;
        lmod_con(1,2)=inde; 
    elseif event==7
        load_bus=[1 3 4  7 8  9  15 16 18 20 21 23 24 25 26 27 28 29 33 37 39 40 41 42 44 45 46 47 48 49 50 51 52];
        inde_c=load_bus(inde);
        inde_bus=load_bus(inde_c);
        line_id =find(line(:,1)==inde_bus);
        if isempty(line_id)
            line_id=find(line(:,2)==inde_c);
            sw_con(2,3)=line(line_id(1),1);
            sw_con(2,2)=line(line_id(1),2);
        else
            sw_con(2,2:3)=line(line_id(1),1:2);
        end

    end
    [bus_v,bus_freq,theta,mac_ang,mac_spd,flag_error]=s_simulwt(sw_con,load_con,lmod_con,bus,line);
    if flag_error==1 
        continue;
    else
        cd(path);
        % fn = sprintf('linetrip_%s-%s.mat',num2str(bus1),num2str(bus2));
        % fn = sprintf('Three_Phase%s-%s.mat',num2str(bus1),num2str(bus2));
    %     fn = sprintf('Line_to_Line_%s-%s.mat',num2str(bus1),num2str(bus2));
    %     fn = sprintf('Capacitor_Bank_%s.mat',num2str(load_bus(inde)));
          fn = sprintf('Motor_Off_%s.mat',num2str(load_bus(inde)));
        save(fn) 
        cd(current_path)
    end
end
% cd(current_path)
 