clear   ;close all;clc;
%% parameters
thres=0.99;gap=10;  
%% input 
event = 1;
% event =1: linetrip; 2:generator trip; 3: short circuit; 4: line to line faults 5:load change;
% 6: capacitor bank switch on; 7: motor off; 
for d=1:7
%% dictionary
dic=d;% 1 is line trip, 2 is generator trip , 3 is short circuit and 4 line to line, 5:load change, 6, capacitor banking, 7, motor off
%% free noise
noise=0;
noise_level=0;
%% output
[rank_event(d,:),rank_dic{1,d},min_angle(d,:), index_small(d,:),E1,E2{1,d}, Vk_event, Vk_dic]=compare_subangle(thres,gap,event,dic,noise,noise_level);
end


for m=1:size(min_angle,2)
    MingAngle(m)= min(min_angle(1:5,m)) ;
    index_min=find(min_angle(1:5,m)==MingAngle(m));
    ind_min(m)=index_min(1);
    if ind_min(m)==1
        if E1(m)>0.6
           ind_min(m)=2;
        end
    elseif ind_min(m)==2
        if E1(m)<=0.6
           ind_min(m)=1;
        end
    end
end

 