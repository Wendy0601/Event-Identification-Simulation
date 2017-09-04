clear   ;close all;clc;
%% parameters
thres=0.99;gap=2;  
% root='H:\identification_data\condition_a';
% root='H:\identification_data\condition_h';
root='H:\identification_data\condition_b';
%% input 
for event =1
    rank_event=[];
    rank_dic={};
    min_angle=[];
    index_small=[];
    E1=[];
    E2=[];
    allnames={};
    ind_min=[];
% event: 1 is line trip, 2 is generator trip , 3 is short circuit and 4 is load change
for d=1 :4
%% dictionary
dic=d;% 1 is line trip, 2 is generator trip , 3 is short circuit and 4 line to line, 5:load change, 6, capacitor banking, 7, motor off
%% free noise
noise=1;
noise_level=40;
%% output
[rank_event(d,:),rank_dic{1,d},min_angle(d,:), index_small(d,:),E1,E2{1,d},allnames]=compare_subangle(root, thres,gap,event,dic,noise,noise_level);
end
 

for m=1:size(min_angle,2)
    MingAngle(m)= min(min_angle(1:d,m)) ;
    index_min=find(min_angle(1:d,m)==MingAngle(m));
    ind_min(m)=index_min(1);
    if ind_min(m)==1 
        if (E1(m)>0.67)  % the max of linetrip events
            ind_min(m)=2;
        end
    elseif ind_min(m)==2
        if E1(m)<=0.5 % the min of generator trip events
           ind_min(m)=1;
        end 
    end
end

IAR(event)=numel(find(ind_min==event))/numel(ind_min) ;
end