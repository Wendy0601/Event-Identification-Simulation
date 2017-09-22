clear  ;close all;clc;
%% parameters
thres=0.99;gap=0;  
retrain=0;
% root='H:\identification_data\condition_a';
root='H:\identification_data\condition_d';% this is test data \condition_c
if retrain==1
    root_dic='H:\identification_data\Retrain_dictionary';
else
    root_dic='H:\identification_data\Dic_current';
end

%% input 
for event =1:4
    rank_event=[];
    rank_dic={};
    min_angle=[];
    index_small=[];
    E1=[];
    E2=[];
    allnames={};
    ind_min=[];
% event: 1 is line trip, 2 is generator trip , 3 is short circuit and 4 is load change
    for d=1:4
        %% dictionary
        dic=d;% 1 is line trip, 2 is generator trip , 3 is short circuit and 4 line to line, 5:load change, 6, capacitor banking, 7, motor off
        %% free noise
        noise=0;
        noise_level=0;
        %% output
        [rank_event(:,d),rank_dic{1,d},min_angle(d,:), index_small(d,:),E1,E2{1,d},allnames]=Current_compare_subangle(root,root_dic, thres,gap,event,dic,noise,noise_level);
    end
 
    for m=1:size(min_angle,2)
        MingAngle(m)= min(min_angle(1:d,m)) ;
        index_min=find(min_angle(1:d,m)==MingAngle(m));
        ind_min(m)=index_min(1);
%         if ind_min(m)==1 
%             if (E1(m)>0.65)  % the max of linetrip events
%                 ind_min(m)=2;
%             end
%         elseif ind_min(m)==2
%             if E1(m)<=300 % the min of generator trip events
%                ind_min(m)=1;
%             end 
%         elseif ind_min(m)==3
%             if E1(m)<=55
%                 ind_min(m)=5;
%             end 
%         elseif ind_min(m)==5
%             if E1(m)>55
%                 ind_min(m)=3;
%             end       
%         end  
    end
 
IAR(event)=numel(find(ind_min==event))/numel(ind_min) ;
end
a=mean(rank_event)

 