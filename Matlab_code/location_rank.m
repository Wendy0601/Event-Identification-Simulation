%%%%%%%% test performance 
clear all  ; close all; clc; 
event=4; 
kmax=50;
thres=0.99;
current_path=pwd;
load('busline83.mat','bus','line');

if event ==2 % generator trip
    t01=51;t02=151; 
    total=16;
elseif (event==3)%short circuit
     t01=71;t02=171; 
     total=68;
elseif (event==4)% load change
     t01=101;t02=201; 
     total=52;
else % event==1: linetrip
     t01=51;t02=151;
     total=58;
end
 
for p=1:total
    if event ==2 
       [bus_v,bus_freq,name ]=generdata(p, current_path) ; 
       [~,strname,~]=fileparts(name);
       GroundTruth(p)=str2num(strname(8:end));
       localbus(p,1:kmax)=location_rate(bus_v,t01,kmax,thres); 
       FromBus=line(GroundTruth(p),1);
       ToBus=line(GroundTruth(p),2);
       rank1=find(localbus(p,1:kmax)==FromBus);
       rank2=find(localbus(p,1:kmax)==ToBus);
    elseif (event==3)
       [bus_v,bus_freq ,name]= shortdata(p,pwd); %t+30
       [~,strname,~]=fileparts(name);
       GroundTruth(p)=str2num(strname(8:end));
       localbus(p,1:kmax)=location_rate(bus_v,t01,kmax,thres); 
       rank1=find(localbus(p,1:kmax)==GroundTruth(p));
       rank2=find(localbus(p,1:kmax)==GroundTruth(p));
    elseif (event==4)
       [bus_v,bus_freq,name ]=loaddata(p, current_path) ;% t+10
       [~,strname,~]=fileparts(name);
       GroundTruth(p)=str2num(strname(8:end));
       localbus(p,1:kmax)=location_rate(bus_v,t01,kmax,thres); 
       rank1=find(localbus(p,1:kmax)==GroundTruth(p));
       rank2=find(localbus(p,1:kmax)==GroundTruth(p));
    else
       [bus_v,bus_freq ,name]= linedata(p,pwd); %t+30
       [~,strname,~]=fileparts(name);
       GroundTruth(p)=str2num(strname(8:end));
       localbus(p,1:kmax)=location_rate(bus_v,t01,kmax,thres); 
       FromBus=line(GroundTruth(p),1);
       ToBus=line(GroundTruth(p),2);
       rank1=find(localbus(p,1:kmax)==FromBus);
       rank2=find(localbus(p,1:kmax)==ToBus);
    end
   
   if isempty(rank1) || isempty(rank2)
       fprintf('Not in the top %d buses \n',kmax);
       continue;
   end
   rank(p)=min(rank1,rank2);
end

AverRank=mean(rank(find(rank~=0)))

% true_num=0; 
%  
% for i=1:total
%     
%     % if localbus(i)==i+52
% % if localbus(i)==i
% % if localbus(i,1)==line(i,1)|| localbus(i,1)==line(i,2) || localbus(i,2)==line(i,1)|| localbus(i,2)==line(i,2)|| localbus(i,3)==line(i,1)|| localbus(i,3)==line(i,2)
% % if localbus(i,1)==line(i,2)|| localbus(i,2)==line(i,2)|| localbus(i,3)==line(i,2)
% 
%     if localbus(i,1)==i  || localbus(i,2)==i|| localbus(i,3)==i || localbus(i,4)==i|| localbus(i,5)==i
%         true_num=true_num+1;
%     end
% end
%  
% 
% 
% 
%  rate=true_num/total


