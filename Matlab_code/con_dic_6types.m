clear all ; close all; clc; 
% load('H:\identification_data\Subspace_Angle\angle_Line_Trip_classic4_86.mat','SubspaceAngles');deg=5;event=1;
% load('H:\identification_data\Subspace_Angle\angle_Generator_trip_classic_86.mat','SubspaceAngles');deg=5;event=2;
% load('H:\identification_data\Subspace_Angle\angle_Three_Phase_classic_86.mat','SubspaceAngles');deg=10;event=3;
% load('H:\identification_data\Subspace_Angle\angle_Line_to_Line_classic_86.mat','SubspaceAngles');deg=10;event=4;
% load('H:\identification_data\Subspace_Angle\angle_Load_change_classic_86.mat','SubspaceAngles');deg=5;event=5;
% load('H:\identification_data\Subspace_Angle\angle_Capacitor_Bank_classic_86.mat','SubspaceAngles');deg=5;event=6;
load('H:\identification_data\Subspace_Angle\angle_Motor_Off_classic_86.mat','SubspaceAngles');deg=5;event=7;

 %% find out how many sub_angle are less then deg
n=size(SubspaceAngles,1); 
for i=1:n
    small_number=find(SubspaceAngles(i,:)<=deg);
    num_small(i)=numel(small_number);
end

%% find out the dictionary if choose those who has more than j small subangle. Finally the best j is selected 
olddic=zeros(1,n);
 for j= 2: max(num_small)-1 
    dic=[];newdic=[];dicangle=[];dic_angle=[];%initalize 
    dic=find(num_small>j);
    m=numel(dic);
    for k=1:n
        for p=1:m
            dicangle(p)=SubspaceAngles(k,dic(p));
        end
            dic_angle(k)=min(dicangle);
    end
    
    adddic=find(dic_angle>deg);
    newdic=[dic adddic];
    if numel(olddic) >= numel( newdic)
        olddic=newdic;% the dictioanry is updated so that the minimum number is created
    else
        
    olddic=olddic;
    end
 end
 
 dic=[];newdic=[];dicangle=[];dic_angle=[];
 m=numel(olddic);
  for k=1:n
        for p=1:m
            dicangle(p)=SubspaceAngles(k,olddic(p));
        end
            dic_angle(k)=min(dicangle);% with the new dictionary, all the sub_angle
    end
 
 newdic=olddic;
 path=pwd;cd('H:\identification_data\Subspace_Angle');load('Dictionary.mat','Dic');cd(path)     
 Dic{1,event}=newdic;   
path=pwd;cd('H:\identification_data\Subspace_Angle');save  Dictionary.mat;cd(path) 