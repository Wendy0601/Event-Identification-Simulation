clear    ;close all;clc;
%% parameter
thres=0.97; % the more serious the fault is, then this threshold should be set 
load('Dic_generline.mat','Dic_real');  %Dic_smooth.mat dic_smooth_gener.mat Dic_real_few.mat Dic_real_few.mat Dic_Vm.mat  Dic_real_unit.mat Dic_real.mat Dic_real{1,j};j=1 faultline; j=2 load shedding; j=3 generator trip 
delay=0;
gap=10;
 
% %% test by three dictionary types
% for dic_number=1 :3
% dic=Dic_real{1,dic_number};%%%%%%%%%%%% choose a dictionary atom 1:outline ;2--load unit trip;3--gener trip
% if dic_number==1  
%           
%         Vk_dic1=cell(1,numel(dic)) ;
%         for k=1:numel(dic)
%               j= dic(k);
%               Vk_dic1{1,k}=data_line_trip_few(j,thres,delay);% construct dictionary
%         end  
%              
% end
% 
% if dic_number==2   
%           
%         Vk_dic2=cell(1,numel(dic)) ;
%         for k=1:numel(dic)
%               j= dic(k);  
%               Vk_dic2{1,k}=data_unit_trip_few(j,thres,delay);% construct dictionary
%         end
%              
% end
% 
% if dic_number==3    
%     
%     Vk_dic3=cell(1,numel(dic)) ;
%      
%         for k=1:numel(dic)
%               j= dic(k);
%               Vk_dic3{1,k}=data_gener_trip_few(j,thres,delay);% construct dictionary
%         end
%              
% end
% end


% load('original_dic.mat','Vk_dic1','Vk_dic2','Vk_dic3');
load('original_smooth_dic.mat','Vk_dic1','Vk_dic2','Vk_dic3');

% load ('dic_delay_unit.mat','Vk_dic1','Vk_dic2','Vk_dic3');
% load('load5gap.mat','Vk_load');
 %% test gener
%  totalcase=12; 
%     for p=1:totalcase  
%          Vk=data_gener_trip_few(p,thres,delay);
%          
%          angles=[]; 
%         for ka=1:size(Vk_dic1,2)
%               angles(ka)=angle0(Vk,Vk_dic1{1,ka});
%         end
%         min_angle_gener(1,p)=min(angles);
%         
%         angles=[]; 
%         for ka=1:size(Vk_dic2,2)
%               angles(ka)=angle0(Vk,Vk_dic2{1,ka});
%         end
%         min_angle_gener(2,p)=min(angles);
%         
%         angles=[]; 
%         for ka=1:size(Vk_dic3,2)
%               angles(ka)=angle0(Vk,Vk_dic3{1,ka});
%         end
%         min_angle_gener(3,p)=min(angles);
%     end
%     min_angle_gener=abs(min_angle_gener);
    %% test load
%      totalcase=8; 
%     for p=1:totalcase  
%          Vk=data_unit_trip_few(p,thres,delay);
%          
% %        
%          angles=[]; 
%         for ka=1:size(Vk_dic1,2)
%               angles(ka)=angle0(Vk,Vk_dic1{1,ka});
%         end
%         whole_angle_unit1(p,:)=angles;
%         whole_angle_unit1=abs(whole_angle_unit1);
%         min_angle_unit(1,p)=min(angles);
%         
%         angles=[]; 
%         for ka=1:size(Vk_dic2,2)
%               angles(ka)=angle0(Vk,Vk_dic2{1,ka});
%         end
%         min_angle_unit(2,p)=min(angles);
%          whole_angle_unit2(p,:)=angles;
%          whole_angle_unit2=abs(whole_angle_unit2);
%          
%         angles=[]; 
%         for ka=1:size(Vk_dic3,2)
%             Vk_dic=Vk_dic3{1,ka};
%               angles(ka)=angle0(Vk,Vk_dic);
%         end
%         min_angle_unit(3,p)=min(angles);
%          whole_angle_unit3(p,:)=angles;
%          whole_angle_unit3=abs(whole_angle_unit3);
%     end
%     min_angle_unit=abs(min_angle_unit);
%     
      
    % test line trip
     totalcase=6; 
    for p=1:totalcase  
         Vk=data_line_trip_few(p,thres,delay);
         
         angles=[]; 
        for ka=1:size(Vk_dic1,2)
              angles(ka)=angle0(Vk,Vk_dic1{1,ka});
        end
        min_angle_line(1,p)= (min(angles));
        whole_angle_line1(p,:)=angles;
        whole_angle_line1=abs(whole_angle_line1);
        min_angle_line(1,p)=min(angles);
        
        angles=[]; 
        for ka=1:size(Vk_dic2,2)
              angles(ka)=angle0(Vk,Vk_dic2{1,ka});
        end
        min_angle_line(2,p)= (min(angles));
         whole_angle_line2(p,:)=angles;
        angles=[]; 
        for ka=1:size(Vk_dic3,2)
              angles(ka)=angle0(Vk,Vk_dic3{1,ka});
        end
        min_angle_line(3,p)=  (min(angles));
         whole_angle_line3(p,:)=angles;
    end
    min_angle_line=abs(min_angle_line);



%%
% V_LT4=Vk_dic1{1,2};figure;plot(V_LT4,'.-');ylabel('V_LT4')
%   figure;plot(Vk_dic2{1,2},'.-'); ylabel('Vk_dic')
%   figure;plot(Vk,'.-'); ylabel('Vk')
% V_G3=Vk_dic3{1,3}; figure;plot(V_G3,'.-');ylabel('V_G3')