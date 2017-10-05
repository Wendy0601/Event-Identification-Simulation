function [k1,k2,min_angle, index_small,E1,E2,allnames]=freq_compare_subangle(root,root_dic,thres,gap,event,dic,noise,noise_level) 

% load('Dic3new.mat','Dic');%  dic_line=Dic{1,dic}; 
% load('H:\identification_data\Subspace_Angle\Dictionary.mat','Dic_new'); Dic=Dic_new;
if event ==1  
    path=fullfile(root,'Line_Trip_classic4_86');
elseif event==2
    path=fullfile(root,'Generator_trip_classic_86');
elseif event==3
    path=fullfile(root,'Three_Phase_classic_86'); 
elseif event==4
    path=fullfile(root,'Load_change_classic_86');
elseif event==5
    path=fullfile(root,'Line_to_ground_classic_86');
%     path=fullfile(root,'Capacitor_Bank_classic_86_new');
elseif event==6
    path=fullfile(root,'Motor_Off_classic_86');     
end
% 
% rootdic='H:\identification_data\Train_NewDic';
% if dic==1  
%     pathdic=fullfile(rootdic,'Line_Trip_classic4_86');
% elseif dic==2
%      pathdic=fullfile(rootdic,'Generator_trip_classic_86');
% elseif dic==3
%      pathdic=fullfile(rootdic,'Three_Phase_classic_86'); 
% elseif dic==4
%      pathdic=fullfile(rootdic,'Load_change_classic_86');
% end
if dic ==1 
    pathdic=fullfile(root_dic, 'Line_trip_dic');    
elseif dic==2
    pathdic=fullfile(root_dic, 'Generator_trip_dic');    
elseif dic==3
    pathdic=fullfile(root_dic, 'Three_Phase_dic');     
elseif dic==4
    pathdic=fullfile(root_dic, 'Load_change_dic');    
elseif dic==5
    pathdic=fullfile(root_dic, 'Line_to_ground_dic');    %Capacitor_Bank_classic_86';
elseif dic==6 
    pathdic=fullfile(root_dic, 'Motor_off_dic');      
end
 
[ allnames, len ] = read_all_file( path,pwd );
[ allnamesDic, lendic ] = read_all_file( pathdic,pwd );
Vk={};
for i=1:len   % p is the number of test data
    for k=1:lendic
%         j= dic_line(k);
        [bus_v1,bus_freq1 ,name1]=readPMUdata(path,allnames,i); 
        if event ==2
           t01=53;t02=153;
           if i<16
               X1=[abs(bus_freq1(1:51+i,:)) ; abs(bus_freq1(min(53+i,68):68,:))]; 
           else
               X1=[abs(bus_freq1(1:67,:))];
           end
        elseif (event==3 | event==5 )
           t01=73;t02=173; 
           X1=abs(bus_freq1(1:68,:));
        elseif (event==4  )
           t01=103;t02=203; 
           X1=abs(bus_freq1(1:68,:));
        else
           t01=53;t02=153; 
           X1=abs(bus_freq1(1:68,:));
        end
        %dataset 2
        [bus_v2,bus_freq2 ,name2]=readPMUdata(pathdic,allnamesDic,k); 
        if dic ==2
           t01=53;t02=153;
           X2=[abs(bus_freq2) ];
        elseif (dic==3 ||dic==5 )
           t01=73;t02=173; 
           X2=abs(bus_freq2(1:68,:));
        elseif (dic==4  )
           t01=103;t02=203; 
           X2=abs(bus_freq2(1:68,:));
        else
           t01=53;t02=153;
           X2=abs(bus_freq2(1:68,:));
        end
        %compute subspace angle
        X1=sub_rowmean(X1,t01); 
        [U1,S1,V1] = svd(X1(:,t01:3:t02) );%%%%%%%%%%%%%%%%%%%%% here the sampling time is 33/seond 
        s1=diag(S1);
        k1(i)=choose_rank(s1,thres,gap);   
        E1(i)=sum(s1(1:k1(i)));
%         figure;plot(V2(:,1:k2(k)));
%         figure; plot(s2,'o'); k2(k)

         X2=sub_rowmean(X2,t01); 
        [U2,S2,V2] = svd(X2(:,t01:3:t02) );%%%%%%%%%%%%%%%%%%%%% here the sampling time is 33/seond 
        s2=diag(S2);
        k2(k)=choose_rank(s2,thres,gap);   
        E2(k)=sum(s2(1:k2(k)));
 
        k12=max(k1(i),k2(k));
        SubspaceAngles(i,k)=angle0(V1(:,1:k12),V2(:,1:k12)); 
        Vk{1,k}=V2(:,1:k2(k));
    end
    smallest=min(SubspaceAngles(i,:)) ;
    min_angle(i)=smallest(1);
    index_s=find(SubspaceAngles(i,:)==min_angle(i));
    index_small(i)=index_s(1);
    Vk_dic{1,i}=Vk{1,index_small};
    Vk_event{1,i}=V1(:,1:k1(i));
end
% if event ==1 
%     save H:\identification_data\Subspace_Angle\angle_Newdic_Line_Trip_classic4_86.mat
% elseif event==2
%     save H:\identification_data\Subspace_Angle\angle_Newdic_Generator_trip_classic_86.mat
% elseif event==3
%     save H:\identification_data\Subspace_Angle\angle_Newdic_Three_Phase_classic_86.mat
% elseif event==4
%     save H:\identification_data\Subspace_Angle\angle_Newdic_Load_change_classic_86.mat
% end 
end 



% figure;plot(X1(28:32,:)');
% [U1,S1,V1] = svd(X1(:,t01:3:t01+89) );%%%%%%%%%%%%%%%%%%%%% here the sampling time is 33/seond 
% s1=diag(S1);
% xlabel('Time (0.01 second)')
% ylabel ('Voltage Magnitudes (p.u.)');
% figure; plot(s1,'o')
% title('30 singular values of voltage magnitudes under 100dB noise') 
