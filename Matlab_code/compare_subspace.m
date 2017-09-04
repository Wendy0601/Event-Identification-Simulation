function [k1,k2,min_angle, index_small,E1,E2,allnames]=compare_subangle(root,thres,gap,event,dic,noise,noise_level) 
 
% load('H:\identification_data\Subspace_Angle\Dictionary_gap0.mat','Dic_gap0');%  Dic
load('Dic3new.mat','Dic');%  
dic_line=Dic{1,dic}; 
if event ==1  
    path=fullfile(root,'Line_Trip_classic4_86');
elseif event==2
    path=fullfile(root,'Generator_trip_classic_86');
elseif event==3
    path=fullfile(root,'Three_Phase_classic_86'); 
elseif event==4
    path=fullfile(root,'Load_change_classic_86');
elseif event==5
    path=fullfile(root,'Capacitor_Bank_classic_86_new');
elseif event==6
    path=fullfile(root,'Motor_Off_classic_86');     
end

if dic ==1 
    pathdic='H:\identification_data\Line_Trip_dic';%'H:\identification_data\Line_trip_86';    
elseif dic==2
    pathdic='H:\identification_data\Generator_trip_dic';
elseif dic==3
    pathdic='H:\identification_data\Three_Phase_dic'; 
elseif dic==4
    pathdic='H:\identification_data\Load_change_dic';
elseif dic==5
    pathdic='H:\identification_data\Capacitor_Bank_dic';%Capacitor_Bank_classic_86';
elseif dic==6 
    pathdic='H:\identification_data\Motor_off_dic';   
end
 
[ allnames, len ] = read_all_file( path,pwd );
[ allnamesDic, ~ ] = read_all_file( pathdic,pwd );
Vk={};
for i=1:len  % p is the number of test data
    for k=1:numel(dic_line)
        j= dic_line(k);
        [bus_v1]=readPMUdata(path,allnames,i); 
        if event ==2
           t01=51;t02=151;
           X1=[abs(bus_v1(1:51+i,:)) ; abs(bus_v1(min(53+i,68):68,:))];
        elseif (event==3)
           t01=71;t02=171; 
           X1=abs(bus_v1(1:68,:));
        elseif (event==4)
            t01=101;t02=201; 
           X1=abs(bus_v1(1:68,:));
        else
           t01=51;t02=151;
           X1=abs(bus_v1(1:68,:));
        end
        %dictionary
        [bus_v2]=readPMUdata(pathdic,allnamesDic,j); 
        if dic ==2
           t01=51;t02=151;
           X2=[abs(bus_v2(1:51+j,:)) ; abs(bus_v2(min(53+j,68):68,:))];
        elseif (dic==3 )
           t01=71;t02=171; 
           X2=abs(bus_v2(1:68,:));
        elseif (event==4)
            t01=101;t02=201; 
           X2=abs(bus_v2(1:68,:));
        else
           t01=51;t02=151;
           X2=abs(bus_v2(1:68,:));
        end
        %compute subspace angle
        X1=add_noise( X1,noise,noise_level );
        X1=sub_rowmean(X1,t01);
        
        [U1,S1,V1] = svd(X1(:,t01:3:t02) );%%%%%%%%%%%%%%%%%%%%% here the sampling time is 33/seond 
        s1=diag(S1);
        k1(i)=choose_rank(s1,thres,gap);   
        E1(i)=sum(s1(1:k1(i)));

%         X2=add_noise( X2,noise,noise_level ); 
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
    
end 
