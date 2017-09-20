<<<<<<< HEAD
clear all; clc; close all;
%parameters
delay=0;%% time delay test
thres=0.97; 
gap=10;
T=30;
deltat=7;
kmax=1;
% read datasets
event=1;event_num=9;
% event=2;event_num=13;
% event=3;event_num=11;
for data=1:event_num
     if event ==1  
        [V_pm1,t]=Real_linetrip(data,0);
    elseif event==2
        [V_pm1,t]=Real_loadchange(data,0);
    else 
        [V_pm1,t]=Real_faults(data,deltat);
    end   
% data process 
    V_nan=delete_nan(V_pm1);% this function is used to delete the nan vectors
    V_unit=unit_colum(V_nan,t);% this is used to change it into unit value
    V_good= delete_bad( V_unit,t-50,0.1);%% delete bad data  
    V003=dele_nan_normalize(V_good,t);% this is used to delecte nan it further
%% smooth
    Vo=V003(:,:)';
    Ve=expsmooth( Vo  ,30,0.7); % figure;subplot(1,2,1);plot( Ve);subplot(1,2,2); plot(V003(:,:)')
    U_pm01=Ve(t:t+T/3,:)';
    [U,S,V] =svd(U_pm01);
    [row, col]=size(U_pm01);
    s = diag(S);   
    sum_s=0;
    s_k=0;
    while sum_s<thres*sum(s)
        sum_s=sum_s+s(s_k+1);
        s_k=s_k+1;
    end

    ku=s_k; 
    Ukk=U(:,1:ku)*S(1:ku,1:ku);%*U(:,1:ku)';
    wei=zeros(row,1);
    for i=1:row
        wei(i)=norm((Ukk(i,:)));
    end
    wei=wei/sum(wei);
    localbus(data,1:kmax)=find_k_max(wei,kmax);
end
=======
clear all; clc; close all;
%parameters
delay=0;%% time delay test
thres=0.97; 
gap=10;
T=30;
deltat=0;
kmax=2;
% read datasets
% event=1;event_num=9;
% event=2;event_num=13;
event=3;event_num=11;
for data=1:event_num
     if event ==1  
        [V_pm1,t]=Real_linetrip(data,0);
    elseif event==2
        [V_pm1,t]=Real_loadchange(data,0);
    else 
        [V_pm1,t]=Real_faults(data,deltat);
     end   
    
% data process 
    [V_nan, index_nan ]=delete_nan(V_pm1);% this function is used to delete the nan vectors
    % remember the original index of bus and delete those nan buses: sort
    % the index of nan buses and delete them one by one, notice that after
    % delete one, the remaining should be reduce one as index_nan(iBus+1)=index_nan(iBus+1)-1;.
    index_bus=(1:size(V_pm1,2));
    index_nan=sort(unique(index_nan(:)))';
     for iBus=1:numel(index_nan) 
         index_bus(index_nan(iBus ))=[];
         if iBus<numel(index_nan)
              index_nan(iBus+1)=index_nan(iBus+1)-1;
         end
     end
    V_unit=unit_colum(V_nan,t);% this is used to change it into unit value 
    [V_unit,index_bad]=dele_nan_normalize(V_unit,t);% this is used to delete nan it further 
    index_bad=sort(unique(index_bad(:)))'; 
    for iBus=1:numel(index_bad) 
         index_bus(index_bad(iBus))=[];
         if iBus<numel(index_bad)
              index_bad(iBus+1:end)=index_bad(iBus+1:end)-1;
         end
     end
% bus_index of those not nan 
    U_pm01=V_unit(:,t:t+T/3,:);
    [U,S,V] =svd(U_pm01);
    [row, col]=size(U_pm01);
    s = diag(S);   
    sum_s=0;
    s_k=0;
    while sum_s<thres*sum(s)
        sum_s=sum_s+s(s_k+1);
        s_k=s_k+1;
    end

    ku=s_k; 
    Ukk=U(:,1:ku)*S(1:ku,1:ku);%*U(:,1:ku)';
    wei=zeros(row,1);
    for i=1:row
        wei(i)=norm((Ukk(i,:)));
    end
    wei=wei/sum(wei);
    localbus(data,1:kmax)=find_k_max(wei,kmax);
    OriginLocation(data,1:kmax)=index_bus(localbus(data,1:kmax)); % this is the original index of buses
end

load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2015-10-17 Line 328 Tripped and Reclosed 3 times and outage\2015_10_17_Line_328_1st_Trip_and_Reclose.mat','VoltageSignal');
Bus_info=VoltageSignal;
for i=1:numel(OriginLocation(:,1))
result_bus(i,1)=VoltageSignal(OriginLocation(i,1));
end



result_bus(:)


if event==1
    GT={'Sherman11807',81,82,83,84;'Sherman11807',81,82,83,84;'Sherman11807',81,82,83,84;'Sherman11807',81,82,83,84;'northfield_mountain_16r_12218',55,56,[],[];'v_manchester_3a_12232_cx_manchestr_s2_b345_bbus_vp_30',41,42,[],[];'v_west_medway_446_12006_be_wmedway_2_b345_bbus_vs_30',114,115,[],[];'v_wachusett_11808_ne_wachusett_a_l345_308_vp_30',106,107,108,[];'v_keene_road_11051_bh_keeneroad_l345_3015_vs_30',27,228,29,30};
elseif event==2
    GT={'v_northfield_mountain_16r_12218_cx_nfldmt16r_l345_312_vp_30',55,56,[];'v_northfield_mountain_16r_12218_cx_nfldmt16r_l345_312_vp_30',55,56,[];'v_northfield_mountain_16r_12218_cx_nfldmt16r_l345_312_vp_30',55,56,[];'v_northfield_mountain_16r_12218_cx_nfldmt16r_l345_312_vp_30',55,56,[];'v_sandy_pond_11806_ne_sandypond_c_l115_138e_vp_30',72,73,74;'v_sandy_pond_11806_ne_sandypond_c_l115_138e_vp_30',72,73,74;'v_northfield_mountain_16r_12218_cx_nfldmt16r_l345_312_vp_30',55,56,[];'v_northfield_mountain_16r_12218_cx_nfldmt16r_l345_312_vp_30',55,56,[];'v_northfield_mountain_16r_12218_cx_nfldmt16r_l345_312_vp_30',55,56,[];'v_northfield_mountain_16r_12218_cx_nfldmt16r_l345_312_vp_30',55,56,[];'v_sandy_pond_11806_ne_sandypond_c_l115_138e_vp_30',72,73,74;'v_sandy_pond_11806_ne_sandypond_c_l115_138e_vp_30',72,73,74};
elseif event==3
    GT={'v_lake_road_27e_12211_cx_lakerd27e_b345_abus_vp_30',31,32,[],[],[],[];'v_ward_hill_11811_ne_wardhill_b_b345_b2_vp_30',111,113,[],[],[],[];'v_wakefield_junction_11812_ne_wakefldjc_a_b345_b2_vp_30',109,110,[],[],[],[];'v_millbury_11803_ne_millbury_3b_l345_366_vp_30',43,44,45,46,47,48;'v_wakefield_junction_11812_ne_wakefldjc_a_b345_b2_vp_30',109,110,[],[],[],[];'v_millbury_11803_ne_millbury_3b_l345_366_vp_30',43,44,45,46,47,48;'v_southington_4c_12222_cx_sothgtn4c_1_b345_a2bus_vp_30',89,90,[],[],[],[];'v_sherman_road_11807_ne_shermanrd_l345_347_vp_30',81,82,83,84,[],[];'v_sherman_road_11807_ne_shermanrd_l345_347_vp_30',81,82,83,84,[],[];'v_west_medway_446_12005_be_wmedway_1_b345_abus_vs_30',114,115,[],[],[],[];'v_millstone_15g_12217_cx_milstn15g_b345_abus_vp_30',49,50,[],[],[],[]};
end
>>>>>>> 45f284448641ac5cb0b549514919c3638a71950a
