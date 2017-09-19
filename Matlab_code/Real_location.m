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