%% Figure 8 in paper 
clear all;clc;close all; 
load('H:\identification_data\outline\olddata\outline1.mat')
thres=0.99; t01=51;t02=151;X=abs(bus_v(1:68,:));
X=abs(bus_v(1:68,:));
[row,col]=size(X);
          for i=1:row
                X(i,:)=X(i,:)-mean(X(i,1:t01-2))*ones(1,col); 
          end
[Vk01,Vk02,Vk03,Vk04,Vk05,Vk06,Vk07,Vk08,k1]=rowsub(X,t01,t02,thres); 
X1=Vk06/norm(Vk06);
X1=X1';
 figure; plot(X1(1,1:34)','r.-','Marker','o');
    hold on;plot(X1(2,1:34)','k--','Linewidth',2,'Marker','.');
    hold on;plot(X1(3,1:34)','b.-','Linewidth',2,'Marker','square');
    hold on;plot(X1(4,1:34)','c.-','Linewidth',2,'Marker','+');
    hold on;plot(X1(5,1:34)','g.-','Linewidth',2,'Marker','<');
    hold on;plot(X1(6,1:34)','m.-','Linewidth',2,'Marker','x');
     xlabel('Time (second)','FontSize',15)
          ylabel('Vr','FontSize',15);
           set(gca,'xtick',1:11:34)
         set(gca,'xticklabel',{'0.5','0.8','1.1','1.5'})
%% Figure 8 in paper 
clear all;clc;close all;          
load('H:\identification_data\outline\initialline59degree\outline1.mat', 'bus_v')
thres=0.99; t01=51;t02=151;X=abs(bus_v(1:68,:));
X=abs(bus_v(1:68,:));
[row,col]=size(X);
          for i=1:row
                X(i,:)=X(i,:)-mean(X(i,1:t01-2))*ones(1,col); 
          end
[Vk01,Vk02,Vk03,Vk04,Vk05,Vk06,Vk07,Vk08,k1]=rowsub(X,t01,t02,thres); 
X1=Vk06/norm(Vk06);
X1=X1';
% X1=[zeros(4,17) X1];
 figure; plot(X1(1,1:34)','r.-','Marker','o');
    hold on;plot(X1(2,1:34)','k--','Linewidth',2,'Marker','.');
    hold on;plot(X1(3,1:34)','b.-','Linewidth',2,'Marker','square');
    hold on;plot(X1(4,1:34)','c.-','Linewidth',2,'Marker','+');
    hold on;plot(X1(5,1:34)','g.-','Linewidth',2,'Marker','<');
    hold on;plot(X1(6,1:34)','m.-','Linewidth',2,'Marker','x');
    xlabel('Time (second)','FontSize',20)
          ylabel('Vr','FontSize',20);
         set(gca,'xtick',1:11:34)
         set(gca,'xticklabel',{'0.5','0.8','1.1','1.5'})
%% Figure 7 in paper 
clear all;clc;close all; 
load('H:\identification_data\outline\olddata\outline1.mat');
index=[1 25 40 48]; X1=abs(bus_v(index,:)); 
X2=X1(:,1:3:500); ax=0.03*[1:167];
    figure;plot(ax,X2(1,1:3:end),'b.-','Linewidth',2,'Marker','square');
    hold on;plot(ax,X2(2,:),'m.-','Linewidth',2,'Marker','x');
    hold on;plot(ax,X2(3,:),'c.-','Linewidth',2,'Marker','+');
    hold on;plot(ax,X2(4,:),'g.-','Linewidth',2,'Marker','<');
    xlabel('Time (second)','FontSize',20)
          ylabel('Voltage Magnitude (pu)','FontSize',20);
%% Figure 7 in paper 
clear all;clc;close all;          
load('H:\identification_data\outline\initialline59degree\outline1.mat', 'bus_v')
index=[1 25 40 48]; X1=abs(bus_v(index,:)); 
X2=X1(:,1:3:500); ax=0.03*[1:167];
    figure;plot(ax,X2(1,:),'b.-','Linewidth',2,'Marker','square');
    hold on;plot(ax,X2(2,:),'m.-','Linewidth',2,'Marker','x');
    hold on;plot(ax,X2(3,:),'c.-','Linewidth',2,'Marker','+');
    hold on;plot(ax,X2(4,:),'g.-','Linewidth',2,'Marker','<');
    xlabel('Time (second)','FontSize',15)
          ylabel('Voltage Magnitude (pu)','FontSize',15);
          
          
load('H:\Labfile_realdata\ISO_New_England_Files\Pump Storage (Load) Tripping Events\NE_Northfield_Mountain_Trip_8_200MW.mat', 'V_pm')          
v=abs(V_pm(:,55:56));figure;plot(v(:,2),'r','Linewidth',2);hold on;plot(v(:,1),'b.-','Linewidth',2 )
xlabel('Time (0.033 second)')
ylabel('Voltage magnitude (V)')

%load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2015-11-27 Line 368 Trip and Close\2015_11_27_Line_368_Trip.mat')
% v=abs(V_pm(:,41:42));figure;plot(v(:,:),'r','Linewidth',2);hold on;plot(v(:,1),'b.-','Linewidth',2 )

load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2016-07-07 Line T172 Fault and S171 Tripped\2016_07_07_S171_Trip.mat', 'V_pm')
v=abs(V_pm(:,:));figure;plot(v(:,82),'r','Linewidth',2);hold on;plot(v(:,85),'b.-','Linewidth',2 )
 
load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2015-10-17 Line 354 Trip\2015_10_17_Line354_Trip.mat', 'V_pm')
v=abs(V_pm(:,55:56));figure;plot(v(:,2),'r','Linewidth',2);hold on;plot(v(:,1),'b.-','Linewidth',2 )
xlabel('Time (0.033 second)')
ylabel('Voltage magnitude (V)')


deltat=0; 
load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping\Line Tripping\2015-10-17 Line 328 Tripped and Reclosed 3 times and outage\2015_10_17_Line_328_2nd_Trip_and_Reclose.mat', 'V_pm')
t=450+deltat;
v=abs(V_pm(:,81:82));figure;plot(v(:,2),'r','Linewidth',2);hold on;plot(v(:,1),'b.-','Linewidth',2 )
xlabel('Time (0.033 second)')
ylabel('Voltage magnitude (V)')
 
% 
% deltat=0; 
% load('H:\Labfile_realdata\ISONE-disturbance\Line Tripping Events\2016-10-11 240-601 and 336 Trip\Line_240-601_and_336_Trip.mat', 'V_pm')
% t=300+deltat;  % the second line trip
% v=abs(V_pm(:,:));figure;plot(v(:,2),'r','Linewidth',2);hold on;plot(v(:,1),'b.-','Linewidth',2 )
% xlabel('Time (0.033 second)')
% ylabel('Voltage magnitude (V)')


% figure in the response letter
load('H:\identification_data\Line_Trip_classic4_86\linetrip_4-14.mat', 'bus_v')
plot(abs(bus_v(4,:)'),'b.-','Marker','.')
hold on;plot(abs(bus_v(14,:)'),'k--','Linewidth',2,'Marker','.')
hold on;plot(abs(bus_v(5,:)'),'c.-','Linewidth',2,'Marker','+')
hold on;plot(abs(bus_v(6,:)'),'m.-','Linewidth',2,'Marker','x');
xlabel('Time (0.01 second)','FontSize',15)
ylabel('Voltage (p.u)','FontSize',15);

load('H:\identification_data\Capacitor_Bank_classic_86_new\Capacitor_Bank_4.mat', 'bus_v')
plot(abs(bus_v(4,:)'),'b.-','Marker','.')
hold on;plot(abs(bus_v(14,:)'),'k--','Linewidth',2,'Marker','.')
hold on;plot(abs(bus_v(5,:)'),'c.-','Linewidth',2,'Marker','+')
hold on;plot(abs(bus_v(6,:)'),'m.-','Linewidth',2,'Marker','x');
xlabel('Time (0.01 second)','FontSize',15)
ylabel('Voltage (p.u)','FontSize',15);

% subspace 
% thres=0.99;gap=10;
% load('H:\identification_data\Line_Trip_classic4_86\linetrip_4-14.mat', 'bus_v')
% t01=51;t02=151;
% X1=abs(bus_v(1:68,:));
% [U1,S1,V1] = svd(X1(:,t01:3:t02) );%%%%%%%%%%%%%%%%%%%%% here the sampling time is 33/seond 
% s1=diag(S1);
% k1=choose_rank(s1,thres,gap);   
% E1=sum(s1(1:k1));
% 
% load('H:\identification_data\Capacitor_Bank_classic_86_new\Capacitor_Bank_4.mat', 'bus_v');
% t01=51;t02=151;
% X2=abs(bus_v(1:68,:));
% X2=sub_rowmean(X2,t01);
% [U2,S2,V2] = svd(X2(:,t01:3:t02) );%%%%%%%%%%%%%%%%%%%%% here the sampling time is 33/seond 
% s2=diag(S2);
% k2=choose_rank(s2,thres,gap);   
% E2=sum(s2(1:k2));
%         
% k12=max(k1,k2);
% SubspaceAngles=angle0(V1(:,1:k12),V2(:,1:k12))


% 10dB-40dB
% load('H:\identification_data\condition_b\Line_Trip_classic4_86\outline2.mat', 'bus_v')
% noise=1;
% t01=51;t02=151;
% n=1;
% for noise_level=10:10:100 
%     X=abs(bus_v(1:68,:));
%     X1=add_noise(X,noise,noise_level ); 
% %     figure;plot(X1') 
%     [U1,S1,V1] = svd(X1);
%     Xr=U1(:,1:6)* S1(1:6,1:6)*V1(:,1:6)';
%     err(n)=100*norm((X1-Xr), 'fro')/norm(X1,'fro');
%     n=n+1;
% end
