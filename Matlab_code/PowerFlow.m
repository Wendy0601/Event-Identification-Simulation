clear all; close all; clc 
data16em2
tol = 1e-8;   % tolerance for convergence
iter_max = 30; % maximum number of iterations
acc = 1.0;   % acceleration factor
[bus_sol,line,line_flw] = ...
         loadflow(bus,line,tol,iter_max,acc,'n',2);

FromBus = line(:,1);
ToBus = line(:,2);
for i=1:size(line,1)
     index_frombus=find(bus_sol(:,1)==FromBus(i));
     index_tobus=find(bus_sol(:,1)==ToBus(i));
    V1(i)=bus_sol(index_frombus,2);
    V1ang(i) = bus_sol(index_frombus,3)*pi/180;
    V2(i)=bus_sol(index_tobus,2);
    V2ang(i) = bus_sol(index_tobus,3)*pi/180;

end
V1 = V1.*exp(j*V1ang); V1=V1';
V2 = V2.*exp(j*V2ang); V2=V2';
R = line(:,3);
X = line(:,4);
B = line(:,5);
tap = line(:,6);
phi = line(:,7);
[S1,S2] = line_pq(V1,V2,R,X,B,tap,phi);
% subplot(2,1,1);plot(t,real(S1));title('active power /pu')
% subplot(2,1,2);plot(t,imag(S1));title('reactive power /pu')
P=real(S1);
Q=imag(S1);
P1=P;
Q1=Q;
% 
% 
% 
% % load('chl6_2dif.mat');
% % load('H:\Old data\bigbus.mat');% 28
% load('no1267outline10-11');% no5biglines12.mat no561013.mat no5lines.mat no 5 lines
% % % load('bus59.mat');% 58
% busnum1 = bus_sol(:,1);
% FromBus1 = line(:,1);
% ToBus1 = line(:,2);
% % From_idx = com_index(busnum,FromBus);
% % To_idx = com_index(busnum,ToBus);
% for i=1:size(line,1)
%     V1(i)=bus_sol(FromBus1(i),2);
%     V1ang = bus_sol(FromBus1(i),3)*pi/180;
%     V2(i)=bus_sol(ToBus1(i),2);
%     V2ang = bus_sol(ToBus1(i),3)*pi/180;
% 
% end
%  
% V1 = V1.*exp(j*V1ang);
%  
% V2 = V2.*exp(j*V2ang);
% R1 = line(:,3);
% X1 = line(:,4);
% B1 = line(:,5);
% tap1 = line(:,6);
% phi1 = line(:,7);
% [S11,S22] = line_pq(V1,V2,R1,X1,B1,tap1,phi1);
% % subplot(2,1,1);plot(t,real(S1));title('active power /pu')
% % subplot(2,1,2);plot(t,imag(S1));title('reactive power /pu')
% P=real(S11);
% Q=imag(S11);
% P2=P;
% Q2=Q;
% 
% %%%%%%%%%%% the change of initial condition
% 
% %%%%%% power flow difference ratio and its plot
% for i=1:83
%     meanp(i)= abs((P1(i,10)-P2(i,10))/abs(P1(i,10)));
% end
% meanp1=mean(meanp)
% plot(P1(1:4,10:20)','-.*','Linewidth',2)
% figure; plot(P2(5:8,10:20)','--o','Linewidth',2)