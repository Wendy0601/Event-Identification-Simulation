clear all; close all; clc 
data16em2
% do any kind of events and run the non linear simulation to produce V1,V2, and bus_sol
bus(:,6)=bus(:,6)+0.05; % this is power injection change
[bus_v,bus_freq,theta,mac_ang,mac_spd,flag_error,ilf,V1,V2,bus_sol]=s_simulwt(sw_con,load_con,lmod_con,bus,line);
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
