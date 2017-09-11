clear all; close all; clc 
%run the s_simu and produce the voltages V1 V2 of all buses and save the dataset as the following .mat format
% load('4-14topology_change.mat');
% load('8-9topology_change.mat');
% load('33-38topology_change.mat');
% load('3-18topology_change.mat');
load('2-3and8-9topology_change.mat');
% load('NoTopology_change.mat');
 
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
 
% these are the power flow with time, you can choose the power flow after 10 or more steps when the initialization is done.
