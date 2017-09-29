clear ;clc;close all; 
%% this file uses Prony method to test the dominant frequency components, the freq_sec is the main frequency components and their magnitudes are in A; the paramter order can be changed to adjust the resolution and the time period is selected just afrer event 70 to 170;
load('freq_30s_linetrip.mat','bus_v','mac_spd','bus_freq');%line1-2_10s.mat mac_30second.mat
% parameters   
dt=0.03; 
% X1= (bus_freq(1,:)); 
X1=abs(bus_v(2,:));
for order=12% the initial order
p=order;
t1=51;
[row,col]=size(X1);
% normalize
for i=1:row
                X1(i,:)=X1(i,:)-mean(X1(i,1:t1-2))*ones(1,col);               
end
% use one second data
y=X1(1,70:3:170);
y=y';
L=length(y);
% parameters
n=p ;% order
N=L;
x1=-y(n+1:N);
for i=1:n
    x2(:,i)=y(n-i+1:N-i);
end
alfa=x2\x1;
alfa1=[1 ; alfa];
mu=roots(alfa1);
mu=mu';
for k=1:N
    T(k,:)=conj(mu).^(k-1);
end
F=y(1:N); 
C= inv(conj(T')*T)*conj(T')*F;
%% the estimated signal
f_bar=T*C;
%figure;plot(y','g ','lineWidth',2);hold on; plot(real(f_bar),'r.-'); xlabel('Time second');ylabel('Rotor speed pu')
Freq=imag(log(mu))/(2*pi*dt);% frequency
A=abs(C);% magnitude
figure;subplot(2,1,1);stem(A); 
subplot(2,1,2);stem(Freq)
end
%snr =20*log(norm( (y))/norm( (y)-f_bar)) 
%% further fit the signal
select=find(A>0.00002)
C_opt=zeros(n,1);
C_opt(select)=C(select);
y_opt=T*C_opt;
snr =20*log(norm( (y))/norm( (y)-y_opt)) 
figure;plot(y','g ','lineWidth',2);hold on; plot(real(y_opt),'r.-'); xlabel('Time second');ylabel('Rotor speed pu')
final_order=numel(select)
freq_spec=Freq(select)