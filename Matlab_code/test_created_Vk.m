clear  ; clc; close all;
% parameters
thres=0.99;gap=10; 

    load('line5-6.mat','bus_v','mac_ang','mac_spd');delay=1;t1=51+delay;t_step=0.01;T=100;  %line5-6.mat   line5-6.mat 3.11  line15-16.mat 4.41 line26-29.mat 5.74
  load('short56.mat','bus_v','mac_ang','mac_spd');delay=0;t1=71+delay;t_step=0.01;T=100; %short56.mat 26.84  short1516.mat 31.14 short26-29.mat 23.96
%   load('load56.mat ','bus_v','mac_ang','mac_spd');delay=1;t1=101+delay;t_step=0.01;T=100; % load56.mat 6.42 load1516.mat 6.07 load2629.mat 5.97

% normalize the data by subtracting the row mean
[s1,v]=Normalize_vol(bus_v,t1,T);
% choose rank
k1=choose_rank(s1,thres,gap);

%% construct Vk1 is the Vr in paper
Vk1=v(:,1:k1);
figure;plot((Vk1/norm(Vk1)),'o-');
    
%% Construct \Psi_r
%create eigenvalues of contitues a_mat
 load('a_line56.mat','a_mat','b_pm','c_v');
% load('a_line2629.mat','a_mat','b_pm','c_v');% simulate the line trip 5-6 event  a_line2629.mat a_line1516.mat a_line56.mat a_line2629.mat a_line1516.mat %load('a_load56.mat','a_mat','b_pm','c_v');%   simulate the load change evnet  
 a1=a_mat(1:32,1:32);b1=b_pm(1:32);  c1=c_v(:,1:32);               

 % choose the initial condition x0
 x0=chosex0(mac_ang,mac_spd,t1);      % dic is the number of dictionary of line trip, t2 is the time of fault

% calculate eigenvectors of a_mat
[Q eig_conti]=eig(a1);
% the discrete eigenvavlues and eigenvectors Q
eig_a1=expm(eig_conti*t_step);
eig_a1=diag(eig_a1); 
invQ=inv(Q);
sigma=eye(size(Q));
CQ=c1*Q;     
% construct \Psi= Vk_make   
jj=1;
n_steps=T ;
for i=1:numel(eig_a1)
   for j=1:n_steps
       Vk_make(i,jj)=eig_a1(i)^(j);
       jj=jj+1;
   end
   jj=1;
end
mainVk=[];
main_vcol=[];

% the ranks of \Psi_r
 k2=k1;  
    
% select the k1 dominant parts by the coefficient   
for i=1:size(Q,1)
    sigma(i,i)=abs(invQ(i,:)*x0) *norm(CQ(:,i))*norm(Vk_make(i,:));
end
sigma1=diag(sigma);
sigma2=abs(sigma1);  
max_sig=sort((  (sigma2 )),'descend');
main_vcol=max_sig(1:k2);
for i=1:k2
    index=find(abs (sigma2)==max_sig(i)); 
    indexVk(i) =index(1);
    sigma2(index(1))=0;
end

% choose the dominant rows of \Psi_r
for i=1:numel(indexVk)
     mainVk(i,:)=Vk_make(indexVk(i),:);
end
indexVk  % the chosen rows  
mainVk=mainVk';   
re_mainVk=real(mainVk);im_mainVk=imag(mainVk);
d=[re_mainVk im_mainVk];
d_sample=d(1:3:end,:); Vk1_sample=(Vk1(1:3:end,:));
figure;plot((d/norm(d) ),'*-');

%% compare the V_r and \Psi_r
subangle_d1=angle0( d_sample,Vk1_sample)  
 
