clear ;close all;clc;
% parameters
gap=10; thres=0.99;  
event=7;
%choose datasets for each type of events
% event =1: linetrip; 2:generator trip; 3: short circuit; 4: line to line faults 5:load change;
% 6: capacitor bank switch on; 7: motor off;  
if event ==1 
    path='H:\identification_data\Line_Trip_classic4_86';%'H:\identification_data\Line_trip_86';
elseif event==2
    path='H:\identification_data\Generator_trip_classic_86';
elseif event==3
    path='H:\identification_data\Three_Phase_86';
elseif event==4
    path='H:\identification_data\Line_to_Line_classic_86';
elseif event==5
    path='H:\identification_data\Load_change_classic_86';
elseif event==6
    path='H:\identification_data\Capacitor_Bank_classic_86';%Capacitor_Bank_classic_86';
elseif event==7
    path='H:\identification_data\Motor_Off_classic_86';   
end

% compute subspace angle
[ allnames, len ] = read_all_file( path,pwd );
for i=1:len
    for j=1:len
        %dataset 1
        [bus_v1,bus_freq1 ,name1]=readPMUdata(path,allnames,i); 
        if event ==2
           t01=51;t02=151;
           X1=[abs(bus_v1(1:51+i,:)) ; abs(bus_v1(min(53+i,68):68,:))];
        elseif (event==3 || event==4)
           t01=72;t02=172; 
           X1=abs(bus_v1(1:68,:));
        else
           t01=51;t02=151;  
           X1=abs(bus_v1(1:68,:));
        end
        %dataset 2
        [bus_v2,bus_freq2 ,name2]=readPMUdata(path,allnames,j); 
        if event ==2
           t01=51;t02=151;
           X2=[abs(bus_v2(1:51+j,:)) ; abs(bus_v2(min(53+j,68):68,:))];
        elseif (event==3 || event==4)
           t01=71;t02=171; 
           X2=abs(bus_v2(1:68,:));
        else
           t01=51;t02=151;  
           X2=abs(bus_v2(1:68,:));
        end
        %compute subspace angle
        X1=sub_rowmean(X1,t01); 
        [U1,S1,V1] = svd(X1(:,t01:3:t02) );%%%%%%%%%%%%%%%%%%%%% here the sampling time is 33/seond 
        s1=diag(S1);
        k1(i)=choose_rank(s1,thres,gap);   
        E1(i)=sum(s1(1:k1(i)));

         X2=sub_rowmean(X2,t01); 
        [U2,S2,V2] = svd(X2(:,t01:3:t02) );%%%%%%%%%%%%%%%%%%%%% here the sampling time is 33/seond 
        s2=diag(S2);
        k2(j)=choose_rank(s2,thres,gap);   
        E2(j)=sum(s2(1:k2(j)));
 
        k12=max(k1(i),k2(j));
        SubspaceAngles(i,j)=angle0(V1(:,1:k1(i)),V2(:,1:k2(j)));  
    end
end


path=pwd;
cd('H:\identification_data\Subspace_Angle');
save  angle_Capacitor_Bank_5pu_classic_86.mat;
cd(path)