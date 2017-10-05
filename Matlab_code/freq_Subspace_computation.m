clear ;close all;clc;
% parameters
gap=10; thres=0.99;   
T=100;
%choose datasets for each type of events
% event =1: linetrip; 2:generator trip; 3: short circuit; 4: line to line faults 5:load change;
% 6: capacitor bank switch on; 7: motor off;  
for event=3
    SubspaceAngles=[];
    if event ==1 
        path='H:\identification_data\Line_Trip_classic4_86';%'H:\identification_data\Line_trip_86';
    elseif event==2
        path='H:\identification_data\Generator_trip_classic_86';
    elseif event==3
        path='H:\identification_data\Three_Phase_classic_86';
    elseif event==4
        path='H:\identification_data\Load_change_classic_86';
    elseif event==5
        path='H:\identification_data\Line_to_ground_classic_86';%Capacitor_Bank_classic_86';%Capacitor_Bank_classic_86';
    end

    % compute subspace angle
    [ allnames, len ] = read_all_file( path,pwd );     
for i=1:len
    for j=1:len
        %dataset 1
        [bus_v1,bus_freq1 ,name1]=readPMUdata(path,allnames,i); 
        if event ==2
           t01=53;t02=153;
           if i<16
               X1=[abs(bus_freq1(1:51+i,:)) ; abs(bus_freq1(min(53+i,68):68,:))]; 
           else
               X1=[abs(bus_freq1(1:67,:))];
           end
        elseif (event==3 | event==5 )
           t01=73;t02=173; 
           X1=abs(bus_freq1(1:68,:));
        elseif (event==4  )
           t01=103;t02=203; 
           X1=abs(bus_freq1(1:68,:));
        else
           t01=53;t02=153; 
           X1=abs(bus_freq1(1:68,:));
        end
        %dataset 2
        [bus_v2,bus_freq2 ,name2]=readPMUdata(path,allnames,j); 
        if event ==2
           t01=53;t02=153;
           if j<16 
               X2=[abs(bus_freq2(1:51+j,:)) ; abs(bus_freq2(min(53+j,68):68,:))]; 
           else
               X2=[abs(bus_freq2(1:67,:)) ];
           end
        elseif (event==3 ||event==5 )
           t01=73;t02=173; 
           X2=abs(bus_freq2(1:68,:));
        elseif (event==4  )
           t01=103;t02=203; 
           X2=abs(bus_freq2(1:68,:));
        else
           t01=53;t02=153;
           X2=abs(bus_freq2(1:68,:));
        end
        %compute subspace angle
        X1=sub_rowmean(X1,t01); 
        [U1,S1,V1] = svd(X1(:,t01:3:t02) );%%%%%%%%%%%%%%%%%%%%% here the sampling time is 33/seond 
        s1=diag(S1);
        k1(i)=choose_rank(s1,thres,gap);   
        E1(i)=sum(s1(1:k1(i)));
%         figure;plot(V1(:,1:k1(i)));
%         figure; plot(s1,'o'); k1(i)

         X2=sub_rowmean(X2,t01); 
        [U2,S2,V2] = svd(X2(:,t01:3:t02) );%%%%%%%%%%%%%%%%%%%%% here the sampling time is 33/seond 
        s2=diag(S2);
        k2(j)=choose_rank(s2,thres,gap);   
        E2(j)=sum(s2(1:k2(j)));
 
        k12=max(k1(i),k2(j));

        SubspaceAngles(i,j)=angle0(V1(:,1:k12),V2(:,1:k12));  
        end
    end

% % 
%     if event ==1 
%         save H:\identification_data\Subspace_Angle\freq_Line_Trip_classic4_86.mat
%     elseif event==2
%         save H:\identification_data\Subspace_Angle\freq_Generator_trip_classic_86.mat
%     elseif event==3
%         save H:\identification_data\Subspace_Angle\freq_Three_Phase_classic_86.mat
%     elseif event==4
%         save H:\identification_data\Subspace_Angle\freq_Load_change_classic_86.mat
%      elseif event==5
%         save H:\identification_data\Subspace_Angle\freq_Line_to_ground.mat    
%     end 

end


%plot(V1(:,1:k1));figure; plot(X1');figure ; plot(X2'); figure; plot(abs(ilf1'))