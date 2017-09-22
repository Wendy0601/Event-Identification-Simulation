clear ;close all;clc;
% parameters
gap=0; thres=0.99;   
%choose datasets for each type of events
% event =1: linetrip; 2:generator trip; 3: short circuit; 4: line to line faults 5:load change;
% 6: capacitor bank switch on; 7: motor off;  
for event=4
    SubspaceAngles=[];
    if event ==1 
        path='H:\identification_data\Line_Trip_classic4_86';%'H:\identification_data\Line_trip_86';
    elseif event==2
        path='H:\identification_data\Generator_trip_classic_86';
    elseif event==3
        path='H:\identification_data\Three_Phase_classic_86';
    elseif event==4
        path='H:\identification_data\condition_c\Load_change_classic_86';
    elseif event==5
        path='H:\identification_data\Line_to_ground_classic_86';%Capacitor_Bank_classic_86';%Capacitor_Bank_classic_86';
    end

    % compute subspace angle
    [ allnames, len ] = read_all_file( path,pwd );
    load('H:\identification_data\Generator_trip_classic_86\linetrip_2-53.mat','line');
    bus_line=line;       
    for i=1:len
        for j=1:len
            %dataset 1
            [bus_v1,bus_freq1 ,name,ilf1]=Current_readPMUdata(path,allnames,i);
            if event ==2
               t01=51;t02=151;
               lossbus=find(bus_line(:,2)==52+i);
               if i>16
                   i=i-16;
                   X1=[abs(ilf1(1:lossbus-1,:)) ; abs(ilf1(min(lossbus+1,83):83,:))];
                   i=i+16;
               else
                   X1=[abs(ilf1(1:lossbus-1,:)) ; abs(ilf1(min(lossbus+1,83):83,:))];
               end
            elseif (event==3 | event==5 )
               t01=71;t02=171; 
               X1=abs(ilf1(1:83,:));
            elseif (event==4  )
               t01=101;t02=201; 
               X1=abs(ilf1(1:83,:));
            else
               t01=51;t02=151; 
               X1=abs(ilf1(1:83,:));
            end
            %dataset 2
            [bus_v2,bus_freq2 ,name2,ilf2]=Current_readPMUdata(path,allnames,j);
            if event ==2
               t01=51;t02=151;
               lossbus=find(bus_line(:,2)==52+j);
               if j>16
                   j=j-16;
                   X2=[abs(ilf2(1:lossbus-1,:)) ; abs(ilf2(min(lossbus+1,83):83,:))];
                   j=j+16;
               else
                   X2=[abs(ilf2(1:lossbus-1,:)) ; abs(ilf2(min(lossbus+1,83):83,:))];
               end
            elseif (event==3 ||event==5 )
               t01=71;t02=171; 
               X2=abs(ilf2(1:83,:));
            elseif (event==4  )
               t01=101;t02=201; 
               X2=abs(ilf2(1:83,:));
            else
               t01=51;t02=151;
               X2=abs(ilf2(1:83,:));
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


    if event ==1 
        save H:\identification_data\Subspace_Angle\CurrentMagnitudes_Line_Trip_classic4_86.mat
    elseif event==2
        save H:\identification_data\Subspace_Angle\CurrentMagnitudes_Generator_trip_classic_86.mat
    elseif event==3
        save H:\identification_data\Subspace_Angle\CurrentMagnitudes_Three_Phase_classic_86.mat
    elseif event==4
        save H:\identification_data\Subspace_Angle\condition_c_CurrentMagnitudes_Load_change_classic_86.mat
     elseif event==5
        save H:\identification_data\Subspace_Angle\CurrentMagnitudes_Line_to_ground.mat    
    end 

end


%plot(V1(:,1:k1));figure; plot(X1');figure ; plot(X2'); figure; plot(abs(ilf1'))