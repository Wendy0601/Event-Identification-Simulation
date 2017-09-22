function [bus_v,bus_freq ,name,ilf]=Current_readPMUdata(path,allnames,num_bus) 
    name=fullfile(path,allnames{1,num_bus});
    data=load(name);
    bus_v= data.bus_v;
    bus_freq=data.bus_freq;
    bus_theta=data.theta; 
    ilf=data.ilf;
end







% function [bus_v,bus_freq]=linedata(j)
%             
%     if j==1
%     load('/identification_data/outline/outline1.mat','theta','bus_v','bus_freq');
%    
%     end
%     if j==2
%     load('/identification_data/outline/outline2.mat','theta','bus_v','bus_freq');
%     
%     end
%     if j==3
%     load('/identification_data/outline/outline3.mat','theta','bus_v','bus_freq');
%     end
%     if j==4
%     load('/identification_data/outline/outline4.mat','theta','bus_v','bus_freq');
%     end    
%     if j==5
%     load('/identification_data/outline/outline5.mat','theta','bus_v','bus_freq');
%     end
%      if j==6
%     load('/identification_data/outline/outline6.mat','theta','bus_v','bus_freq');
%     end
%     if j==7
%     load('/identification_data/outline/outline7.mat','theta','bus_v','bus_freq');
%     end
%     if j==8
%     load('/identification_data/outline/outline8.mat','theta','bus_v','bus_freq');
%     end
%     if j==9
%     load('/identification_data/outline/outline9.mat','theta','bus_v','bus_freq');
%     end    
%     if j==10
%     load('/identification_data/outline/outline10.mat','theta','bus_v','bus_freq');
%     end
%      if j==11
%     load('/identification_data/outline/outline11.mat','theta','bus_v','bus_freq');
%     end
%     if j==12
%     load('/identification_data/outline/outline12.mat','theta','bus_v','bus_freq');
%     end
%     if j==13
%     load('/identification_data/outline/outline13.mat','theta','bus_v','bus_freq');
%     end
%     if j==14
%     load('/identification_data/outline/outline14.mat','theta','bus_v','bus_freq');
%     end    
%     if j==15
%     load('/identification_data/outline/outline15.mat','theta','bus_v','bus_freq');
%     end
%      if j==16
%     load('/identification_data/outline/outline16.mat','theta','bus_v','bus_freq');
%     end
%     if j==17
%     load('/identification_data/outline/outline17.mat','theta','bus_v','bus_freq');
%     end
%     if j==18
%     load('/identification_data/outline/outline18.mat','theta','bus_v','bus_freq');
%     end
%     if j==19
%     load('/identification_data/outline/outline19.mat','theta','bus_v','bus_freq');
%     end    
%     if j==20
%     load('/identification_data/outline/outline20.mat','theta','bus_v','bus_freq');
%     end
%      if j==21
%     load('/identification_data/outline/outline21.mat','theta','bus_v','bus_freq');
%     end
%     if j==22
%     load('/identification_data/outline/outline22.mat','theta','bus_v','bus_freq');
%     end
%     if j==23
%     load('/identification_data/outline/outline23.mat','theta','bus_v','bus_freq');
%     end
%     if j==24
%     load('/identification_data/outline/outline24.mat','theta','bus_v','bus_freq');
%     end    
%     if j==25
%     load('/identification_data/outline/outline25.mat','theta','bus_v','bus_freq');
%     end
%      if j==26
%     load('/identification_data/outline/outline26.mat','theta','bus_v','bus_freq');
%     end
%     if j==27
%     load('/identification_data/outline/outline27.mat','theta','bus_v','bus_freq');
%     end
%     if j==28
%     load('/identification_data/outline/outline28.mat','theta','bus_v','bus_freq');
%     end
%     if j==29
%     load('/identification_data/outline/outline29.mat','theta','bus_v','bus_freq');
%     end    
%     if j==30
%     load('/identification_data/outline/outline30.mat','theta','bus_v','bus_freq');
%     end
% % %     
% 
% 
%     if j==31
%     load('/identification_data/outline/outline31.mat','theta','bus_v','bus_freq');
%     end
%     if j==32
%     load('/identification_data/outline/outline32.mat','theta','bus_v','bus_freq');
%     end
%     if j==33
%     load('/identification_data/outline/outline33.mat','theta','bus_v','bus_freq');
%     end
%     if j==34
%     load('/identification_data/outline/outline34.mat','theta','bus_v','bus_freq');
%     end    
%     if j==35
%     load('/identification_data/outline/outline35.mat','theta','bus_v','bus_freq');
%     end
%      if j==36
%     load('/identification_data/outline/outline36.mat','theta','bus_v','bus_freq');
%     end
%     if j==37
%     load('/identification_data/outline/outline37.mat','theta','bus_v','bus_freq');
%     end
%     if j==38
%     load('/identification_data/outline/outline38.mat','theta','bus_v','bus_freq');
%     end
%     if j==39
%     load('/identification_data/outline/outline39.mat','theta','bus_v','bus_freq');
%     end    
%     if j==40
%     load('/identification_data/outline/outline40.mat','theta','bus_v','bus_freq');
%     end
%      if j==41
%     load('/identification_data/outline/outline41.mat','theta','bus_v','bus_freq');
%     end
%     if j==42
%     load('/identification_data/outline/outline42.mat','theta','bus_v','bus_freq');
%     end
%     if j==43
%     load('/identification_data/outline/outline43.mat','theta','bus_v','bus_freq');
%     end
%     if j==44
%     load('/identification_data/outline/outline44.mat','theta','bus_v','bus_freq');
%     end    
%     if j==45
%     load('/identification_data/outline/outline45.mat','theta','bus_v','bus_freq');
%     end
%      if j==46
%     load('/identification_data/outline/outline46.mat','theta','bus_v','bus_freq');
%     end
%     if j==47
%     load('/identification_data/outline/outline47.mat','theta','bus_v','bus_freq');
%     end
%     if j==48
%     load('/identification_data/outline/outline48.mat','theta','bus_v','bus_freq');
%     end
%     if j==49
%     load('/identification_data/outline/outline49.mat','theta','bus_v','bus_freq');
%     end    
%    if j==50
%     load('/identification_data/outline/outline50.mat','theta','bus_v','bus_freq');
%     end
%     if j==51
%     load('/identification_data/outline/outline51.mat','theta','bus_v','bus_freq');
%     end
%     if j==52
%     load('/identification_data/outline/outline52.mat','theta','bus_v','bus_freq');
%     end
%     if j==53
%     load('/identification_data/outline/outline53.mat','theta','bus_v','bus_freq');
%     end
%     if j==54
%     load('/identification_data/outline/outline54.mat','theta','bus_v','bus_freq');
%     end    
%     if j==55
%     load('/identification_data/outline/outline55.mat','theta','bus_v','bus_freq');
%     end
%      if j==56
%     load('/identification_data/outline/outline56.mat','theta','bus_v','bus_freq');
%     end
%     if j==57
%     load('/identification_data/outline/outline57.mat','theta','bus_v','bus_freq');
%     end
%     if j==58
%     load('/identification_data/outline/outline58.mat','theta','bus_v','bus_freq');
%     end
%     if j==59
%     load('/identification_data/outline/outline59.mat','theta','bus_v','bus_freq');
%     end  
%     if j==60
%     load('/identification_data/outline/outline60.mat','theta','bus_v','bus_freq');
%     end  
%      if j==61
%     load('/identification_data/outline/outline61.mat','theta','bus_v','bus_freq');
%     end
%     if j==62
%     load('/identification_data/outline/outline62.mat','theta','bus_v','bus_freq');
%     end
%     if j==63
%     load('/identification_data/outline/outline63.mat','theta','bus_v','bus_freq');
%     end
%     if j==64
%     load('/identification_data/outline/outline64.mat','theta','bus_v','bus_freq');
%     end    
%     if j==65
%     load('/identification_data/outline/outline65.mat','theta','bus_v','bus_freq');
%     end
%      if j==66
%     load('/identification_data/outline/outline66.mat','theta','bus_v','bus_freq');
%     end
%     if j==67
%     load('/identification_data/outline/outline67.mat','theta','bus_v','bus_freq');
%     end
%     if j==68
%     load('/identification_data/outline/outline68.mat','theta','bus_v','bus_freq');
%     end
%     if j==69
%     load('/identification_data/outline/outline69.mat','theta','bus_v','bus_freq');
%     end    
%     if j==70
%     load('/identification_data/outline/outline70.mat','theta','bus_v','bus_freq');
%     end
%      if j==71
%     load('/identification_data/outline/outline71.mat','theta','bus_v','bus_freq');
%     end
%     if j==72
%     load('/identification_data/outline/outline72.mat','theta','bus_v','bus_freq');
%     end
%     if j==73
%     load('/identification_data/outline/outline73.mat','theta','bus_v','bus_freq');
%     end
%     if j==74
%     load('/identification_data/outline/outline74.mat','theta','bus_v','bus_freq');
%     end    
%     if j==75
%     load('/identification_data/outline/outline75.mat','theta','bus_v','bus_freq');
%     end
%      if j==76
%     load('/identification_data/outline/outline76.mat','theta','bus_v','bus_freq');
%     end
%     if j==77
%     load('/identification_data/outline/outline77.mat','theta','bus_v','bus_freq');
%     end
%     if j==78
%     load('/identification_data/outline/outline78.mat','theta','bus_v','bus_freq');
%     end
%     if j==79
%     load('/identification_data/outline/outline79.mat','theta','bus_v','bus_freq');
%     end    
%     if j==80
%     load('/identification_data/outline/outline80.mat','theta','bus_v','bus_freq');
%     end
%     if j==81
%     load('/identification_data/outline/outline81.mat','theta','bus_v','bus_freq');
%     end
%     if j==82
%     load('/identification_data/outline/outline82.mat','theta','bus_v','bus_freq');
%     end
%     if j==83
%     load('/identification_data/outline/outline83.mat','theta','bus_v','bus_freq');
%     end
% 
% end
%             
            
            
