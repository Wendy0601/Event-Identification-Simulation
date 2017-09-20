% %%%%%%%%%%% the change of initial condition
% 
% %%%%%% power flow difference ratio and its plot
load('C:\07_research\condi_b_power.mat','P');
P1=P(:,10);
load('C:\07_research\condi_d_power.mat','P'); Bus_loss=12; %eta=3.3
% load('C:\07_research\condi_e_power.mat','P'); Bus_loss=19;%eta=28.9
% load('C:\07_research\condi_f_power.mat','P'); Bus_loss=58;%eta=13.4
% load('C:\07_research\condi_g_power.mat','P'); Bus_loss=10;%eta=20.7
% load('C:\07_research\condi_h_power.mat','P'); Bus_loss1=6;Bus_loss2=19;
% P2=P(:,10);
% for i=1:82
%     if i<Bus_loss
%          meanp(i)= abs((P1(i,1)-P2(i,1))/abs(P1(i,1)));
%     elseif i>Bus_loss
%         meanp(i)= abs((P1(i+1,1)-P2(i,1))/abs(P1(i+1,1)));
%     else
%        meanp(i)=1;% continue%
%     end
%         
% end
% meanp1=mean(meanp(find(meanp<2)))
% plot(P1(1:4,10:20)','-.*','Linewidth',2)
% figure; plot(P2(5:8,10:20)','--o','Linewidth',2)

P2=P(:,10);
for i=1:81
    if i<Bus_loss1
         meanp(i)= abs((P1(i,1)-P2(i,1))/abs(P1(i,1)));
    elseif i>Bus_loss1 && i < Bus_loss2
        meanp(i)= abs((P1(i+1,1)-P2(i,1))/abs(P1(i+1,1)));
    elseif i > Bus_loss2
       meanp(i)=abs((P1(i+2,1)-P2(i,1))/abs(P1(i+2,1)));
    else
        meanp(i)=1;
    end
     
        
end
meanp1=mean(meanp(find(meanp<2)))
