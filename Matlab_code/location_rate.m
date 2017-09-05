function [location]=location_rate(bus_v,t1,kmax,threshold )
X1= [abs(bus_v(1:68,:))];
[row,col]=size(X1);
for i=1:row       
    X1(i,:)=X1(i,:)-mean(X1(i,1:t1-2))*ones(1,col);            
end 
[U,S,V] = svd(X1(:,t1:3:t1+10) );
s = diag(S);
sum_s=0;
s_k=0;
while sum_s<threshold*sum(s)
    sum_s=sum_s+s(s_k+1);
    s_k=s_k+1;
end

ku=s_k; 
Ukk=U(:,1:ku)*S(1:ku,1:ku);%*U(:,1:ku)';
wei=zeros(row,1);
for i=1:row
    wei(i)=norm((Ukk(i,:)));
end
wei=wei/sum(wei);
location=find_k_max(wei,kmax);
% location(1,1)=find(wei==max(wei))  ;%%%%%%% pick 
% wei(location(1,1))=0;
% maxw=find(wei==max(wei)) ;
% location(1,2:1+numel(maxw))=maxw;
% wei(location(1,2:1+numel(maxw)))=0;
% location(1,3)=find(wei==max(wei)) ;
% wei(location(1,3))=0;
% location(1,4)=find(wei==max(wei)) ;
% wei(location(1,4))=0;
% location(1,5)=find(wei==max(wei)) ;
% figure;bar(wei)
end




