%clc; clear all; close all;
%a=[4 8 6 -1 -2 -3 -1 3 4 5];
%K=3;
function b=MedianFilter(a,K1,K2)
[row, col]=size(a);
for i=1:row
    for j=1:col
        if j<=K1
            b(i,j)= mean(a(i,1:K1) );
        elseif (j+K2)<col 
            b(i,j)=mean(a(i,j:j+K2));
        else
%             b(i,j)=mean(a(i,j:j+K));
            b(i,j)=mean(a(i,j-K2+1:j));
        end
    end
end
return 