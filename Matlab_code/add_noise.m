function [ X_noise ] = add_noise( X,noise,noise_level )
 if noise==1
        [row,col]=size(X);
        for s=1:row
            X_noise(s,:)=awgn(X(s,:),noise_level,'measured');
        end
 else 
    X_noise=X;
 end
end

