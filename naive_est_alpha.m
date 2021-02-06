function q_T = naive_est_alpha(k,d,T,xi)

%estimator theta_T

q_T=zeros(k,d);
for i=1:k 
    % idx=[]; %missing observation
    for j=1:d
    % count=zeros(d,1);
    sum1=0;
    
    for t=1:(T)
        sum1=sum1+(xi(i,t)==j);
    end
    if sum1==0
        q_T(i,j)=1/10^5;
    else
    q_T(i,j)=sum1/(T);
    end        
end
end

