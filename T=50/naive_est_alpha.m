function q_T = naive_est_alpha(d,m,T,xi)

%estimator theta_T
q_T=zeros(d,m);
for i=1:d 
for k=1:m
        sum1=0;
        for t=1:(T)
            sum1=sum1+(xi(i,t)==k-1);
        end
        if sum1==0
            q_T(i,k)=10^(-9);
        else
        q_T(i,k)=sum1/(T);
        end        
end
end
end