function [val0,q_all] = est_alpha_from_xi(d,m,T,xi)

%estimator theta_T
theta_T=zeros(m,m,d);
for i=1:d 
theta_q=zeros(m);
for k=1:m
    for j=1:m
        sum1=0;
        for t=1:(T-1)
            sum1=sum1+(xi(i,t)==k-1)*(xi(i,t+1)==j-1);
        end
        
        if sum1==0
            theta_q(k,j)=10^(-19);
        else
        theta_q(k,j)=sum1/(T-1);
        end        
    end
end
theta_T(:,:,i)=theta_q;
end

%estimate transition matrix q_ij
q_all=zeros(m,m,d);
for i=1:d 
q=zeros(m);
for k=1:m
    for j=1:m
%         if (sum(theta_T(k,1:m,i))==0)
%             error('division by zero computing q'); 
%             disp(theta_T(k,1:m,i)); disp(theta_T);return
%             q(i,j)=theta_q(i,j)/0.00001;
%         else
            q(k,j)=theta_T(k,j,i)/sum(theta_T(k,1:m,i));
%         end
    end
end
q_all(:,:,i)=q;
end

val0=zeros(m,1,d);
for i=1:d
    val=zeros(m,1);
for k=1:m
    val(k)=sum(theta_T(k,1:m,i));
end
val0(:,:,i)=val;
end
end