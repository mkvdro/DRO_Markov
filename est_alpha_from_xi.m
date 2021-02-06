function [alpha0,q_all] = est_alpha_from_xi(k,d,T,xi) % k segments, d brands, T time span, xi sample data (size: d* T * k)
    %estimator theta_T
    theta_T=zeros(d,d,k);
    for i=1:k
    theta_q=zeros(d);
    for l=1:d
        for j=1:d
            sum1=0;
            for t=1:(T-1)
                sum1=sum1+(xi(i,t)==l)*(xi(i,t+1)==j); % sequence of observation
            end
            
            if sum1==0
                theta_q(l,j)=1/10^5; % puts equal mass on unobserved sequence
            else
            theta_q(l,j)=sum1/(T-1);
            end        
        end
    end
    theta_T(:,:,i)=theta_q;
    end

    %estimate transition matrix q hat
    q_all=zeros(d,d,k);
    for i=1:k
    q=zeros(d);
    q(:,:)=theta_T(:,:,i)./sum(theta_T(:,1:d,i)); % marginalizing out the 
    q_all(:,:,i)=q;
    end

    alpha0=zeros(d,k); % d number of brands (state), k number of customer segment (different chains)
    % stationary distribution is a column vector

    for i=1:k
        val=zeros(d,1);
        for j=1:d
            val(j)=sum(theta_T(j,1:d,i));
        end
        % mc=dtmc(q_all(:,:,i));
        % buf=asymptotics(mc); % return stationary distribution (row vector) of the corresponding transition matrix
        alpha0(:,i)=val;
    end
end