function cost_val = cost_noM(a,k,x,alpha_all,r,d)
cost_val=zeros(k,1);
cost=-a.*x;
for i=1:k
alpha0=alpha_all(i,:);
% disp(alpha0);
cbar=max(cost); % 
beta0 = 10^7;%(cbar-exp(-r)*dot(cost,alpha0))/(1-exp(-r));
% disp(cbar);
% disp(beta0);
su = @(eta_dual) 0;
for j=1:d
    su =@(eta_dual) su(eta_dual) + alpha0(j) * log(eta_dual-cost(j));
end
lambda =@(eta_dual) exp(su(eta_dual)-r); %compute lambda
    
f = @(eta_dual) eta_dual -lambda(eta_dual);
% opts = optimset('fminbnd','Display','off');
[eta_star,cost_val(i)] = fminbnd(f,cbar,beta0);
end
end

