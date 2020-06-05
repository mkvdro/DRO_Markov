function cost = cost_noM(alpha0,r,m)

cbar=m-1; % max(L(x,xi)) = max(xi)=0
beta = (cbar-exp(-r)*dot(0:m-1,alpha0))/(1-exp(-r));

su = @(eta_dual) 0;
for j=1:(m-1)
    su =@(eta_dual) su(eta_dual) + alpha0(j) * log(eta_dual-(j-1));
%             disp('on the fly sum value inside lambda');
%             disp(su(eta_0));
end
lambda =@(eta_dual) exp(su(eta_dual)-r); %compute lambda
    
f = @(eta_dual) eta_dual -lambda(eta_dual);
% disp('initial value for lambda');    
% disp(lambda(eta_0));
% disp('initial value for sum inside lambda');    
% disp(su(eta_0));
% disp('initial value for f');    
% disp(f(eta_0));
% opts = optimset('fminbnd','Display','off');
[eta_star,cost] = fminbnd(f,cbar,beta);

end

