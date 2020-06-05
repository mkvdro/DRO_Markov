function s = linear_sub(alpha0,r,c,m,q)

cbar=zeros(1,m);
beta=zeros(1,m);
% disp('alpha0');
% disp(alpha0);
for i=1:m
    cbar(i)=max(c(i,:)); 
end

% disp('cbar');disp(cbar);
sum2=0;
for i=1:m
    sum2=sum2+cbar(i)-exp(-r)*sum(q(i,:).*c(i,:));
end
for i=1:m
    beta(i)=(sum2)/(1-exp(-r))-sum(cbar)+cbar(i); %upper bound for eta_dual(i)
    if isnan(beta(i))
        beta(i)=10^7;
%         error('beta value is nan');
    end
end
eta_0=zeros(1,m);

for i=1:m
        eta_0(i)=cbar(i)+10^(-19); %sample within the bound; (beta(i)-cbar(i))*rand +
end

su = @(eta_dual) 0;
    for j=1:m
        for k=1:m
%             disp('difference betw
            if alpha0(j)~=0
            su =@(eta_dual) su(eta_dual) + alpha0(j) * log((eta_dual(j)-c(j,k))/alpha0(j)) * q(j,k);
%             disp('on the fly sum value inside lambda');
%             disp(su(eta_0));
            end
        end
    end
    lambda =@(eta_dual) exp(su(eta_dual)-r); %compute lambda
    
    f = @(eta_dual) sum(eta_dual) -lambda(eta_dual);
y=f(eta_0);
if isnan(y)
    disp('nan f');
    disp(cbar);disp(beta);
    s=nan;
    return
end
% disp('initial value');
% disp(y);
% disp('initial value for lambda');    
% disp(lambda(eta_0));
% disp('initial value for sum inside lambda');    
% disp(su(eta_0));
% disp('initial value for f');    
% disp(f(eta_0));
opts = optimoptions('fmincon','Display','off');
eta_star = fmincon(f,eta_0,[],[],[],[],cbar+10^(-20),beta,[],opts);

sum1=0;
for i=1:m
    for j=1:m
        if (alpha0(j)~=0) && ((eta_star(j)-c(j,k))~=0)
        sum1 = sum1 + alpha0(j) * log((eta_star(j)-c(j,k))/alpha0(j)) * q(j,k);            
        end
    end
end
lambda_star=exp(sum1-r); %recover lambda

s=zeros(1,m^2);

for i=1:m
%     disp('pi for q');
%     disp(alpha0(i));
    for j=1:m
        if alpha0(i)==0
            s(m*(i-1)+j)=0;
        else
        s(m*(i-1)+j)=lambda_star*alpha0(i)*q(i,j)/(eta_star(i)-c(i,j));
        if isnan(s(m*(i-1)+j))
         error('optimizer nan');
        end
        end
    end
end

end

