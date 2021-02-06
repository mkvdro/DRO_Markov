function s = linear_sub(alpha0,r,c,m,q)
    %outputs arg max <S, c> where c is the gradient
    % solved using duality

    cbar=max(c,[],2); 
    eta_0=cbar+10;

    su = @(eta_dual) alpha0'*(sum(q.*(log(repmat( eta_dual, [1,m] )-c)-log(repmat( alpha0, [1,m] ))),2));

    lambda =@(eta_dual) exp(su(eta_dual)-r); %compute lambda
        
    f = @(eta_dual) sum(eta_dual) -lambda(eta_dual);
    y=f(eta_0);
    if isnan(y)
        disp('nan f');
        disp(cbar);
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
    options = optimoptions(opts,'MaxFunctionEvaluations', max(1000,round(10^3*r)));
    options = optimoptions(options,'MaxIterations',  max(1000,round(10^3*r)));
    options = optimoptions(options,'OptimalityTolerance', 10^(-3));
    options = optimoptions(options,'FunctionTolerance', 10^(-3));
    options = optimoptions(options,'StepTolerance', 10^(-3));
    eta_star = fmincon(f,eta_0,[],[],[],[],cbar,[],[],options);
    
    % disp(eta_star-cbar);
    if (eta_star-cbar>=ones(size(eta_0)))
        s=nan;
        return 
    end
    sum2 = alpha0'*(sum(q.*(log(repmat( eta_star, [1,m] )-c)-log(repmat( alpha0, [1,m] ))),2));
    lambda_star = exp(sum2-r);
    % disp(lambda_star);disp(eta_star);
    
    s=zeros(m);
    
    
    for i=1:m
    %     disp('pi for q');
    %     disp(alpha0(i));
        for j=1:m
            if alpha0(i)~=0
            s(i,j)=lambda_star*alpha0(i)*q(i,j)/(eta_star(i)-c(i,j));
            if isnan(s(i,j))
                s(i,j)=10^(-2);
            % error('optimizer nan');
            end
            end
        end
    end
    
    end
    
    
    