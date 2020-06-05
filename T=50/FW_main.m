function cost= FW_main(epsilon,r,iter,q,m,alpha0)

%initialize nu

initial_cost=5;approx_cost=4;
udx=1;
while approx_cost<initial_cost
    g=zeros(1,iter);
nu=reshape((sample(alpha0,r,q,m)).',1,[]);
c=zeros(m,m);
nu_best=zeros(1,m^2);

initial_cost=dot(0:(m-1), alpha0);
disp(udx);
for t=1:iter
    fprintf("iter %d",t);
    cv=Grad_Psi(m, nu); %define cost vector     
    for i=1:m
        for j=1:m
            c(i,j)=cv(m*(i-1)+j); %cost vector in matrix form
        end
    end
    stt=linear_sub_SGD(alpha0,r,c,m,q); %subproblem
    if isnan(stt)
        break
    end
    dir=stt-nu;
    g(t)=dot(dir,(cv));
    
    if g(t)<epsilon
        nu_best=nu; %disp(nu_best);
        break
    end
    
    buf_lin_search=-10^9;
    gammat=-1;
    for gammax = 0:0.0001:1
        nu_buf=nu+gammax*dir;
        nu_mat=zeros(m);
        for i=1:m
            for j=1:m
                nu_mat(i,j)=nu_buf(m*(i-1)+j);

            end
        end
        if sum(isinf(nu_mat(:)))>=1 || sum(isnan(nu_mat(:)))>=1
           udx=udx+1; break
        end
        nu_mat=nu_mat./sum(nu_mat,2);
        nu_buf=reshape(nu_mat.',1,[]);
        f_lin_search = Psi(m, nu_buf);
        if f_lin_search > buf_lin_search
            gammat = gammax;
            buf_lin_search = f_lin_search;
        end
    end
    if gammat==-2
        break
    end
    nu=nu+gammat*dir;
    if sum(isinf(nu(:)))>=1 || sum(isnan(nu(:)))>=1
            break
    end
end
if nu_best==zeros(1,m^2)
   nu_best=nu;
end

nu_best_mat=zeros(m);
for i=1:m
    for j=1:m
        nu_best_mat(i,j)=nu_best(m*(i-1)+j);
    end
end
nu_best_mat=nu_best_mat./sum(nu_best_mat,2);


I = eye(size(nu_best_mat));
if sum(isinf(nu_best_mat(:)))>=1 || sum(isnan(nu_best_mat(:)))>=1
    initial_cost=5;approx_cost=4;udx=udx+1;
else
Y = null(nu_best_mat'-I);
pi_approx = Y./(sum(Y));

approx_cost=dot(0:(m-1), pi_approx);
cost=approx_cost;
end
end
end