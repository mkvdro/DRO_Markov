function cost=FW_main(k,a,x_cur,epsilon,r,iter,q_all,alpha_all)
    %returns the cost vector k*1, each row represents 
    %the unnormalized cost value for the respective customer group

    %parameters
    d=length(q_all); % number of available brands

    % initialize cost
    cost=zeros(k,1); %negative profit
    % if r<=10^(-4)
    %     iter=iter*10;
    % end
    for i=1:k
        % fprintf('%dth chain',i);
        %specify the dynamics of ith markov chain
        q=q_all(:,:,i); %transition kernel
        alpha0=alpha_all(:,i); % stationary distribution
        nu=q;
        nu_best=zeros(d); % best worst-case transition kernel

        for t=1:iter % iterate within the specified iteration number
            % fprintf('%dth iter',t);
            c = grad_psi(a.*x_cur,nu); %define cost vector: notice that c is the gradient of the psi function  
            c_copy=c;
            c=reshape(c_copy,[d,d]);
            % disp(c-c');
            stt=linear_sub(alpha0,r,c,d,q); %subproblem, determining argmax for descent direction
            if isnan(stt)
                break
            end

            dir=stt-nu;
            g=trace(dir*c);
        
            if g<epsilon
                % fprintf('%d gap',g);
                nu_best=nu; %disp(nu_best);
                break
            end
            buf_lin_search=-10^9;
            gammat=-1;
            for gammax = 0:0.05:1
                nu_mat=nu+gammax*dir;
                nu_mat=nu_mat./sum(nu_mat,2);
                f_lin_search = Psi(a.*x_cur, nu_mat);
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
        
        if nu_best==zeros(d)
            nu_best=nu;
        end

        nu_best_mat=nu_best./sum(nu_best,2);

        approx_cost=Psi(a.*x_cur,nu_best);

        % disp(nu_best_mat);
        cost(i)=approx_cost;
        % disp(approx_cost);
    end
end