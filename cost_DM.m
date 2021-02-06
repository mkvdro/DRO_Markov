function cost = cost_DM(a,x_cur,r,q_all,k,d)
cost=zeros(k,1);
for l=1:k
    q=q_all(:,:,l);
    h=10^(-3); %grid size
    p=q;
    N=10^2;
    
    N_mon=10^2;
    cost_lower=-10^6;
    % Discretize the space of all probability simplex
    for sample=1:N_mon
 %       disp('r');disp(r);
        for n=1:N % N: number of all possible realizations of hat P := Q
            % generate 
            ind=zeros(d,1);
            % fprintf('sample order %d',n);
            
            if r>=1
            %     p=rand(d);
                for j=1:d
                    for kk=1:d
                        p(j,kk) = max(q(j,kk)+0.05*rand,0.01); 
                    end
                end
            elseif r>=0.1
                for j=1:d
                    for kk=1:d
                        p(j,kk) = max(q(j,kk)+0.0001*rand,0.001);
                    end
                end
        
            elseif r>=10^(-2)
                for kk=1:d
                    row=randi([1,d]);
                    p(row,kk) = max(q(row,kk)+0.0001*rand,0.001);
                end

            else
%                 for j=1:d
                    % for kk=1:d
                        row=randi(d);
                        col=randi(d);
                        diff=rand;
                        p(row,col) = q(row,col)+0.0001*r*diff;
                        p(col,row)=q(col,row)-0.0001*r*diff;
                    % end
%                 end
            end
            p = p./sum(p,2);
%            disp(p-q);
            for i=1:d %compute wasserstein distance for each row of Q


            %     disp(Aeq);disp(A);
            %     error('rest');
            %     opts = optimoptions('linprog','Display','off'); 
            %opts=optimoptions('fmincon','Display','off'); 
            dw=ws_distance(p,q,2);
            %   [opt_pm,dw]=linprog(f,A,ones(1,2*d),Aeq,beq,zeros(1,d^2),ones(1,d^2));%,opts);
            %   [opt_pm,dw,exitflag] = fmincon(@(x) dot(f',x),pm,A,ones(1,2*d),Aeq,beq,zeros(1,d^2),ones(1,d^2),[],opts);

            %     disp(opt_p);
                
            %    fprintf('dw %d',dw);
                
                if dw<=r
                    ind(i)=1;
                end
            end
    
            if ind==1
                break
            end
    
        end
    
        if min(abs(ind))==0
            error('no correct sample');
        end
    
        % Optimize Psi(p)
        buf=Psi(a.*x_cur,p);
        
        if buf>cost_lower
            cost(l)=buf;
            p_t=p;
        end
       % disp(cost(l));
    end
        
%    disp('fina');
%    disp(Psi(a.*x_cur,p_t));
%    disp(p_t);
end
end    
    