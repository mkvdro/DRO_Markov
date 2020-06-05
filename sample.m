function val = sample(alpha0,r,q,d)
% alpha0=alpha0(:,:,1);
% d=4;
% q=q(:,:,1);
p=q;
N_samples = 10^7;
for n=1:N_samples
    fprintf('sample order %d',n);
    % randomly sample p
%     if r>=0.5
%         for j=1:d
%             for k=1:d
%             p(j,k) = max(q(j,k)+1.5*r*rand-1.5*r,10^(-11)); 
%             end
%         end
%     else
    if r>=0.1
        for j=1:d
        for k=1:d
            p(j,k) = max(q(j,k)+0.1*r*rand-r,10^(-15)); 
        end
        end
    elseif r>=0.07
    for j=1:d
    for k=1:d
        p(j,k) = max(q(j,k)+0.01*r*rand,0.00001);
    end
    end
    elseif r>=10^(-4)
    for j=1:d
        for k=1:d
            p(j,k) = max(q(j,k)+0.002*r*rand,0);
        end
    end
    else
        for j=1:d
        for k=1:d
            p(j,k) = max(q(j,k)+r*rand,0.00000000000001);
        end
        end
    end

    p = p./sum(p,2);
    Dc = 0;
    for i=1:d
        sum0=0;
        for j=1:d
            if (q(i,j)~=0)&&(p(i,j)~=0)
                 sum0 = sum0 + q(i,j)*log(q(i,j)/p(i,j));
            end
        end
        Dc = Dc + alpha0(i)*sum0;
    end
    disp(Dc);
    if Dc <= r
        val = p; 
        break
    end
    if isnan(Dc)
%         dbstop if naninf 
        disp('p');
        disp(p);
%         return
    end
end

    if isnan(Dc)
%         dbstop if naninf 
        disp('p');
        disp(p);
        error('no right sample');
    end
end
