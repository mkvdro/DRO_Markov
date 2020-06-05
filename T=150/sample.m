function val = sample(alpha0,r,q,d)

p=q;
N_samples = 10^7;
for n=1:N_samples
    fprintf('sample order %d',n);
    % randomly sample p
    if r>=0.8
%         p=rand(d);
        for j=1:d
            for k=1:d
            p(j,k) = max(q(j,k)+r*rand-0.5*r,0.000001); 
            end
        end
    elseif r>=0.1
        for j=1:d
        for k=1:d
            p(j,k) = max(q(j,k)+0.2*r*rand-0.2*r,10^(-15)); 
        end
        end
    elseif r>=0.07
    for j=1:d
    for k=1:d
        p(j,k) = max(q(j,k)+0.2*r*rand,10^(-15));
    end
    end
    elseif r>=10^(-4)
    for j=1:d
        for k=1:d
            p(j,k) = max(q(j,k)+0.25*r*rand,10^(-15));
        end
    end
    else
        for j=1:d
        for k=1:d
            p(j,k) = max(q(j,k)+.5*r*rand,0.00000000000001);
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
