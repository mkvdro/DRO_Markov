rng(2);
% parameter setting

k=3; % how many customer segements: i.e., how many different markov chain dynamics
d=10; % how many brands
T=500; % length of each xi^(i)/ sample size
n_exper = 100; % number of independent experiments

P=rand(d,d,k); % k  transition matrix, each of dimension d*d
for i=1:k
    P(randi([1,d]),randi([1,d]),i)=5;%randi([1,d])
    P(randi([1,d]),randi([1,d]),i)=4;
    % P(randi([1,d]),randi([1,d]),i)=3;
    P(:,:,i)=P(:,:,i)./sum(P(:,:,i),2); %normalize to standard row-stochastic matrix
    mc=dtmc(P(:,:,i));
    while (isergodic(mc)==0) && (isreducible(mc)==1) % keep rolling the dice until getting an ergodic and irreducible chain dynamics
        P(:,:,i)=rand(d,d);
        P(:,:,i)=P(:,:,i)./sum(P(:,:,i),2); %normalize to standard row-stochastic matrix
        mc=dtmc(P(:,:,i));
    end
end


for n=1:n_exper
    % generating pseudo-data from predefined markov chain dynamics: Random Walk

    xi = ones(k,T); % for each customer segment, generate a different trajectory with regard to its jumping dynamics
    for i=1:k
        X=simulate(dtmc(P(:,:,i)),T-1);
        % disp(X);
        xi(i,:)=X';
    end
    xi=reshape(xi.',1,[]); %flatten xi matrix

    if n==1
        save('xi_new.mat','xi');
        matObj = matfile('xi_new.mat','Writable',true); %save
    else
        matObj.xi(n,:)=xi; % otherwise append
    end
end
save('transition_dynamics.mat','P')  % function form
