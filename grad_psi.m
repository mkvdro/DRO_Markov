function val0 = Grad_Psi(L,P)
    % P input TRANSITION matrix
    % L: negative of the vector of constant that appears in the loss function: i.e., Psi= sum(L(x,i)*pi_i,i)

    % this is done using chain rule and formula provided in "sensitivity analysis of discrete markov chains" by Hal Caswell

    % first we compute d pi/ d P using the formula
    % to do that we use the fundamental matrix of ergodic MC: Z
    mc=dtmc(P);
    pi_0=asymptotics(mc); % stationary distribution (row vector)
    pi_0=pi_0';
    d = length(P); %size of the matrix
    Z=inv(eye(d)-P+pi_0*ones(1,d)); % fundamental matrix formula
    d_pi = kron(pi_0',(Z-pi_0*ones(1,d))); % the given formula
    val0=-L'*d_pi; % dot product following from the chain rule
        
end