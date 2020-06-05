function [x,sumco, val] = test_real(d,m,p_i,q_i)

P_real=zeros(m,m,d);

p=p_i;
q=q_i;

for i=1:d
    P_real(1,1,i)=1-p(i);
    P_real(1,2,i)=p(i);
    for k=2:m-1
        P_real(k,k-1,i)=q(i);
        P_real(k,k,i)=1-p(i)-q(i);
        P_real(k,k+1,i)=p(i);
    end
    P_real(m,m-1,i)=q(i);
    P_real(m,m,i)=1-q(i);
end

disp(P_real);

%calculating stationary distr.
alpha_real=zeros(m,1,d);
for i=1:d
I = eye(m);

Y = null(P_real(:,:,i)'-I);
pi_approx = Y./(sum(Y));
alpha_real(:,:,i)=pi_approx;
end

sumco=zeros(1,d);
c_exp=zeros(1,d);

for i=1:d
    sumco(i) = dot(0:m-1,alpha_real(:,:,i));
    c_exp(i)=0.2;
end

disp(sumco);
disp('sumco');

opts = optimoptions('quadprog','Display','off');
[x,val]=quadprog(eye(d)*(1/c_exp(1)),sumco,[],[],ones(1,d),1,zeros(1,d),ones(1,d),[],opts);
% val0=-sumco*x+1/2*x'*(eye(d).*(1/c_exp(1)))*x;
disp('true x');
disp(x);
end
