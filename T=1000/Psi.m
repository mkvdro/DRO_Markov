function val = Psi(d,p)
%p input vector
pm=zeros(d);

for i=1:d
    for j=1:d
        pm(i,j)=p(d*(i-1)+j);
    end
end

I = eye(size(pm));

Y = null(pm'-I);
pi0 = Y./(sum(Y));
% [V, ~] = eig( pm.' ); %// note the transpose .' - we are looking for the **left** EV
% st = V(:,1).';
% pi0=st./sum(st);

% L=zeros(1,d);
val=0;
for j=1:d
   val=val+(j-1)*pi0(j); %L(j)=h*eta*j
end
% disp('pi0');disp(pi0);
% val=dot(L,pi0);
if isnan(val)
    disp(pi0);
    error('Psi undef');
end
end