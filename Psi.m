function val = Psi(x,p)

mc=dtmc(p);
pi0=asymptotics(mc);

val=-pi0*x;
if isnan(val)
    disp(pi0);
    error('Psi undef');
end
end