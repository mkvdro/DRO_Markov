clear all
close all;
% system('caffeinate &');
warning('off');
tic;
% reliability plot
d=5; % how many assets
m=10; % cap for values of asset
T=150; % length of each xi^(i)
N_sample=10; 

p_i=[0.86 0.74 0.43 0.32 0.88];
q_i=[0.12 0.05 0.02 0.03 0.12];

    
[true_x,hold_real,cost_real] = test_real(d,m,p_i,q_i);

matObj = matfile('xi_new.mat','Writable',true);
[nrows, ~] = size(matObj, 'xi');
xi=zeros(d,T,N_sample);
for ns=1:N_sample
    xi(:,:,ns)=(reshape(matObj.xi(ns,:),[T,d])).';
end

epsilon=0.001;
iter=10;

cost_fin=zeros(1,N_sample); 
cost_out=zeros(1,N_sample); 
% 
cost_fin_iid=zeros(1,N_sample);
cost_out_iid=zeros(1,N_sample);
% 
cost_fin_saa=zeros(1,N_sample);
cost_out_saa=zeros(1,N_sample);

%%%%
xr=logspace(-4,0,10); %range r
[~,N_r]=size(xr);
y=zeros(1,N_r);
y1=zeros(1,N_r);
y2=zeros(1,N_r);
reliability=zeros(1,N_r);
reliability_iid=zeros(1,N_r);

z=zeros(1,N_r);
z1=zeros(1,N_r);
z2=zeros(1,N_r);

w=zeros(1,N_r);
w1=zeros(1,N_r);
w2=zeros(1,N_r);

i=1;
for r=xr
    fprintf('order r %d ',i);
%     disp(r);
    n_disappt=0;
    
for nsample=1:N_sample
    fprintf('exp order %d ',nsample);
    [alpha0,q]=est_alpha_from_xi(d,m,T,xi(:,:,nsample));
    hold_cost_m=zeros(1,d);
    c_exp=zeros(1,d);
    for j=1:d
        hold_cost_m(j)=FW_main(epsilon,r,iter,q(:,:,j),m,alpha0(:,:,j));
        c_exp(j)=0.2;
    end
    opts = optimoptions('quadprog','Display','off');
    [x,cost_fin(nsample)] = quadprog(eye(d)*(1/c_exp(1)),hold_cost_m,[],[],ones(1,d),1,zeros(1,d),ones(1,d),[],opts);
    cost_out(nsample) = hold_real*x+1/2*x'*eye(d)*(1/c_exp(1))*x;
    
    if cost_out(nsample)>cost_fin(nsample)
        n_disappt=n_disappt+1;
    end
end
reliability(i)=1-n_disappt/N_sample;

z(i)=mean(cost_out);
z1(i)=z(i)-2*std(cost_out);
z2(i)=z(i)+2*std(cost_out);
i=i+1;
end

i=1;
for r=xr
    n_disappt=0;
for nsample=1:N_sample
    alpha0=naive_est_alpha(d,m,T,xi(:,:,nsample));
    hold_cost_i=zeros(1,d);
    for j=1:d
        hold_cost_i(j)=cost_noM(alpha0(j,:),r,m);
    end
    
    opts = optimoptions('quadprog','Display','off');
    [xiid,cost_fin_iid(nsample)]=quadprog(eye(d)*(1/c_exp(1)),hold_cost_i,[],[],ones(1,d),1,zeros(1,d),ones(1,d),[],opts);
    cost_out_iid(nsample) = hold_real*xiid+1/2*xiid'*eye(d)*(1/c_exp(1))*xiid;
    
    if cost_out_iid(nsample)>cost_fin_iid(nsample)
        n_disappt=n_disappt+1;
    end
end
reliability_iid(i)=1-n_disappt/N_sample;


w(i)=mean(cost_out_iid); 
w1(i)=w(i)-2*std(cost_out_iid);
w2(i)=w(i)+2*std(cost_out_iid);
i=i+1;
end


n_disappt=0;
for nsample=1:N_sample
    alpha0=naive_est_alpha(d,m,T,xi(:,:,nsample));
    hold_cost_saa=zeros(1,d);
    for j=1:d
        hold_cost_saa(j)=dot(0:(m-1),alpha0(j,:));
    end
    
    opts = optimoptions('quadprog','Display','off');
    [xiid,cost_fin_saa(nsample)]=quadprog(eye(d)*(1/c_exp(1)),hold_cost_saa,[],[],ones(1,d),1,zeros(1,d),ones(1,d),[],opts);
    cost_out_saa(nsample) = hold_real*xiid+1/2*xiid'*eye(d)*(1/c_exp(1))*xiid;
    
    if cost_out_saa(nsample)>cost_fin_saa(nsample)
        n_disappt=n_disappt+1;
    end
    
end
reliability_SAA=1-n_disappt/N_sample+0*xr;
saa_out=mean(cost_out_saa)+0*xr;
saa_2=mean(cost_out_saa)-2*std(cost_out_saa)+0*xr;
saa_3=mean(cost_out_saa)+2*std(cost_out_saa)+0*xr;


figure(1)
hold on;
hmeansaa=plot(xr, reliability_SAA, 'LineWidth',2);
hmeansaa.Color='y';
hmeaniid=plot(xr, reliability_iid, 'LineWidth',2);
hmeaniid.Color='r';
hmean_m=plot(xr, reliability,'LineWidth',2);
hmean_m.Color='b';
set(gca,'XScale','log')

hold off;

figure(2)
hold on;
%markov
x3 = [xr, fliplr(xr)];
inBetween = [z1, fliplr(z2)];
h2=fill(x3, inBetween, 'b','Edgecolor', 'none');
set(h2,'FaceAlpha',0.2)
hmeanout=plot(xr,z, 'b', 'LineWidth', 2);

%iid
x4 = [xr, fliplr(xr)];
inBetween2 = [w1, fliplr(w2)];
h3=fill(x4, inBetween2, 'r','Edgecolor', 'none');
set(h3,'FaceAlpha',0.2)
hmeanout_iid=plot(xr,w, 'r', 'LineWidth', 2);

real=plot(xr,cost_real+0*xr, 'g', 'LineWidth', 2);
xlabel('r')
ylabel('cost')
% legend([hline],{'True cost'})
set(gca,'XScale','log')
hold off;

figure(3)
hold on;
%markov
x3 = [xr, fliplr(xr)];
inBetween = [z1, fliplr(z2)];
h2=fill(x3, inBetween, 'b','Edgecolor', 'none');
set(h2,'FaceAlpha',0.2)
hmeanout=plot(xr,z, 'b', 'LineWidth', 2);

%saa
x4 = [xr, fliplr(xr)];
inBetween2 = [saa_2, fliplr(saa_3)];
h3=fill(x4, inBetween2, 'y','Edgecolor', 'none');
set(h3,'FaceAlpha',0.2)
hmeanout_saa=plot(xr,saa_out, 'y', 'LineWidth', 2);

real=plot(xr,cost_real+0*xr, 'g', 'LineWidth', 2);
xlabel('r')
ylabel('cost')
% legend([hline],{'True cost'})
set(gca,'XScale','log')


% FW_main.m
% grad_pi_k.m
% Grad_Psi.m
% linear_sub.m
% main.m
% minor.m
% naive_est_alpha.m
% Psi.m
% sample.m
% so_minor.m
% test_real.m

save('simulation_data')
t=toc;
etime=t./60
% matlab2tikz('T200.tex');

