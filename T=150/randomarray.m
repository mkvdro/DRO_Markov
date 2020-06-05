d=20; % how many assets
m=10; % cap for values of asset
T=200; % length of each xi^(i)
% for n=1:10^7
xi=zeros(d,T);
% for i=1:d
%     p_i(i)=round(0.3+rand*0.6,2);
%     q_i(i)=round((1-p_i(d))*rand,2);
% end
% p_i=[0.86 0.74 0.43 0.32 0.88 0.6 0.67];
% q_i=[0.12 0.05 0.12 0.03 0.12 0.22 0.31];
% p_i=[0.86 0.74 0.33 0.52 0.88 0.6 0.67];
% q_i=[0.12 0.05 0.12 0.03 0.12 0.22 0.31];
p_i=zeros(1,d);
q_i=zeros(1,d);
for i=1:d
    p_i(i)=round(0.3+rand*0.6,2);
    q_i(i)=round((1-p_i(i))*rand,2);
end
p=p_i;q=q_i;