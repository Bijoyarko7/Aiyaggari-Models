%Solving for equilibrium in Aiyagari (1994)

alpha = 0.36;
beta = 0.99;
delta = 0.025;
[zgrid,piz] = rouwenhorst(0.95,0.1,7,0.0);
zgrid = exp(zgrid);
pit = piz^1000;
pinv = pit(1,:);
kgrid = zeros(151,1);
kgrid(1:50) = linspace(0,1,50);
kgrid(51:100) = linspace(1.25,10,50);
kgrid(101:130) = linspace(12,25,30);
kgrid(131:151) = linspace(28,500,21);
f1 = zeros(5000,7);
for i=1:7
    f1(1:1000,i)=pinv(i)/1000;
end

%solve for equilibrium steady state r
options = optimset('fzero');
options = optimset('Display','off','TolX',1e-08);
v2 = zeros(151,7); optk2 = zeros(151,7); optc2 = zeros(151,7);
r = fzero(@(r) aiyagari_projection_egm_7(r,alpha,beta,delta,zgrid,piz,kgrid,f1),[0.032 0.035],options)
[f1,aggk,aggn,v2,optk2,optc2] = aiyagari_projection_egm2_7(r,alpha,beta,delta,zgrid,piz,kgrid,f1,v2,optk2,optc2,0);
kstar = (alpha*beta/(1-beta*(1-delta)))^(1/(1-alpha))*aggn;
[aggk kstar]
save('aiyagari_dec.mat','v2','optk2','optc2')
save('aiyagari_r.mat','r','aggk','aggn')
save('aiyagari_f1.mat','f1')

%construct figure of equilibrium r
% rpts = linspace(0.032,0.035,101); supply = zeros(101,1); demand = zeros(101,1);
% for i=1:101
%     [f1,aggk,aggn,v2,optk2,optc2] = aiyagari_projection_egm2_7(rpts(i),alpha,beta,delta,zgrid,piz,kgrid,f1,v2,optk2,optc2,1);
%     supply(i) = aggk;
%     demand(i) = (rpts(i)/alpha)^(1/(alpha-1))*aggn;
% end
% plot(supply,rpts,'LineWidth',2)
% hold
% plot(demand,rpts,'LineWidth',2)
% plot(supply,ones(101,1)*(1/beta-1+delta),'k--','LineWidth',2)
% save('aiyagari_eq.mat','rpts','supply','demand')

%compute dynamic equilibrium transition
% [f1,aggk,aggn,v,optk,optc] = aiyagari_projection_egm2_7(0.032,alpha,beta,delta,zgrid,piz,kgrid,f1,v2,optk2,optc2,0);
% rvals2 = linspace(0.032,r,700);
% factor = 0.9;
% for iter = 1:100
%     rvals = rvals2;
%     rvals2 = aiyagari_projection_egm_7_transition(rvals,alpha,beta,delta,zgrid,piz,kgrid,f1,700,v2,optk2,optc2);
%     rdiff = norm(rvals2-rvals)
%     if rdiff < 1e-07
%         break
%     end
%     rvals2 = factor*rvals+(1-factor)*rvals2;
% end
% save('aiyagari_trans.mat','rvals2')

%compute stochastic equilibrium
a1 = [(1-0.95)*log(aggk) 0.95]; b1 = a1;
us = rand(10000,1);
Kgrid = linspace(0.9*aggk,1.1*aggk,10);
Agrid = [1.01;0.99]; pia = [0.95 0.05;0.05 0.95];
factor = 0.9;
hot = 0;
va = zeros(151,7,10,2); optka = zeros(151,7,10,2); optca = zeros(151,7,10,2);
for i=1:1000
    a0 = a1; b0 = b1;
    [a1,b1,va,optka,optca]=aiyagari_projection_egm_7_dynamic(a0,b0,alpha,beta,delta,zgrid,piz,kgrid,f1,10000,...
        v2,optk2,optc2,Agrid,pia,Kgrid,aggk,aggn,us,hot,va,optka,optca);
    hot = 1;
    if max([abs(a1-a0),abs(b1-b0)]) < 1e-06
        break
    end
    [a0' a1;b0' b1]
    a1 = factor*a0+(1-factor)*a1';
    b1 = factor*b0+(1-factor)*b1';
end
save('aiyagari_law.mat','a1','b1')
save('aiyagari_agg.mat','va','optka','optca')