%Solving for equilibrium in Aiyagari (1994)

alpha = 0.36;
beta = 0.99;
delta = 0.025;
theta = 2.0;
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
v2 = zeros(151,7); optk2 = zeros(151,7); optc2 = zeros(151,7); opth2 = zeros(151,7);
r = fzero(@(r) aiyagari_projection_egm_7_labor(r,alpha,beta,delta,theta,zgrid,piz,kgrid,f1),[0.032 0.035],options)
[f1,aggk,aggn,v2,optk2,optc2,opth2] = aiyagari_projection_egm2_7_labor(r,alpha,beta,delta,theta,zgrid,piz,kgrid,f1,v2,optk2,optc2,opth2,0);
hstar = (1-alpha)/(theta+1-alpha-delta*theta*(1/beta-1+delta)/alpha);
kstar = (alpha*beta/(1-beta*(1-delta)))^(1/(1-alpha))*hstar;
[aggk kstar]
save('aiyagari_dec_labor.mat','v2','optk2','optc2','opth2')
save('aiyagari_r_labor.mat','r','aggk','aggn')
save('aiyagari_f1_labor.mat','f1')
fstart = f1; vs = v2;

%construct figure of equilibrium r
% rpts = linspace(0.032,0.035,101); supply = zeros(101,1); demand = zeros(101,1);
% for i=1:101
%     [f1,aggk,aggn,v2,optk2,optc2,opth2] = aiyagari_projection_egm2_7_labor(rpts(i),alpha,beta,delta,theta,zgrid,piz,kgrid,f1,v2,optk2,optc2,opth2,1);
%     supply(i) = aggk/aggn;
%     demand(i) = (rpts(i)/alpha)^(1/(alpha-1));
% end
% plot(supply,rpts,'LineWidth',2)
% hold
% plot(demand,rpts,'LineWidth',2)
% plot(supply,ones(101,1)*(1/beta-1+delta),'k--','LineWidth',2)
% save('aiyagari_eq_labor.mat','rpts','supply','demand')

%compute dynamic equilibrium transition
% [f1,aggk,aggn,v,optk,optc,opth] = aiyagari_projection_egm2_7_labor(0.032,alpha,beta,delta,theta,zgrid,piz,kgrid,f1,v2,optk2,optc2,opth2,0);
% rvals2 = linspace(0.032,r,700);
% factor = 0.9;
% for iter = 1:100
%     rvals = rvals2;
%     rvals2 = aiyagari_projection_egm_7_transition_labor(rvals,alpha,beta,delta,theta,zgrid,piz,kgrid,f1,700,v2,optk2,optc2,opth2);
%     rdiff = norm(rvals2-rvals)
%     if rdiff < 1e-07
%         break
%     end
%     rvals2 = factor*rvals+(1-factor)*rvals2;
% end
% save('aiyagari_trans_labor.mat','rvals2')

%compute stochastic equilibrium
a1 = [(1-0.964874602777629)*log(aggk)*1.01 0.964874602777629]; 
b1 = [(1-0.964874602777629)*log(aggk)*0.99 0.964874602777629];
c1 = [log(aggn) 0];
d1 = [log(aggn) 0];
us = rand(10000,1);
Kgrid = linspace(0.9*aggk,1.1*aggk,10);
Agrid = [1.01;0.99]; pia = [0.95 0.05;0.05 0.95];
Ngrid = linspace(0.9*aggn,1.1*aggn,40);
factor = 0.9;
%v2 = zeros(151,7,10,2); 
%optka = zeros(151,7,10,2,40); optca = zeros(151,7,10,2,40); optha = zeros(151,7,10,2,40);
%for i = 1:10
%    for j = 1:2
%       v2(:,:,i,j) = vs;
%    end
%end
load aiyagari_agg_labor.mat
load aiyagari_law_labor.mat
vdiff = 100;
for i=1:9999
    a0 = a1; b0 = b1; c0 = c1; d0 = d1;
    v = v2;
    if vdiff > 1e-03
        sim = 0;
    else
        sim = 0;
    end
    [a1,b1,c1,d1,v2,optka,optca,optha,aggks,aggns]=aiyagari_projection_egm_7_dynamic_labor(a0,b0,c0,d0,alpha,beta,delta,theta,zgrid,piz,...
                                      kgrid,fstart,10000,Agrid,pia,Kgrid,Ngrid,aggk,us,v,sim);
    [a0 a1]
    [b0 b1]
    a1 = factor*a0+(1-factor)*a1;
    b1 = factor*b0+(1-factor)*b1;
    vdiff = norm(v2-v,"fro")
    save('aiyagari_tlaw.mat','a0','b0','a1','b1')
    if norm([a1-a0';b1-b0']) < 1e-05 && vdiff < 1e-07
        break
    end
end
save('aiyagari_law_labor.mat','a1','b1','c1','d1')
save('aiyagari_agg_labor.mat','v2','optka','optca','optha')
save('aiyagari_sim_labor.mat','aggks','aggns')