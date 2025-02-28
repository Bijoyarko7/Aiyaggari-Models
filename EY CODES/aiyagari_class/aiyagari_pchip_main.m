%Solving for equilibrium in Aiyagari (1994)

tic;
alpha = 0.36;
beta = 0.99;
delta = 0.025;
zgrid = [1.25;0.75];
piz = [0.9 0.1;0.1 0.9];
m = 51;
kgrid(1:20) = linspace(0,2,20);
kgrid(21:40) = linspace(2.5,25,20);
kgrid(41:51) = linspace(30,200,11);
fstart = zeros(500,2); fstart(1:100,1)=0.5/100; fstart(1:100,2)=0.5/100;
options = optimoptions('fsolve','Display','none');

r = fzero(@(r) aiyagari_projection(r,alpha,beta,delta,zgrid,piz,m,kgrid),[0.0348 0.035])
xout(1:m) = kgrid; xout(m+1:2*m) = kgrid; xout(2*m+1:3*m) = -sqrt(kgrid); xout(3*m+1:4*m) = -sqrt(kgrid);
[aggk1,aggn1,~,fstart] = aiyagari_projection_2(0.034,alpha,beta,delta,zgrid,piz,m,kgrid,xout,fstart);
f1 = fstart;
[aggk2,aggn2,xouts,~] = aiyagari_projection_2(r,alpha,beta,delta,zgrid,piz,m,kgrid,xout,f1);
optks(:,1)=xouts(1:m); optks(:,2)=xouts(m+1:2*m);
toc

tic;
rs = linspace(0.025,0.035,101); supply = zeros(101,1); demand = zeros(101,1);
xout(1:m) = kgrid; xout(m+1:2*m) = kgrid; xout(2*m+1:3*m) = -sqrt(kgrid); xout(3*m+1:4*m) = -sqrt(kgrid);
f1=zeros(500,2);
f1(1:100,1:2)=0.5/100;
for i=1:101
    [aggk,aggn,xout,f1] = aiyagari_projection_2(rs(i),alpha,beta,delta,zgrid,piz,m,kgrid,xout,f1);
    supply(i) = aggk;
    demand(i) = (rs(i)/alpha)^(1/(alpha-1))*aggn;
end
figure(3)
plot(supply,rs,'LineWidth',2)
hold on
plot(demand,rs,'LineWidth',2)
plot(supply,(1/beta-1+delta)*ones(101,1),'k--','LineWidth',2)
toc
