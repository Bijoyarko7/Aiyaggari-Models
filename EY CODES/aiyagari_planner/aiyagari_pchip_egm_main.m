%Solving for equilibrium in Aiyagari (1994)

tic;
alpha = 0.36;
beta = 0.96;
delta = 0.08;
zgrid = [1.0 0.01];
piz = [0.95 0.05;0.5 0.5];
pit = piz^1000;
pinv = pit(1,:);

m = 151; n = length(pinv);
kgrid(1:50) = linspace(0,2,50);
kgrid(51:100) = linspace(2.5,25,50);
kgrid(101:150) = linspace(27,1000,50);
kgrid(151) = 2500;
options = optimset('Display','none');

r = fzero(@(r) aiyagari_egm(r,alpha,beta,delta,zgrid,piz,m,n,kgrid,pinv),[0.1182 0.1184],options)
[aggk2,aggn2,optks,optcs,f1] = aiyagari_egm2(r,alpha,beta,delta,zgrid,piz,m,n,kgrid,pinv);
save('aiyagari_unemp.t1','optks','optcs','f1','r','aggk2','aggn2')

ng = 50000;
optkp = optks; optcp = optcs; jks = zeros(ng,1); cks = zeros(ng,1);
kgrid2=linspace(kgrid(1),1500,ng);
for i=1:ng
    jk=binarySearch(kgrid,kgrid2(i));
    jks(i)=max(1,min(jk,length(kgrid)-1));
    cks(i)=(kgrid2(i)-kgrid(jks(i)))/(kgrid(jks(i)+1)-kgrid(jks(i)));
end
g2 = -0.005235824355089;
r2 = 0.124727866226197;
for t = 1:1000
    r = r2;
    g = g2;
    [optkp,optcp,fx,f1,r2,g2,aggk2,aggn2]=aiyagari_egm_planner2(r,g,alpha,beta,delta,zgrid,piz,m,n,kgrid,optkp,optcp,pinv,f1,ng,kgrid2,cks,jks);
    if abs(g2-g)<1e-07 && abs(r2-r)<1e-07
        break
    else
        [r g abs(r2-r) abs(g2-g)]
        r2 = 0.99*r+0.01*r2;
        g2 = 0.99*g+0.01*g2;
    end
end
save('aiyagari_unemp.t2','optkp','optcp','f1','r','g','aggk2','aggn2')
