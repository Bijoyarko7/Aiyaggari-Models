function y = aiyagari_projection(r,alpha,beta,delta,zgrid,piz,m,kgrid)

w = (1-alpha)*(r/alpha)^(alpha/(alpha-1));
%aggn = 0.5*zgrid(1)+0.5*zgrid(2); 
options = optimoptions('fsolve','Display','none');

x0(1:m) = kgrid;
x0(m+1:2*m) = kgrid;
x0(2*m+1:3*m) = -sqrt(kgrid);
x0(3*m+1:4*m) = -sqrt(kgrid);
[x,f]=fsolve(@(x) aiyagari_pchip(x,alpha,beta,delta,r,w,kgrid,zgrid,piz),x0,options);
optk2(1:m,1) = x(1:m);
optk2(1:m,2) = x(m+1:2*m);

kgrid2=linspace(kgrid(1),kgrid(end),500); optkg=zeros(500,2);
for i=1:500
    jk=binarySearch(kgrid,kgrid2(i));
    jk=max(1,min(jk,length(kgrid)-1));
    ck=(kgrid2(i)-kgrid(jk))/(kgrid(jk+1)-kgrid(jk));
    optkg(i,:)=(1-ck)*optk2(jk,:)+ck*optk2(jk+1,:);
end
f1=zeros(500,2);
f1(1:100,1:2)=0.5/100;
[f1,aggk,aggn]=get_tmat(kgrid2,zgrid,piz,optkg,f1);

 figure(1)
 plot(kgrid2,f1(:,1),'LineWidth',2);
 hold on
 plot(kgrid2,f1(:,2),'LineWidth',2);
 hold off
 figure(2)
 plot(kgrid(1:20),optk2(1:20,1),'LineWidth',2)
 hold on
 plot(kgrid(1:20),optk2(1:20,2),'LineWidth',2)
 plot(kgrid(1:20),kgrid(1:20),'k--','LineWidth',2)
 hold off

y = alpha*aggk^(alpha-1)*aggn^(1-alpha)-r;
[r aggk y]
