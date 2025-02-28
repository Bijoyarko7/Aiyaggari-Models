function [a1,b1,c1,d1,v2,optk,optc,opth,aggks,aggns] = aiyagari_projection_egm_7_dynamic_labor(a0,b0,c0,d0,alpha,beta,delta,theta,zgrid,piz,...
    kgrid,fstart,capt,Agrid,pia,Kgrid,Ngrid,aggk,us,v,sim)

r = zeros(10,2,40); w = zeros(10,2,40);
for l=1:2
    for n=1:40
        r(:,l,n)=alpha*Agrid(l)*Kgrid.^(alpha-1)*Ngrid(n)^(1-alpha);
        w(:,l,n)=(1-alpha)*Agrid(l)*Kgrid.^(alpha)*Ngrid(n)^(-alpha);
    end
end
v2 = zeros(151,7,10,2); optk = zeros(151,7,10,2,40); optc = zeros(151,7,10,2,40); opth = zeros(151,7,10,2,40); vhat = zeros(151,7,10,2,40);
chat = zeros(151,7); khat = zeros(151,7); hhat = zeros(151,7); Kp = zeros(10,2); vstar = zeros(151,7,10,2,2);
%vstar is organized as k',z',K,a,a'

    for k = 1:10
        Kp(k,1)=exp(a0(1)+a0(2)*log(Kgrid(k)));
        Kp(k,2)=exp(b0(1)+b0(2)*log(Kgrid(k)));
        jK = binarySearch(Kgrid,Kp(k,1));
        A = (Kp(k,1)-Kgrid(jK))/(Kgrid(jK+1)-Kgrid(jK));
        vstar(1:151,1:7,k,1,1) = (1-A)*v(1:151,1:7,jK,1)+A*v(1:151,1:7,jK+1,1);
        vstar(1:151,1:7,k,1,2) = (1-A)*v(1:151,1:7,jK,2)+A*v(1:151,1:7,jK+1,2);
        jK = binarySearch(Kgrid,Kp(k,2));
        A = (Kp(k,2)-Kgrid(jK))/(Kgrid(jK+1)-Kgrid(jK));
        vstar(1:151,1:7,k,2,1) = (1-A)*v(1:151,1:7,jK,1)+A*v(1:151,1:7,jK+1,1);
        vstar(1:151,1:7,k,2,2) = (1-A)*v(1:151,1:7,jK,2)+A*v(1:151,1:7,jK+1,2);
        for n = 1:40
            for l = 1:2
                for j = 1:7
                    ev = zeros(151,1);
                    for jj = 1:7
                        for ll = 1:2
                            ev=ev+piz(j,jj)*pia(l,ll)*vstar(:,jj,k,l,ll);
                        end
                    end
                    pvv = pchip(kgrid,ev);
                    [breaks,coefs,la,ka,da] = unmkpp(pvv);
                    pvvp = mkpp(breaks,repmat(ka-1:-1:1,da*la,1).*coefs(:,1:ka-1),da);
                    chat(:,j) = 1/beta./ppval(pvvp,kgrid);
                    hhat(:,j) = max(0,1-theta*chat(:,j)/(w(k,l,n)*zgrid(j)));
                    khat(:,j) = (chat(:,j)+kgrid-w(k,l,n)*zgrid(j)*hhat(:,j))/(1+r(k,l,n)-delta);
                    pc = pchip(khat(:,j),chat(:,j));
                    for i = 1:151
                        if kgrid(i) > khat(1,j)
                            optc(i,j,k,l,n) = ppval(pc,kgrid(i));
                            opth(i,j,k,l,n) = max(0,1-theta*optc(i,j,k,l,n)/(w(k,l,n)*zgrid(j)));
                            optk(i,j,k,l,n) = (1+r(k,l,n)-delta)*kgrid(i)+w(k,l,n)*zgrid(j)*opth(i,j,k,l,n)-optc(i,j,k,l,n);
                        else
                            optk(i,j,k,l,n) = 0.0;
                            optc(i,j,k,l,n) = ((1+r(k,l,n)-delta)*kgrid(i)+w(k,l,n)*zgrid(j))/(1+theta);
                            opth(i,j,k,l,n) = (optc(i,j,k,l,n)-(1+r(k,l,n)-delta)*kgrid(i))/(w(k,l,n)*zgrid(j));
                            if opth(i,j,k,l,n)<0
                                opth(i,j,k,l,n) = 0;
                                optc(i,j,k,l,n) = (1+r(k,l,n)-delta)*kgrid(i);
                            end
                        end
                    end
                    vhat(:,j,k,l,n) = log(optc(:,j,k,l,n))+theta*log(1-opth(:,j,k,l,n))+beta*ppval(pvv,optk(:,j,k,l,n));
                end
            end
        end
    end

kgrid2=linspace(kgrid(1),400,5000); optkg=zeros(5000,7,10,2,40); opthg=zeros(5000,7,10,2,40); 
for i=1:7
    for j=1:10
        for k=1:2
            for n=1:40
                optkg(:,i,j,k,n) = interp1(kgrid,optk(:,i,j,k,n),kgrid2);
                opthg(:,i,j,k,n) = interp1(kgrid,opth(:,i,j,k,n),kgrid2);
            end
        end
    end
end
f1 = zeros(5000,7,capt+1);
f1(:,:,1) = fstart;
aggks = zeros(capt+1,1); nas = zeros(capt,1); aggns = zeros(capt,1);
nas(1) = 1;
for t=1:capt
    if us(t)<pia(nas(t),1)
        nas(t+1)=1;
    else
        nas(t+1)=2;
    end
end
aggks(1)=aggk;
options = optimset('fzero');
options = optimset('Display','off','TolX',1e-08);
if sim == 1
for t=1:capt
    jK = binarySearch(Kgrid,aggks(t));
    jA = nas(t);
    cK = (aggks(t)-Kgrid(jK))/(Kgrid(jK+1)-Kgrid(jK));
    aggns(t) = fzero(@(n) getn(n,opthg,f1(:,:,t),jK,cK,jA,Ngrid,zgrid),[1.01*Ngrid(1) 0.99*Ngrid(end)],options);
    jN = binarySearch(Ngrid,aggns(t));
    cN = (aggns(t)-Ngrid(jN))/(Ngrid(jN+1)-Ngrid(jN));
    optkgg = (1-cK)*(1-cN)*optkg(:,:,jK,jA,jN)+cK*(1-cN)*optkg(:,:,jK+1,jA,jN)+...
        (1-cK)*cN*optkg(:,:,jK,jA,jN+1)+cK*cN*optkg(:,:,jK+1,jA,jN+1);
    [f1(:,:,t+1),aggks(t+1)]=get_tmat_dynamic(kgrid2,zgrid,piz,optkgg,f1(:,:,t));
end
n1 = 0; n2 = 0; 
for t = 1:capt
    if nas(t) == 1
        n1 = n1+1;
    else
        n2 = n2+1;
    end
end
X1 = zeros(n1,2); Y1 = zeros(n1,1); Z1 = zeros(n1,1);
X2 = zeros(n2,2); Y2 = zeros(n2,1); Z2 = zeros(n2,1);
n = 0;
for t = 1:capt
    if nas(t) == 1
        n = n+1;
        X1(n,1)=1;
        X1(n,2)=log(aggks(t));
        Y1(n,1)=log(aggks(t+1));
        Z1(n,1)=log(aggns(t));
    end
end
a1 = (X1'*X1)\(X1'*Y1);
c1 = (X1'*X1)\(X1'*Z1);
n = 0;
for t = 1:capt
    if nas(t) == 2
        n = n+1;
        X2(n,1)=1;
        X2(n,2)=log(aggks(t));
        Y2(n,1)=log(aggks(t+1));
        Z2(n,1)=log(aggns(t));
    end
end
b1 = (X2'*X2)\(X2'*Y2);
d1 = (X2'*X2)\(X2'*Z2);
a1 = a1'; b1 =  b1';
else
    a1 = a0;
    b1 = b0;
    c1 = c0;
    d1 = d0;
end
for k = 1:10
    aggn1 = exp(c1(1)+c1(2)*log(Kgrid(k)));
    jN = binarySearch(Ngrid,aggn1);
    cN = (aggn1-Ngrid(jN))/(Ngrid(jN+1)-Ngrid(jN));
    v2(:,:,k,1)=(1-cN)*vhat(:,:,k,1,jN)+cN*vhat(:,:,k,1,jN+1);
    aggn2 = exp(d1(1)+d1(2)*log(Kgrid(k)));
    jN = binarySearch(Ngrid,aggn2);
    cN = (aggn2-Ngrid(jN))/(Ngrid(jN+1)-Ngrid(jN));
    v2(:,:,k,2)=(1-cN)*vhat(:,:,k,2,jN)+cN*vhat(:,:,k,2,jN+1);
end

