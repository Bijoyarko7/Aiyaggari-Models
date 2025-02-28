function [a1,b1,v,optk,optc] = aiyagari_projection_egm_7_dynamic(a0,b0,alpha,beta,delta,zgrid,piz,...
    kgrid,fstart,capt,vs,optks,optcs,Agrid,pia,Kgrid,aggk,aggn,us,hot,v2,optk2,optc2)

r = zeros(10,2); w = zeros(10,2);
for l=1:2
    r(:,l)=alpha*Agrid(l)*Kgrid.^(alpha-1)*aggn^(1-alpha);
    w(:,l)=(1-alpha)*Agrid(l)*Kgrid.^(alpha)*aggn^(-alpha);
end

if hot == 0
    v2 = zeros(151,7,10,2); optk2 = zeros(151,7,10,2); optc2 = zeros(151,7,10,2);
    for k=1:10
        for l=1:2
            v2(1:151,1:7,k,l)=vs;
            optk2(1:151,1:7,k,l)=optks;
            optc2(1:151,1:7,k,l)=optcs;
        end
    end
end
chat = zeros(151,7); khat = zeros(151,7); Kp = zeros(10,2); vstar = zeros(151,7,10,2,2);
%vstar is organized as k,z,K,a,a'

for iter = 1:9999
    v = v2; optk = optk2; optc = optc2;
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
                khat(:,j) = (chat(:,j)+kgrid-w(k,l)*zgrid(j))/(1+r(k,l)-delta);
                pc = pchip(khat(:,j),chat(:,j));
                for i = 1:151
                    if kgrid(i) > khat(1,j)
                        optc2(i,j,k,l) = ppval(pc,kgrid(i));
                        optk2(i,j,k,l) = (1+r(k,l)-delta)*kgrid(i)+w(k,l)*zgrid(j)-optc2(i,j,k,l);
                    else
                        optk2(i,j,k,l) = 0.0;
                        optc2(i,j,k,l) = (1+r(k,l)-delta)*kgrid(i)+w(k,l)*zgrid(j);
                    end
                end
                v2(:,j,k,l) = log(optc2(:,j,k,l))+beta*ppval(pvv,optk2(:,j,k,l));
            end
        end
    end
    vdiff = norm(abs(v2-v),"fro");
    if vdiff < 1e-07
        break
    end
end
kgrid2=linspace(kgrid(1),400,5000); optkg=zeros(5000,7,10,2);
for i=1:7
    for j=1:10
        for k=1:2
            optkg(:,i,j,k) = interp1(kgrid,optk2(:,i,j,k),kgrid2);
        end
    end
end
f1 = zeros(5000,7,capt+1);
f1(:,:,1) = fstart;
aggks = zeros(capt+1,1); nas = zeros(capt,1);
nas(1) = 1;
for t=1:capt
    if us(t)<pia(nas(t),1)
        nas(t+1)=1;
    else
        nas(t+1)=2;
    end
end
aggks(1)=aggk;
for t=1:capt
    jK = binarySearch(Kgrid,aggks(t));
    jA = nas(t);
    cK = (aggks(t)-Kgrid(jK))/(Kgrid(jK+1)-Kgrid(jK));
    optkgg = (1-cK)*optkg(:,:,jK,jA)+cK*optkg(:,:,jK+1,jA);
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
X1 = zeros(n1,2); Y1 = zeros(n1,1);
X2 = zeros(n2,2); Y2 = zeros(n2,1);
n = 0;
for t = 1:capt
    if nas(t) == 1
        n = n+1;
        X1(n,1)=1;
        X1(n,2)=log(aggks(t));
        Y1(n,1)=log(aggks(t+1));
    end
end
a1 = (X1'*X1)\(X1'*Y1);
n = 0;
for t = 1:capt
    if nas(t) == 2
        n = n+1;
        X2(n,1)=1;
        X2(n,2)=log(aggks(t));
        Y2(n,1)=log(aggks(t+1));
    end
end
b1 = (X2'*X2)\(X2'*Y2);
