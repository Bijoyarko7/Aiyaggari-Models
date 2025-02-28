function res = aiyagari_pchip(x,alpha,beta,delta,r,w,kgrid,zgrid,piz)

m = length(kgrid);
pp1 = pchip(kgrid,x(1:m));
pp2 = pchip(kgrid,x(m+1:2*m));
res = zeros(length(x),1);
for i = 1:m
    kp1 = x(i);
    c = (1+r-delta)*kgrid(i)+w*zgrid(1)-kp1;
    kpp1 = ppval(pp1,kp1);
    cp1 = (1+r-delta)*kp1+w*zgrid(1)-kpp1;
    kpp2 = ppval(pp2,kp1);
    cp2 = (1+r-delta)*kp1+w*zgrid(2)-kpp2;
    mu1 = x(i+2*m);
    res(i) = c^(-1)-max(mu1,0)^2-beta*(1+r-delta)*(piz(1,1)*cp1^(-1)+piz(1,2)*cp2^(-1));
    kp2 = x(i+m);
    c = (1+r-delta)*kgrid(i)+w*zgrid(2)-kp2;
    kpp1 = ppval(pp1,kp2);
    cp1 = (1+r-delta)*kp2+w*zgrid(1)-kpp1;
    kpp2 = ppval(pp2,kp2);
    cp2 = (1+r-delta)*kp2+w*zgrid(2)-kpp2;
    mu2 = x(i+3*m);
    res(i+m) = c^(-1)-max(mu2,0)^2-beta*(1+r-delta)*(piz(2,1)*cp1^(-1)+piz(2,2)*cp2^(-1));
    res(i+2*m) = max(-mu1,0)^2-kp1;
    res(i+3*m) = max(-mu2,0)^2-kp2;
end

