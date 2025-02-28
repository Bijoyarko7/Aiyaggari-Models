function y = aiyagari_projection_egm_7_labor(r,alpha,beta,delta,theta,zgrid,piz,kgrid,f1)

w = (1-alpha)*(r/alpha)^(alpha/(alpha-1));

v2 = zeros(151,7); optk2 = zeros(151,7); optc2 = zeros(151,7); opth2 = zeros(151,7);
for i=1:7
    v2(1:151,i)=log((r-delta)*kgrid+w*zgrid(i)*0.33)/(1-beta);
    optk2(1:151,i)=kgrid;
    optc2(1:151,i)=(r-delta)*kgrid+w*zgrid(i)*0.33;
    opth2(1:151,i)=0.33;
end
chat = zeros(151,7); khat = zeros(151,7); hhat = zeros(151,7);

for iter = 1:9999
    v = v2;
    optk = optk2;
    optc = optc2;
    opth = opth2;
    for j = 1:7
        ev = zeros(151,1);
        for k = 1:7
            ev=ev+piz(j,k)*v(:,k);
        end
        pvv = pchip(kgrid,ev);
        [breaks,coefs,l,k,d] = unmkpp(pvv);
        pvvp = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
        chat(:,j) = 1/beta./ppval(pvvp,kgrid);
        hhat(:,j) = max(0,1-theta*chat(:,j)/(w*zgrid(j)));
        khat(:,j) = (chat(:,j)+kgrid-w*zgrid(j)*hhat(:,j))/(1+r-delta);
        pc = pchip(khat(:,j),chat(:,j));
        for i = 1:151
            if kgrid(i) > khat(1,j)
                optc2(i,j) = ppval(pc,kgrid(i));
                opth2(i,j) = max(0,1-theta*optc2(i,j)/(w*zgrid(j)));
                optk2(i,j) = (1+r-delta)*kgrid(i)+w*zgrid(j)*opth2(i,j)-optc2(i,j);
            else
                optk2(i,j) = 0.0;
                optc2(i,j) = ((1+r-delta)*kgrid(i)+w*zgrid(j))/(1+theta);
                opth2(i,j) = (optc2(i,j)-(1+r-delta)*kgrid(i))/(w*zgrid(j));
                if opth2(i,j)<0
                    opth2(i,j) = 0;
                    optc2(i,j) = (1+r-delta)*kgrid(i);
                end
            end
        end
        v2(:,j) = log(optc2(:,j))+theta*log(1-opth2(:,j))+beta*ppval(pvv,optk2(:,j));
    end
    vdiff = norm(v2-v);
    if vdiff < 1e-08
        break
    end
end
kgrid2=linspace(kgrid(1),400,5000); optkg=zeros(5000,7); opthg=zeros(5000,7);
for i=1:7
    optkg(:,i) = interp1(kgrid,optk2(:,i),kgrid2);
    opthg(:,i) = interp1(kgrid,opth2(:,i),kgrid2);
end
[f1,aggk,aggn]=get_tmat_labor(kgrid2,zgrid,piz,optkg,opthg,f1);

y = alpha*aggk^(alpha-1)*aggn^(1-alpha)-r;
[r aggk aggn y]