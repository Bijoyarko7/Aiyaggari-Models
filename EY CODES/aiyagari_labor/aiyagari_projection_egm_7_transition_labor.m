function y = aiyagari_projection_egm_7_transition_labor(r,alpha,beta,delta,theta,zgrid,piz,kgrid,fstart,capt,v,optk,optc,opth)

w = (1-alpha)*(r/alpha).^(alpha/(alpha-1));
y = zeros(1,capt);

v2 = zeros(151,7,capt+1); optk2 = zeros(151,7,capt+1); optc2 = zeros(151,7,capt+1); opth2 = zeros(151,7,capt+1);
v2(1:151,1:7,capt+1)=v;
optk2(1:151,1:7,capt+1)=optk;
optc2(1:151,1:7,capt+1)=optc;
opth2(1:151,1:7,capt+1)=opth;
chat = zeros(151,7); khat = zeros(151,7); hhat = zeros(151,7);

for iter = capt:-1:1
    for j = 1:7
        ev = zeros(151,1);
        for k = 1:7
            ev=ev+piz(j,k)*v2(:,k,iter+1);
        end
        pvv = pchip(kgrid,ev);
        [breaks,coefs,l,k,d] = unmkpp(pvv);
        pvvp = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
        chat(:,j) = 1/beta./ppval(pvvp,kgrid);
        hhat(:,j) = max(0,1-theta*chat(:,j)/(w(iter)*zgrid(j)));
        khat(:,j) = (chat(:,j)+kgrid-w(iter)*zgrid(j)*hhat(:,j))/(1+r(iter)-delta);
        pc = pchip(khat(:,j),chat(:,j));
        for i = 1:151
            if kgrid(i) > khat(1,j)
                optc2(i,j,iter) = ppval(pc,kgrid(i));
                opth2(i,j,iter) = max(0,1-theta*optc2(i,j,iter)/(w(iter)*zgrid(j)));
                optk2(i,j,iter) = (1+r(iter)-delta)*kgrid(i)+w(iter)*zgrid(j)*opth2(i,j,iter)-optc2(i,j,iter);
            else
                optk2(i,j,iter) = 0.0;
                optc2(i,j,iter) = ((1+r(iter)-delta)*kgrid(i)+w(iter)*zgrid(j))/(1+theta);
                opth2(i,j,iter) = (optc2(i,j,iter)-(1+r(iter)-delta)*kgrid(i))/(w(iter)*zgrid(j));
                if opth2(i,j,iter)<0
                    opth2(i,j,iter) = 0;
                    optc2(i,j,iter) = (1+r(iter)-delta)*kgrid(i);
                end
            end
        end
        v2(:,j,iter) = log(optc2(:,j,iter))+theta*log(1-opth2(:,j,iter))+beta*ppval(pvv,optk2(:,j,iter));
    end
end
kgrid2=linspace(kgrid(1),400,5000); optkg=zeros(5000,7,capt); opthg=zeros(5000,7,capt);
for i=1:7
    for j=1:capt
        optkg(:,i,j) = interp1(kgrid,optk2(:,i,j),kgrid2);
        opthg(:,i,j) = interp1(kgrid,opth2(:,i,j),kgrid2);
    end
end
f1 = zeros(5000,7,capt+1);
f1(:,:,1) = fstart;
for t=1:capt
    [f1(:,:,t+1),aggk,aggn]=get_tmat_transition_labor(kgrid2,zgrid,piz,optkg(:,:,t),opthg(:,:,t),f1(:,:,t));
    y(t) = alpha*aggk^(alpha-1)*aggn^(1-alpha)-r(t);
end
