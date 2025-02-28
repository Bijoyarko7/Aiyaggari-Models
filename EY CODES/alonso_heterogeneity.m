%planning problem with fixed heterogeneous beliefs

a2=1.05; a1=0.95;
A1=[a2;a2;a1;a1]; A2=[a2;a1;a2;a1];

pi = kron([0.75 0.25;0.25 0.75],[0.75 0.25;0.25 0.75]);
v1 = 0.05; v2 = 0.05;
pi1=pi+[v1 -v1 v1 -v1;v1 -v1 v1 -v1;v1 -v1 v1 -v1;v1 -v1 v1 -v1];
pi2=pi+[v2 v2 -v2 -v2;v2 v2 -v2 -v2;v2 v2 -v2 -v2;v2 v2 -v2 -v2];
beta = 0.96;

thetas = linspace(0.01,0.99,51);
pit = zeros(51,4,4);
for i = 1:51
    pit(i,:,:) = thetas(i)*pi1+(1-thetas(i))*pi2;
end
C2 = zeros(51,4);
for iter = 1:9999
    C = C2;
    for j = 1:4
        pC(j) = pchip(thetas,C(:,j));
    end
    for i = 1:51
        for j = 1:4
            ecsum = 0.0;
            for k = 1:4
                thetap = thetas(i)*pi1(j,k)/(thetas(i)*pi1(j,k)+(1-thetas(i))*pi2(j,k));
                termC = ppval(pC(k),thetap);
                kappa1 = beta*(a2*(pi1(2,k)*thetas(i)+pi2(2,k)*(1-thetas(i)))-a1*(pi1(3,k)*thetas(i)+pi2(3,k)*(1-thetas(i))))/...
                    ((a2-a1)*(pi1(2,k)*thetas(i)+pi2(2,k)*(1-thetas(i)))+pi1(3,k)*thetas(i)+pi2(3,k)*(1-thetas(i)));
                kappa2 = beta-kappa1;
                ecsum = ecsum+pit(i,j,k)*(log(A1(k)*kappa1+A2(k)*kappa2)/(1-beta)+termC);
            end
            C2(i,j) = thetas(i)*log(thetas(i))+(1-thetas(i))*log(1-thetas(i))+log(1-beta)+beta*ecsum;
        end
    end
    Cdiff = norm(C2-C)
    if Cdiff < 1e-08
        break
    end
end

figure(1)
plot(thetas,C2(:,1),'LineWidth',2)
hold on
plot(thetas,C2(:,2),'LineWidth',2)
plot(thetas,C2(:,3),'LineWidth',2)
plot(thetas,C2(:,4),'LineWidth',2)