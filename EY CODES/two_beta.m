%solving dynamics for two-beta economy

thetas = linspace(0.001,0.999,51);
A1 = 1;
A2 = 1;
beta1 = 0.95;
beta2 = 0.98;
alpha = 1/3;

options = optimset('fsolve');
options = optimset('Display','off','TolF',1e-08);
pistar = fsolve(@(x) getpistar(x,thetas,A1,A2,beta1,beta2,alpha),alpha*(thetas*beta1+(1-thetas)*beta2),...
    options);
c1star = thetas.*(1-(1+(A1/A2)^(1/(alpha-1)))*pistar);
c2star = (1-thetas)./thetas.*c1star;
%pistar = fsolve(@(x) getpistar(x,thetas,A1,A2,beta1,beta2,alpha),pistar,options);
figure(1)
plot(thetas,pistar,'LineWidth',2)
title('\pi(\theta)','FontSize',13)
xlabel('\theta','FontSize',13)
figure(2)
plot(thetas,c1star,'LineWidth',2)
hold
plot(thetas,c2star,'LineWidth',2)
title('c_i/(A_1k_1^{\alpha} + A_2k_2^{\alpha})','FontSize',13)
xlabel('\theta','FontSize',13)
legend('Low \beta','High \beta','FontSize',13)
figure(3)
plot(thetas,beta1*thetas./(beta1*thetas+beta2*(1-thetas))-thetas,'LineWidth',2)
title('\theta^{\prime}','FontSize',13)
xlabel('\theta','FontSize',13)
hold
plot(thetas,0*thetas,'k--','LineWidth',2)
