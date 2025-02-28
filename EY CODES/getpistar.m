function y = getpistar(x,thetas,A1,A2,beta1,beta2,alpha)

thetap = (beta1*thetas)./(beta1*thetas+beta2*(1-thetas));
pip = pchip(thetas,x,thetap);
bhat = beta1*thetas+beta2*(1-thetas);
prod = A1^(1/(alpha-1))/(A1^(1/(alpha-1))+A2^(1/(alpha-1)));
kappa1 = 1-(1+(A1/A2)^(1/(alpha-1)))*x;
kappa2 = 1-(1+(A1/A2)^(1/(alpha-1)))*pip;
y = bhat*prod.*kappa1-kappa2.*x;

sum(y.*y)