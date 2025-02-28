function [aggk,aggn,x,f1] = aiyagari_projection_2(r,alpha,beta,delta,zgrid,piz,m,kgrid,x0,f1)

w = (1-alpha)*(r/alpha)^(alpha/(alpha-1));
aggn = 0.5*zgrid(1)+0.5*zgrid(2); 
options = optimoptions('fsolve','Display','none');

%x0(1:m) = kgrid;
%x0(m+1:2*m) = kgrid;
%x0(2*m+1:3*m) = -sqrt(kgrid);
%x0(3*m+1:4*m) = -sqrt(kgrid);
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
[f1,aggk,aggn]=get_tmat(kgrid2,zgrid,piz,optkg,f1);

% kp = 11;
% zt = 1;
% rng(538);
% T = 1000000;
% e = rand(T,1);
% pp1 = pchip(kgrid,optk2(:,1));
% pp2 = pchip(kgrid,optk2(:,2));
% kts = zeros(T,1);
% aggk = 0;
% for kkk = 1:T
%     kts(kkk) = kp;
%     if zt == 1
%         kp = ppval(pp1,kts(kkk));
%         if (e(kkk) < piz(1,1))
%             zt = 1;
%         else
%             zt = 2;
%         end
%     else
%         kp = ppval(pp2,kts(kkk));
%         if (e(kkk) < piz(2,1))
%             zt = 1;
%         else
%             zt = 2;
%         end
%     end
%     aggk = aggk+kts(kkk)/T;
% end

y = alpha*aggk^(alpha-1)*aggn^(1-alpha)-r;
[r aggk y]