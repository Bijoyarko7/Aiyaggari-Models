function [optk2,optc2,fx,f1,mpk,gterm,aggk,aggn] = aiyagari_egm_planner2(r,g,alpha,beta,delta,zgrid,piz,m,n,kgrid,optk2,optc2,pinv,f1,ng,kgrid2,cks,jks)

w = (1-alpha)*(r/alpha)^(alpha/(alpha-1));
chat = zeros(m,n); khat = zeros(m,n);

for iter = 1:10000
    optc = optc2;
    optk = optk2;
    emuc = zeros(m,n);
    for j=1:n
        for k=1:n
            for i=1:m
                emuc(i,j)=emuc(i,j)+piz(j,k)/optc(i,k);
            end
        end
    end
    for i=1:m
        for j=1:n
            chat(i,j)=1/(beta*(1+r-delta)*emuc(i,j)+beta*g);
            khat(i,j)=(chat(i,j)+kgrid(i)-w*zgrid(j))/(1+r-delta);
        end
    end
    for j=1:n
        ppc = pchip(khat(:,j),chat(:,j));
        for i=1:m
            if kgrid(i) > khat(1,j)
                optc2(i,j) = ppval(ppc,kgrid(i));
                optk2(i,j) = (1+r-delta)*kgrid(i)+w*zgrid(j)-optc2(i,j);
            else
                optk2(i,j) = 0.0;
                optc2(i,j) = (1+r-delta)*kgrid(i)+w*zgrid(j);
            end
        end
    end
    %[norm(optc2-optc) norm(optk2-optk)]
    if norm(optc2-optc)<1e-09
        break
    end
end

optkg=zeros(ng,n); optcg=zeros(ng,n);
for i=1:ng
    jk=jks(i);
    ck=cks(i);
    optkg(i,:)=(1-ck)*optk2(jk,:)+ck*optk2(jk+1,:);
    optcg(i,:)=(1-ck)*optc2(jk,:)+ck*optc2(jk+1,:);
end
[f1,aggk,aggn,gterm]=get_tmat_planner(kgrid2,zgrid,piz,optkg,optcg,f1,alpha);

mpk = alpha*(aggk/aggn)^(alpha-1);
fx(1) = mpk-r;
fx(2) = gterm-g;
end