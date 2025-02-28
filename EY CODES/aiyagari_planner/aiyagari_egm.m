function y = aiyagari_egm(r,alpha,beta,delta,zgrid,piz,m,n,kgrid,pinv)

w = (1-alpha)*(r/alpha)^(alpha/(alpha-1));
optc2 = zeros(m,n); chat = zeros(m,n); khat = zeros(m,n); optk2 = zeros(m,n);

for i=1:n
    optc2(:,i) = (r-delta)*kgrid+w*zgrid(i);
    optk2(:,i) = kgrid(i);
end
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
            chat(i,j)=1/(beta*(1+r-delta)*emuc(i,j));
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

kgrid2=linspace(kgrid(1),1500,50000); optkg=zeros(50000,n);
for i=1:50000
    jk=binarySearch(kgrid,kgrid2(i));
    jk=max(1,min(jk,length(kgrid)-1));
    ck=(kgrid2(i)-kgrid(jk))/(kgrid(jk+1)-kgrid(jk));
    optkg(i,:)=(1-ck)*optk2(jk,:)+ck*optk2(jk+1,:);
end
f1=zeros(50000,n);
for i=1:1000
    f1(i,1:n)=pinv/1000;
end
[f1,aggk,aggn]=get_tmat(kgrid2,zgrid,piz,optkg,f1);

mpk = alpha*(aggk/aggn)^(alpha-1);
y = mpk-r;
[r mpk aggk aggn y]
end