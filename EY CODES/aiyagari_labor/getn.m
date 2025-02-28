function z = getn(n,opthg,f1,jK,cK,jA,Ngrid,zgrid)

jN = binarySearch(Ngrid,n);
cN = (n-Ngrid(jN))/(Ngrid(jN+1)-Ngrid(jN));
[nk,nz] = size(f1);

opthgg = (1-cK)*(1-cN)*opthg(:,:,jK,jA,jN)+cK*(1-cN)*opthg(:,:,jK+1,jA,jN)+...
        (1-cK)*cN*opthg(:,:,jK,jA,jN+1)+cK*cN*opthg(:,:,jK+1,jA,jN+1);
aggn = 0;
for i=1:nk
    for j=1:nz
        aggn = aggn+zgrid(j)*opthgg(i,j)*f1(i,j);
    end
end
z = n-aggn;

end