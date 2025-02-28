function [fout,kbar] = get_tmat_dynamic(kgrid,zgrid,piz,optk,f0)

nk = length(kgrid); [nz,~] = size(piz);
jlos=zeros(nk,nz); wgts=zeros(nk,nz);
for j=1:nz
    for i=1:nk
        jlos(i,j)=binarySearch(kgrid,optk(i,j));
        if jlos(i,j)<1
            jlos(i,j)=1;
            wgts(i,j)=1.0;
        elseif jlos(i,j)>=nk
            jlos(i,j)=nk-1;
            wgts(i,j)=0.0;
        else
            wgts(i,j)=1-(optk(i,j)-kgrid(jlos(i,j)))/(kgrid(jlos(i,j)+1)-kgrid(jlos(i,j)));
        end
    end
end
f1 = zeros(nk,nz);
for i=1:nk
    for j=1:nz
        jlo=jlos(i,j);
        wgt=wgts(i,j);
        for k=1:nz
            f1(jlo,k) = f1(jlo,k)+piz(j,k)*wgt*f0(i,j);
            f1(jlo+1,k) = f1(jlo+1,k)+piz(j,k)*(1-wgt)*f0(i,j);
        end
    end
end
fout = f1;
kbar = 0;
for i=1:nk    
    for j=1:nz
        kbar = kbar+optk(i,j)*f0(i,j);
    end
end

end
