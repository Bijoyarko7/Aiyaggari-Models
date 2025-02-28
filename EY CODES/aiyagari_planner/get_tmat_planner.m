function [fout,kbar,nbar,gterm] = get_tmat_planner(kgrid,zgrid,piz,optk,optc,f1,alpha)

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
for t=1:100000
    f0 = f1;
    f1 = zeros(nk,nz);
    for i=1:nk
        for j=1:nz
            jlo=jlos(i,j);
            wgt=wgts(i,j);
            for k=1:nz
                f1(jlo,k)=f1(jlo,k)+piz(j,k)*wgt*f0(i,j);
                f1(jlo+1,k)=f1(jlo+1,k)+piz(j,k)*(1-wgt)*f0(i,j);
            end
        end
    end
    if (norm(f1-f0))<1e-08
        break
    end
end
fout = f1;
kbar=0;
nbar=0;
for i=1:nk
    kbar=kbar+kgrid(i)*sum(f1(i,:));
    for j=1:nz
        nbar=nbar+zgrid(j)*f1(i,j);
    end
end
gterm=0;
for i=1:nk
    for j=1:nz
        kp=optk(i,j);
        jp=jlos(i,j);
        wp=1-wgts(i,j);
        cp1=(1-wp)*optc(jp,1)+wp*optc(jp+1,1);
        cp2=(1-wp)*optc(jp,2)+wp*optc(jp+1,2);
        gterm=gterm+piz(j,1)/cp1*(alpha*(alpha-1)*kbar^(alpha-2)*kp+alpha*(1-alpha)*kbar^(alpha-1)*zgrid(1))*f1(i,j)+...
            piz(j,2)/cp2*(alpha*(alpha-1)*kbar^(alpha-2)*kp+alpha*(1-alpha)*kbar^(alpha-1)*zgrid(2))*f1(i,j);
    end
end

end
