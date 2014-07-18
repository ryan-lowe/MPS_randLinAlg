function [mps] = compressSVD(mps,k)

N=length(mps);

[mps,~]=prepare(mps,'lr');

for i=1:N-1
    [D1,D2,d]=size(mps{i});
    mpstemp=permute(mps{i},[3,1,2]);
    mpstemp=reshape(mpstemp,[d*D1,D2]);
    [U,S,V]=svd2(mpstemp);
    if (length(S)>k)
        S=S(1:k,1:k);
        S=S/norm(S);
        U=U(:,1:k);
        V=V(1:k,:);
        mps{i}=U;
        mps{i}=reshape(mps{i},[d,D1,k]);
        mps{i+1}=contracttensors(S*V,2,2,mps{i+1},3,1);
    else
        mps{i}=reshape(mps{i},[d,D1,D2]);
    end
    mps{i}=permute(mps{i},[2,3,1]);
end
[D1,D2,d]=size(mps{N});
mpstemp=permute(mps{N},[3,1,2]);
mpstemp=reshape(mpstemp,[d*D1,D2]);
[U,~,~]=svd2(mpstemp);
if (length(S)>k)
    U=U(:,1:k);
    mps{N}=U;
    mps{N}=reshape(mps{N},[d,D1,k]);
else
    mps{N}=reshape(mps{N},[d,D1,D2]);
end
mps{N}=permute(mps{N},[2,3,1]);
