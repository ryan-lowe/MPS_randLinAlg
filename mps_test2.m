% THIS IS A TEST SCRIPT
function [mps,V]=mps_test2(coeff,d,N)
% % d=2;N=10;
% % coeff=rand(d^N,1);
% % coeff=coeff/norm(coeff);

mps=cell(N,1);
V=coeff;r=1;
for ind=1:N
    ind
    V=reshape(V,[r*d,d^(N-ind)]);
    [U,S,V]=svd2(V);
    size(U)
    size(S)
    size(V)
    rold=r;r=size(S,1);
    V=S*V;
    U=reshape(U,[rold,d,r]);
    mps{ind}=permute(U,[1,3,2]);
end
% 
% U=coeff;r=1;
% for ind=1:N
%     U=reshape(U,[d^(N-ind),r*d]);
%     [U,S,V]=svd2(U);
%     rold=r;r=size(S,1);
%     U=U*S;
%     V=reshape(V,[r,d,rold]);
%     mps{ind}=permute(V,[1,3,2]);
% end