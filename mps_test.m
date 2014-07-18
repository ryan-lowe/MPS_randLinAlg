
L=10;
B=zeros(2,1,2,2);
A=zeros(L-2,2,2,2);
C=rand(2,2,2,2,2,2,2,2,2,2);
a=2;

C=reshape(C,[2,2^(L-1)]);
[U S V]=svd(C);
%B(1,:,:,:)=reshape(U,[1,2,2]);
C=S*V';

C=reshape(C,[a*2,2^(L-2)]);

[U S V]=svd(C,0);

%U=reshape(U,[a,2,4]);
%A(1,:,:,:)=U;
%C=S*V';
%C=reshape(C,[4*2,2^(L-3)]);

% [U S V]=svd(C);
% 
% U=reshape(U,[4,2,8]);

% for i=2:L-1
%     C=reshape(C,[a*2,2^(L-i)]);
%     [U S V]=svd(C);
%     A(i-1,:,:,:)=reshape(U,[a,2,a]);
%     C=S*V';
% end
% 
% C=reshape