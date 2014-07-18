
n=6;
m=6;
k=4;
p=5;
A=[ 0 0 0 1 0 0; 
    0 0 0 0 0 1;
    0 0 0 0 1 0;
    0 0 0 0 2 1;
    0 0 0 0 0 0;
    0 0 1 0 0 0];
    
omega = randn(n,k+p);
Y = A*omega;
[Q, ~] = qr(Y,0);
% Q=zeros(m,k+p);
% %Use Gram-Schmidt
% Q(:,1)=Y(:,1);
% Q(:,1)=Q(:,1)/norm(Q(:,1));
% for i = 2:k+p
%     Q(:,i)=Y(:,i)-(dot(Q(:,i-1),Y(:,i)))/(dot(Y(:,i),Y(:,i)))*Q(:,i-1);
%     Q(:,i)=Q(:,i)/norm(Q(:,i));   
% end
