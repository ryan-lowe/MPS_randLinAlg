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
q=1;

omega=randn(n,2*k);
Y=(A*A')^q*A*omega;
[Q R]=qr(Y);

B=Q'*A;
[U1, sigma, V1]=svd(B);
V=V1';
U=Q*U1;
