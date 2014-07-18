n=6;
m=6;
k=4;
p=5;
A=[ 0 0 0 1 0 0; 
    0 0 0 0 0 0;
    0 0 2 0 1 0;
    1 0 0 0 0 1;
    0 0 1 0 0 0;
    0 0 0 1 0 0];
q=1;

omega=randn(n,l);
Y=(A*A')^q*A*omega;
[Q, R]=qr(Y);

B=Q'*A*Q;
[V, lambda]=eig(B);
U=Q*V;