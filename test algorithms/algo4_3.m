n=6;
m=6;
l=4;
q=1;
A=[ 0 0 0 1 0 0; 
    0 0 0 0 0 1;
    0 0 0 0 1 0;
    0 0 0 0 2 1;
    0 0 0 0 0 0;
    0 0 1 0 0 0];

omega=randn(n,l);
Y=(A*A')^q*A*omega;
[Q, ~]=qr(Y,0);


