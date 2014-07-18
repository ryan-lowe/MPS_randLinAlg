
%testing increasing size of A
k=2;
%q=1;
%m=100;

data1=zeros(5,1);
data2=zeros(5,1);

for z=1:5
     m=z*100;
     n=m;
    %k=2*z;
    
    A=rand(m);

    epsilon=1e-4;
    %works for r<=3
    r=10;

    omega=randn(n,r);
    y=zeros(m,r);
    for i=1:r
        y(:,i)=A*omega(:,i);
    end

    j=0;
    q=zeros(m,r);

    %find max of ||y(j+1)|| to ||y(j+r)||
    max=norm(y(:,1));
    for k=j+1:j+r
       if norm(y(:,k))>max
           max=norm(y(:,k));
       end
    end
    while max>epsilon/(10*sqrt(2/pi))
       j=j+1;
       if j==1
           q(:,j)=y(:,j)/norm(y(:,j));
           %3-d (ie. a sequence of 2-d matrices) Q matrix not used, as 
           %matrices in a sequence must be of same dimension
           Q=q(:,j);
           omega2=randn(n,1);
           y(:,j+r)=A*omega2;
       else
           y(:,j)=y(:,j)-Q*Q'*y(:,j);
           q(:,j)=y(:,j)/norm(y(:,j));
           Q=[Q q(:,j)];
           omega3=randn(n,1);
           y(:,j+r)=A*omega3-Q*Q'*A*omega3;
       end
       for i=j+1:j+r-1
           y(:,i)=y(:,i)-q(:,j)*dot(q(:,j),y(:,i));
       end

       %find max of ||y(j+1)|| to ||y(j+r)||
       max=norm(y(:,j+1));
       for k=j+1:j+r
           if norm(y(:,k))>max
                max=norm(y(:,k));
           end
       end
    end
    
%     A=rand(m);
% 
%     omega=randn(n,l);
%     Y=(A*A')^q*A*omega;
%     [Q, ~]=qr(Y,0);
    
    d1=0;
    for i=1:3
    tic
    B=Q'*A;
    [U1, sigma, V1]=svd(B);
    V=V1';
    U=Q*U1;
    d1=d1+toc;
    end;
    data1(z)=d1/3;
    
    d2=0;
    for j=1:3
    tic
    [U2, sigma2, V2]=svd(A);
    d2=d2+toc;
    end
    data2(z)=d2/3;

end
