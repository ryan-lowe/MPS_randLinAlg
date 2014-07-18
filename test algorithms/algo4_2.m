n=100;
m=100;
% A=[ 0 0 0 1 0 0; 
%     0 0 0 0 0 1;
%     0 0 0 0 1 0;
%     0 0 0 0 2 1;
%     0 0 0 0 0 0;
%     0 0 1 0 0 0];

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

   
   
   
   

