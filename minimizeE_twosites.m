% ************************ two-site optimization **************************

function [A1,A2,E,Heff] = minimizeE_twosites(hsetj,Hleft,Hright,dir)
%This file serves as the one-site optimization for minimizeE2.m
%INPUT: Hleft and Hright tensors, as calculated with updateHleft.m and
%updateHright.m, the current mpo site hsetj, and the current direction dir
%that the minimizeE2.m file is traversing
%OUTPUT: Heff, the contraction of Hleft, Hright, and hsetj, and the
%eigenvector A and eigenvalue E of Heff.

DAl = size(Hleft,1); 
DAr = size(Hright,1); 
d = 2;

% calculation of Heff
sizeH=size(hsetj)
% hsetj=contracttensors(hsetj1,3,3,hsetj2,3,2);
% hsetj=reshape(hsetj,[

Heff = contracttensors(Hleft,3,2,hsetj,6,1); 
Heff = contracttensors(Heff,7,3,Hright,3,2); 

%sizeH=size(Heff)

%Heff is initially in the arrangement [a1,a2,o1,o2,b1,b2]
sizeHeff=size(Heff)
Heff = permute(Heff,[1,7,3,4,2,8,5,6]); 
Heff = reshape(Heff,[DAl*DAr*d*d,DAl*DAr*d*d]); 


%
% optimization
%
options.disp = 0; 
[A1,E] = eigs(Heff,1,'sr',options); 

sizeA1=size(A1)
A1=reshape(A1,[DAl*d,DAr*d]);
[U,S,V]=svd(A1);


%dimensions after doing svd were not correct for rearrangement, so
%'inner dimension' constructed ( = D from original mpo)

%probably isn't correct, though (and still doesn't work for boundaries)
innerD=8;
innerDU=min(innerD,size(U,2));
innerDV=min(innerD,size(V,2));

U=U(:,1:innerDU);
S=S(1:innerDU,1:innerDV);
V=V(:,1:innerDV);

sizeU=size(U)
sizeS=size(S)
sizeV=size(V)
if strcmp(dir,'lr')
    A1=U;
    A2=S*V';   
else
    A1=V';
    A2=U*S;
end
% A1=U;
% A2=S*V';

A1=reshape(A1,[DAl,innerD,d]);
A2=reshape(A2,[innerD,DAr,d]);