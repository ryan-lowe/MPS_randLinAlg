% ************************ two-site optimization **************************

function [A1,A2,E,Heff] = minimizeE_twosites(hsetj,Hleft,Hright,dir,mpsJ)
%This file serves as the one-site optimization for minimizeE2.m
%INPUT: Hleft and Hright tensors, as calculated with updateHleft.m and
%updateHright.m, the current mpo site hsetj, and the current direction dir
%that the minimizeE2.m file is traversing
%OUTPUT: Heff, the contraction of Hleft, Hright, and hsetj, and the
%eigenvector A and eigenvalue E of Heff.

DAl = size(Hleft,1); 
DAr = size(Hright,1); 
d = 2;
Hl=size(hsetj,1);
Hr=size(hsetj,2);

% calculation of Heff

hsetj=reshape(hsetj,[Hl,Hr,d,d,d,d]);
Heff = contracttensors(Hleft,3,2,hsetj,6,1); 
Heff = contracttensors(Heff,7,3,Hright,3,2); 

Heff=reshape(Heff,[DAl,DAl,d,d,d,d,DAr,DAr]);

%Heff is initially in the arrangement [a1,a2,o1,o2,b1,b2]

Heff = permute(Heff,[1,3,5,7,2,4,6,8]); 
Heff = reshape(Heff,[DAl*DAr*d*d,DAl*DAr*d*d]); 

% if (norm(Heff-Heff')<1e-10)
%     display('conj works');
% else
%     display('conj doesnt work');
% end

%initial vector used for eigs calculation (as long as mpsJ isn't empty)
if ~isempty(mpsJ)
    [mpsJx mpsJy mpsJz mpsJw mpsJq]=size(mpsJ);
    initV=reshape(mpsJ,[mpsJx*mpsJy*mpsJz*mpsJw*mpsJq,1]);
    options.v0=initV;
end

options.disp = 0; 
[A1,E] = eigs(Heff,1,'sr',options); 

A1=reshape(A1,[DAl*d,DAr*d]);
[U,S,V]=svd(A1);


if strcmp(dir,'lr')
    A1=U;
    A2=S*V';   
else
    A2=V';
    A1=U*S;
end
innerD1=(size(A1,1)*size(A1,2))/(d*DAl);
innerD2=(size(A2,1)*size(A2,2))/(d*DAr);


A1=reshape(A1,[DAl,innerD1,d]);
A2=reshape(A2,[innerD2,DAr,d]);
% sizeA1=size(A1)
% sizeA2=size(A2)