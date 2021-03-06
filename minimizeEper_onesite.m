% ************************ one-site optimization **************************

function [A,E,Heff] = minimizeEper_onesite(hsetj,Hleft,Hright,Hleft2,Hright2,hsetId,mpsJ)
%This file serves as the one-site optimization for minimizeEper.m
%INPUT: Hleft and Hright tensors, as calculated with updateHleft.m and
%updateHright.m, the current mpo site hsetj, and the previous mps mpsJ,
%used as an 'iniital vector' for the eigs calculation. Hleft2 is the Hleft
%for the mps-mps contraction (which is non-zero because of periodic
%boundary conditions)
%OUTPUT: Heff, the contraction of Hleft, Hright, and hsetj, and the
%eigenvector A and eigenvalue E of Heff.

DAl = size(Hleft,1); 
DAr = size(Hright,1); 
d = 2;

% calculation of Heff

Heff = contracttensors(Hleft,3,2,hsetj,4,1); 
Heff = contracttensors(Heff,5,3,Hright,3,2); 

Heff2=contracttensors(Hleft2,3,2,hsetId,4,1);
Heff2=contracttensors(Heff2,5,3,Hright2,3,2);

Heff2 = permute(Heff2,[1,5,4,2,6,3]);
Heff2 = reshape(Heff2,[DAl*DAr*d,DAl*DAr*d]); 

%Heff is initially in the arrangement [a1,a2,o2,o1,b1,b2]

Heff = permute(Heff,[1,5,4,2,6,3]); 
Heff = reshape(Heff,[DAl*DAr*d,DAl*DAr*d]); 


if (norm(Heff-Heff')<1e-10)
    display('conj works');
else
    display('conj doesnt work');
end

% optimization
options.disp = 0; 

%initial vector used for eigs calculation (as long as mpsJ isn't empty)
if ~isempty(mpsJ)
    [mpsJx mpsJy mpsJz]=size(mpsJ);
    initV=reshape(mpsJ,[mpsJx*mpsJy*mpsJz,1]);
    options.v0=initV;
end

[A,E] = eigs(Heff,Heff2,1,'sr',options);


A=reshape(A,[DAl,DAr,d]);
