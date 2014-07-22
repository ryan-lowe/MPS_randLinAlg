% ************************ one-site optimization **************************

function [A,E,Heff] = minimizeE_onesite(hsetj,Hleft,Hright,mpsJ)%add previous vector value to use as initial vector in eigs

DAl = size(Hleft,1); 
DAr = size(Hright,1); 
d = 2;

% calculation of Heff

Heff = 0; 
sizeHleft=size(Hleft)
sizehsetj=size(hsetj)
Heff = contracttensors(Hleft,3,2,hsetj,4,1); 
sizeHeff=size(Heff)
sizeHright=size(Hright)
Heff = contracttensors(Heff,5,3,Hright,3,2); 

sizeH=size(Heff)

%Heff is initially in the arrangement [a1,a2,o1,o2,b1,b2]

Heff = permute(Heff,[1,3,5,2,4,6]); %[1,3,2,4,6,5]
Heff = reshape(Heff,[DAl*DAr*d,DAl*DAr*d]); 

if (Heff==Heff')
    display('conj works');
else
    display('conj doesnt work');
end

Heff2=Heff;
%
% optimization
%
options.disp = 0; 

%initial vector used for eigs calculation
if size(mpsJ,1)>0
    %sizeJ=size(mpsJ)
    mpsJx=size(mpsJ,1);
    mpsJy=size(mpsJ,2);
    mpsJz=size(mpsJ,3);
    initV=reshape(mpsJ,[mpsJx*mpsJy*mpsJz,1]);
    options.v0=initV;
end

[A,E] = eigs(Heff,1,'sr',options); 

%A==A'

A=reshape(A,[DAl,DAr,d]);
%sizeA=size(A)