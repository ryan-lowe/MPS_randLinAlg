% ************************ one-site optimization **************************

function [A,E,Heff] = minimizeE_onesite(hsetj,Hleft,Hright)%add previous vector value to use as initial vector in eigs

DAl = size(Hleft,1); 
DAr = size(Hright,1); 
d = 2;

% calculation of Heff

Heff = 0; 
%sizeHleft=size(Hleft)
%sizehsetj=size(hsetj)
Heff = contracttensors(Hleft,3,2,hsetj,4,1); 
%sizeHeff=size(Heff)
%sizeHright=size(Hright)
Heff = contracttensors(Heff,5,3,Hright,3,2); 
Heff = permute(Heff,[1,3,5,2,4,6]); 
Heff = reshape(Heff,[DAl*DAr*d,DAl*DAr*d]); 

%
% optimization
%
options.disp = 0; 
[A,E] = eigs(Heff,1,'sr',options); 

A=reshape(A,[DAl,DAr,d]);
%sizeA=size(A)