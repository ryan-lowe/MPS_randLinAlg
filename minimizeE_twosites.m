% ************************ two-site optimization **************************

function [A1,A2,E,Heff] = minimizeE_twosites(hsetj,Hleft,Hright,dir)
DAl = size(Hleft,1); 
DAr = size(Hright,1); 
d = 2;

% calculation of Heff

% hsetj=contracttensors(hsetj1,3,3,hsetj2,3,2);
% hsetj=reshape(hsetj,[

Heff = 0; 
Heff = contracttensors(Hleft,3,2,Hright,3,2); 
Heff = contracttensors(Heff,5,5,hsetj,3,3); 
Heff = permute(Heff,[1,3,5,2,4,6]); 
Heff = reshape(Heff,[DAl*DAr*d,DAl*DAr*d]); 

%
% optimization
%
% options.disp = 0; 
[A1,E] = eigs(Heff,1,'sr',options); 

[U,S,V]=svd(Heff);
if strcmp(dir,'lr')
    A2=S*V';
else
    A2=U*S;
end
% A1=U;
% A2=S*V';

A1=reshape(A1,[DAl,DAr,d]);
A2=reshape(A2,[DAl,DAr,d]);