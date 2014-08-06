function [E,mps,vals,energy]=minimizeEper(mpo,D,precision)
%This file finds the value of E that minimizes <y|mpo|y>/<y|y>, where y is
%any mps, and mpo is the (open boundary) mpo form of some Hamiltonian H. 
%This is done via a variational search in the MPS space, for the purpose 
%of finding the ground state of H.

%In this code, the input MPO has periodic boundary conditions.

%INPUT: an mpo (mpo form of a Hamiltonian H), matrix dimension D of the
%mpo, and the desired precision (code will stop once E value changes less
%than this amount in subsequent iterations)
%OUTPUT: minimum value E, and the corresponding minimum mps

N = size(mpo,2);
d = 2; 
mps = createrandommps(N,D,d); 
mps = prepare(mps,'rl');

mpoId=cell(1,N);
for i=1:N
    mpoId{i}=eye(d,d);
    mpoId{i}=reshape(mpoId{i},[1,1,d,d]);
end

% storage-initialization
Hstorage = initHstorage(mps,mpo,d);
Hstorage2 = initHstorage(mps,mpoId,d);
 sizeHS2=size(Hstorage2)
 Hstorage2=Hstorage2
% optimization sweeps 
count = 0;
vals=[];
energy=[];
while 1
    Evalues = [];

    % ****************** cycle 1: j -> j+1 (from 1 to N-1) **************** 
    for j = 1:(N-1)
        %count=count+1
        %j=j
        % optimization
        Hleft = Hstorage{j};
        Hright = Hstorage{j+1};
        hsetj = mpo{j};
        
        Hleft2 = Hstorage2{j};
        Hright2 = Hstorage2{j+1};
        hsetId = mpoId{j};
        
        %for first iteration, Hleft and Hright dimension modifications
        %by prepare.m cause initialization vector mps{j} for eigs
        %to be wrong dimension
        if count==0
            [A,E] = minimizeEper_onesite(hsetj,Hleft,Hright,Hleft2,Hright2,hsetId,[]);             
        else
            [A,E] = minimizeEper_onesite(hsetj,Hleft,Hright,Hleft2,Hright2,hsetId,mps{j}); 
        end
        [A,U] = prepare_onesite(A,'rl');
        mps{j} = A; 
        Evalues = [Evalues,E];
 
        % storage-update 
        Hstorage{j+1} = updateCleft(Hleft,A,hsetj,A); 
    end
    count=count+1
    val=std(Evalues)/abs(mean(Evalues));
    vals=[vals val];
    energy=[energy real(E)];
    if (std(Evalues)/abs(mean(Evalues))<precision || count>=100) 
        mps{1}=contracttensors(mps{1},3,2,U,2,1); 
        mps{1}=permute(mps{1},[1,3,2]);
        break;
    end
    

end