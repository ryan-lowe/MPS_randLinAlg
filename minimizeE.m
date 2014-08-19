function [E,mps, vals, energy]=minimizeE(mpo,D,precision)
%This file finds the value of E that minimizes <y|mpo|y>/<y|y>, where y is
%any mps, and mpo is the (open boundary) mpo form of some Hamiltonian H. 
%This is done via a variational search in the MPS space, for the purpose 
%of finding the ground state of H.
%INPUT: an mpo (mpo form of a Hamiltonian H), matrix dimension D of the
%mpo, and the desired precision (code will stop once E value changes less
%than this amount in subsequent iterations)
%OUTPUT: minimum value E, and the corresponding minimum mps


N = size(mpo,2);
d = 2; 
mps = createrandommps(N,D,d); 

mps = prepare(mps,'rl');

% storage-initialization
Hstorage = initHstorage(mps,mpo,d);
vals=[];
energy=[];
count =0;
% optimization sweeps 
while 1
    Evalues = [];

    % ****************** cycle 1: j -> j+1 (from 1 to N-1) **************** 
    for j = 1:(N-1)
        
        % optimization
        Hleft = Hstorage{j};
        Hright = Hstorage{j+1};
        hsetj = mpo{j};
        
        %for first iteration, Hleft and Hright dimension modifications
        %by prepare.m cause initialization vector mps{j} for eigs
%         to be wrong dimension
%         if count==0
%             [A,E] = minimizeE_onesite(hsetj,Hleft,Hright,[]);             
%         else
            [A,E] = minimizeE_onesite(hsetj,Hleft,Hright,mps{j}); 
%         end
        [A,U] = prepare_onesite(A,'lr');
        
        mps{j} = A; 
        mps{j+1}=contracttensors(U,2,2,mps{j+1},3,1);
        Evalues = [Evalues,E];
 
        % storage-update 
        Hstorage{j+1} = updateCleft(Hleft,A,hsetj,A); 
        
    end
    count=count+1
   
% ****************** cycle 2: j -> j-1 (from N to 2) ****************** 
    for j = N:(-1):2
        
        % minimization
        Hleft = Hstorage{j};
        Hright = Hstorage{j+1};
        hsetj = mpo{j}; 
        [A,E] = minimizeE_onesite(hsetj,Hleft,Hright,mps{j}); 
        
        [A,U] = prepare_onesite(A,'rl');
        mps{j-1}=contracttensors(mps{j-1},3,2,U,2,1);
        mps{j-1}=permute(mps{j-1},[1,3,2]);
        
        mps{j} = A; 
        Evalues = [Evalues,E];

        % storage-update 
        Hstorage{j}=updateCright(Hright,A,hsetj,A); 
    end
    
    val = std(Evalues)/abs(mean(Evalues));
    energy=[energy real(E)];
    vals=[vals val];
    if (std(Evalues)/abs(mean(Evalues))<precision || count>=100) 
        mps{1}=contracttensors(mps{1},3,2,U,2,1); 
        mps{1}=permute(mps{1},[1,3,2]);
        break;
    end
end