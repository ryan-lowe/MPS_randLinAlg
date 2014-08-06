function [E,mps,vals, energy]=minimizeE2(mpo,D,precision)
%This file finds the value of E that minimizes <y|mpo|y>/<y|y>, where y is
%any mps, and mpo is the (open boundary) mpo form of some Hamiltonian H. 
%This is done via a variational search in the MPS space, for the purpose 
%of finding the ground state of H.

%This code is slightly altered to compute 2 sites at a time, for increased
%stability benefits

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
energy=[];
vals=[];
count=0;
% optimization sweeps 
while 1
    Evalues = [];

    % ****************** cycle 1: j -> j+1 (from 1 to N-1) **************** 
    for j = 1:(N-1)
        
        % optimization
        Hleft = Hstorage{j};
        Hright = Hstorage{j+2};
        hsetj1 = mpo{j};
        hsetj2 = mpo{j+1};
        
        hsetj=contracttensors(hsetj1,4,2,hsetj2,4,1);
        %sizeH1=size(hsetj)
        hsetj=permute(hsetj,[1,4,2,3,5,6]);
        HS1=size(hsetj,1);
        HS2=size(hsetj,2);
        %hsetj=reshape(hsetj,[HS1,HS2,4,4]); %probably not right
        
        %sizem1=size(mps{1})
        %sizem2=size(mps{2})
        
        [A1,A2,E] = minimizeE_twosites(hsetj,Hleft,Hright,'lr'); 
        [A1,U] = prepare_onesite(A1,'lr');
        mps{j} = A1; 
        %[A2,U] = prepare_onesite(A2,'lr');
        mps{j+1} = A2;
        Evalues = [Evalues,E];
 
        % storage-update 
        Hstorage{j+1} = updateCleft(Hleft,A1,hsetj1,A1); 
        
    end
   count=count+1

    % ****************** cycle 2: j -> j-1 (from N to 2) ****************** 
    for j = N:(-1):2

        % minimization
        Hleft = Hstorage{j-1};
        Hright = Hstorage{j+1};
        hsetj1 = mpo{j};
        hsetj2 = mpo{j-1};
%         h1=size(hsetj1)
%         h2=size(hsetj2)
        hsetj=contracttensors(hsetj2,4,2,hsetj1,4,1);
        hsetj=permute(hsetj,[1,4,2,3,5,6]);
        HS1=size(hsetj,1);
        HS2=size(hsetj,2);
        hsetj=reshape(hsetj,[HS1,HS2,4,4]); %probably not right
        
        [A1,A2,E] = minimizeE_twosites(hsetj,Hleft,Hright,'rl'); 
        [A2,U] = prepare_onesite(A2,'rl');
        mps{j} = A2; 
        %[A2,U] = prepare_onesite(A2,'rl');
        mps{j-1} = A1;
        Evalues = [Evalues,E];

        % storage-update 
        Hstorage{j}=updateCright(Hright,A2,hsetj1,A2); 

    end
    
    val=std(Evalues)/abs(mean(Evalues));
    vals=[vals val];
    energy=[energy real(E)];    
    if (std(Evalues)/abs(mean(Evalues))<precision || count>=100) 
        mps{1}=contracttensors(mps{1},3,2,U,2,1); 
        mps{1}=permute(mps{1},[1,3,2]);
        break;
    end
end