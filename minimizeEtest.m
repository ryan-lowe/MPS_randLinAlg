function [E,mps1, vals, energy]=minimizeEtest(mpo,D,precision)
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
mps1 = createrandommps(N,D,d); 
mps1 = prepare(mps1,'rl');
mps2=mps1;

% storage-initialization
Hstorage1 = initHstorage(mps1,mpo,d);
Hstorage2 = Hstorage1;
vals=[];
energy=[];
count =0;
% optimization sweeps 
while 1
    Evalues = [];

    % ****************** cycle 1: j -> j+1 (from 1 to N-1) **************** 
    for j = 1:(N-1)
        j
        % optimization
        Hleft1 = Hstorage1{j};
        Hright1 = Hstorage1{j+1};
        hsetj1 = mpo{j};
        
        hsetj2 = mpo{j+1};
        Hleft2=Hstorage2{j};
        Hright2=Hstorage2{j+2};
        

        %%%%%%%%%%%% 1 site
        if count==0
            [A,E] = minimizeE_onesite(hsetj1,Hleft1,Hright1,[]);             
        else
            [A,E] = minimizeE_onesite(hsetj1,Hleft1,Hright1,mps1{j}); 
        end
        [A,U] = prepare_onesite(A,'lr');
        mps1{j} = A; 
        Evalues = [Evalues,E];
 
        Hstorage1{j+1} = updateCleft(Hleft1,A,hsetj1,A); 
        
        %%%%%%%%%%%% 2 site
        hsetj=contracttensors(hsetj1,4,2,hsetj2,4,1);
        hsetj=permute(hsetj,[1,4,2,3,5,6]);
    
        if count==0
          [A1,A2,E] = minimizeE_twosites(hsetj,Hleft2,Hright2,'lr',[]); 
        else
          mpsJ=contracttensors(mps2{j},3,2,mps2{j+1},3,1);
          [A1,A2,E] = minimizeE_twosites(hsetj,Hleft2,Hright2,'lr',mpsJ); 
        end        
        [A1,U] = prepare_onesite(A1,'lr');
        mps2{j} = A1; 
        mps2{j+1} = A2;
        Evalues = [Evalues,E];
 
        Hstorage2{j+1} = updateCleft(Hleft2,A1,hsetj1,A1);
        
        sizeA=size(A)
        sizeA1=size(A1)
       % Adiff=A-A1
    end
    count=count+1

% ****************** cycle 2: j -> j-1 (from N to 2) ****************** 
    for j = N:(-1):2
        j
        % optimization
        Hleft1 = Hstorage1{j};
        Hright1 = Hstorage1{j+1};
        hsetj1 = mpo{j}; 
        [A,E] = minimizeE_onesite(hsetj1,Hleft1,Hright1,mps1{j}); 
        [A,U] = prepare_onesite(A,'rl');
        mps1{j} = A; 
        Evalues = [Evalues,E];

        % storage-update 
        Hstorage1{j}=updateCright(Hright1,A,hsetj1,A);     
        
        %%%%%%%%%%%% 2 site
        Hleft2 = Hstorage2{j-1};
        Hright2 = Hstorage2{j+1};
        hsetj2 = mpo{j-1};

        hsetj=contracttensors(hsetj2,4,2,hsetj1,4,1);
        hsetj=permute(hsetj,[1,4,2,3,5,6]);
        mpsJ=contracttensors(mps2{j-1},3,2,mps2{j},3,1);
        
        [A1,A2,E] = minimizeE_twosites(hsetj,Hleft2,Hright2,'rl',mpsJ); 
        [A2,U] = prepare_onesite(A2,'rl');
        mps2{j} = A2; 
        mps2{j-1} = A1;
        Evalues = [Evalues,E];

        Hstorage2{j}=updateCright(Hright2,A2,hsetj1,A2); 
        
        sizeA=size(A)
        sizeA2=size(A2)
        %Adiff=A-A2
    end
    val = std(Evalues)/abs(mean(Evalues));
    energy=[energy real(E)];
    vals=[vals val];
    %Evalues = Evalues
    if (std(Evalues)/abs(mean(Evalues))<precision || count>=5) 
%         mps{1}=contracttensors(mps{1},3,2,U,2,1); 
%         mps{1}=permute(mps{1},[1,3,2]);
        break;
    end
end