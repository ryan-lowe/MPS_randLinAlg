function [E,mps]=minimizeE(mpo,D,precision)

N = size(mpo,2);
d = 2; 
mps = createrandommps(N,D,d); 
mps = prepare(mps,'rl');


% storage-initialization
Hstorage = initHstorage(mps,mpo,d);
 
% optimization sweeps 
count =0;
while 1
    Evalues = [];
clc
    % ****************** cycle 1: j -> j+1 (from 1 to N-1) **************** 
    for j = 1:(N-1)
        j=j
        % optimization
        Hleft = Hstorage{j};
        Hright = Hstorage{j+1};
        hsetj = mpo{j};
        [A,E] = minimizeE_onesite(hsetj,Hleft,Hright); 
        [A,U] = prepare_onesite(A,'lr');
        mps{j} = A; 
        Evalues = [Evalues,E];
 
        % storage-update 
        Hstorage{j} = updateCleft(Hleft,A,hsetj,A); 
     end
clc

% ****************** cycle 2: j -> j-1 (from N to 2) ****************** 
    for j = N:(-1):2
        j=j
        % minimization
        Hleft = Hstorage{j};
        Hright = Hstorage{j+1};
        hsetj = mpo{j}; 
        [A,E] = minimizeE_onesite(hsetj,Hleft,Hright); 
        [A,U] = prepare_onesite(A,'rl');
        mps{j} = A; 
        Evalues = [Evalues,E];

        % storage-update 
        Hstorage{j+1}=updateCright(Hright,A,hsetj,A); 
    end
    if (std(Evalues)/abs(mean(Evalues))<precision) 
        mps{1}=contracttensors(mps{1},3,2,U,2,1); 
        mps{1}=permute(mps{1},[1,3,2]);
        break;
    end
end