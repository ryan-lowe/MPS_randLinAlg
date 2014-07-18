function [E,mps]=minimizeE2(mpo,D,precision,mpsB)

N = size(mpo,2); 
d = 2; 
mps = createrandommps(N,D,d); 
mps = prepare(mps,'rl');

% storage-initialization
Hstorage = initHstorage(mps,mpo,d);
if ~isempty(mpsB), Cstorage = initCstorage(mps,[],mpsB,N); end

% optimization sweeps 
while 1
    Evalues = [];

    % ****************** cycle 1: j -> j+1 (from 1 to N-1) **************** 
    for j = 1:(N-1)
        % projector-calculation 
        if ~isempty(mpsB)
            B = mpsB{j};
            Cleft = Cstorage{j};
            Cright = Cstorage{j+1}; 
        end

        % optimization
        Hleft = Hstorage{j};
        Hright = Hstorage{j+2};
        hsetj1 = mpo{j};
        hsetj2 = mpo{j+1};
        hsetj=contracttensors(hsetj1,3,3,hsetj2,3,2);
        DAl = size(Hleft,1); 
        DAr = size(Hright,1);
        HS1 = size(hsetj1,1);
        hsetj=contracttensors(hsetj1,3,3,hsetj2,3,2);
        hsetj=reshape(hsetj,[DAl,HS1*d,DAr]); %probably not right
        [A1,A2,E] = minimizeE_twosites(hsetj,Hleft,Hright,'lr'); 
        [A1,U] = prepare_onesite(A1,'lr');
        mps{j} = A1; 
        %[A2,U] = prepare_onesite(A2,'lr');
        mps{j+1} = A2;
        Evalues = [Evalues,E];
 
        % storage-update 
        Hstorage{j+1} = updateCleft(Hleft{m},A1,hsetj1,A1); 
        if ~isempty(mpsB) 
            Cstorage{j+1} = updateCleft(Cleft,A1,[],B);
        end
    end
   
    % ****************** cycle 2: j -> j-1 (from N to 2) ****************** 
    for j = N:(-1):2
        % projector-calculation 
        if ~isempty(mpsB)
            B = mpsB{j};
            Cleft = Cstorage{j};
            Cright = Cstorage{j+1}; 
        end
        
        % minimization
        Hleft = Hstorage{j-1};
        Hright = Hstorage{j+1};
        hsetj1 = mpo{j};
        hsetj2 = mpo{j-1};
        DAl = size(Hleft,1); 
        DAr = size(Hright,1);
        HS1 = size(hsetj1,1);
        hsetj=contracttensors(hsetj1,3,3,hsetj2,3,2);
        hsetj=reshape(hsetj,[DAl,HS1*d,DAr]); %probably not right
        [A1,A2,E] = minimizeE_twosites(hsetj,Hleft,Hright,'rl'); 
        [A1,U] = prepare_onesite(A1,'rl');
        mps{j} = A1; 
        %[A2,U] = prepare_onesite(A2,'rl');
        mps{j-1} = A2;
        Evalues = [Evalues,E];

        % storage-update 
        Hstorage{j}=updateCright(Hright{m},A1,hsetj1,A1); 
        if ~isempty(mpsB) 
            Cstorage{j}=updateCright(Cright,A,[],B);
        end
    end
    if (std(Evalues)/abs(mean(Evalues))<precision) 
        mps{1}=contracttensors(mps{1},3,2,U,2,1); 
        mps{1}=permute(mps{1},[1,3,2]);
        break;
    end
end