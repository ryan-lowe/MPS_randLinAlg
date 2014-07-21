function [mpo]=mpo_ising(N,J,h)

% CREATE THE MPO OPERATOR FOR THE ISING HAMILTONIAN WITH N QUBITS AND
% GIVEN PARAMETERS;
% H=\SUM_sites J*SX_{i}*SX_{i+1}-H*SZ_{i}

mpo=cell(1,N);
d=2;
sx=[0,1;1,0]; sy=[0,-1i;1i,0]; sz=[1,0;0,-1]; id=eye(2); 

%this part correspond to left side of Eq.183 in Schollwock review
aux=zeros(1,3,d,d);
for d_index2=1:d
    for d_index1=1:d
        aux(1,:,d_index1,d_index2)=[-h*sz(d_index1,d_index2), J*sx(d_index1,d_index2),...
            id(d_index1,d_index2)];
    end
end
mpo{1}=aux;
aux2=aux;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this part correspond to right side of Eq.183 in Schollwock review
aux=zeros(3,1,d,d);
for d_index2=1:d
    for d_index1=1:d
        aux(:,1,d_index1,d_index2)=[id(d_index1,d_index2), sx(d_index1,d_index2),...
            -h*sz(d_index1,d_index2)];
    end
end
mpo{N}=aux;
aux3=aux;
%%%%%%%%%%%%%%%%%%%%%%%

%this part correspond to generic non-boundary term in Eq.182 of Schollwock
%review
for indN=2:N-1
    aux=zeros(3,3,d,d);
    aux(3,:,:,:)=aux2;
    aux(:,1,:,:)=aux3;
    mpo{indN}=aux;
end