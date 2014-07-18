function mpsout =applympo(mps,mpo)

%evaluates <MPS1|MPS2> [hence the norm is sqrt(overlap)] for open boundary conditions
%the complexity is O(N*D^3*d)
%!!!!REMEMBER TO PROVIDE INPUT STATE THAT ARE NORMALIZED!!!!

N1=length(mps);N2=length(mpo);

if N1~=N2, error('the MPS and MPO do not have the same number of qubits'); end
N=N1;
mpsout=cell(N);

for i=N:-1:1
    [mpoD1,mpoD2,~,d]=size(mpo{indN});
    [mpsD1,mpsD2,~]=size(mps{indN});
    
    [mpsout{i} ~] = contracttensors(mpo{i},4,4,mps{i},3,3);  
    
    mpsout{indN}=permute(mpsout{indN},[1 4 2 5 3]);
    mpsout{indN}=reshape(mpsout{indN},[mpoD1*mpsD1,mpoD2*mpsD2,d]);
end


