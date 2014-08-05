function [mpo]=createrandommpo(N,D,d)

%Ryan's code for creating a random mpo.
%INPUT: number of desired sites (qubits) N, physical dimension d (usually
%2), matrix dimension D
%OUPUT: a random MPO

mpo=cell(1,N); 
mpo{1}=randn(1,D,d,d)/sqrt(D); 
mpo{N}=randn(D,1,d,d)/sqrt(D); 
for i=2:(N-1)
    mpo{i}=randn(D,D,d,d)/sqrt(D); 
end