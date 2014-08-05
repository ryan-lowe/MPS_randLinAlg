function [mps]=createrandommps(N,D,d)

%Code for creating a random mps.
%INPUT: number of desired sites (qubits) N, physical dimension d (usually
%2), matrix dimension D
%OUPUT: a random MPOS

mps=cell(1,N); 
mps{1}=randn(1,D,d)/sqrt(D); 
mps{N}=randn(D,1,d)/sqrt(D); 
for i=2:(N-1)
    mps{i}=randn(D,D,d)/sqrt(D); 
end