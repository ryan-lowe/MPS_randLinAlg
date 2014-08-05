function [Cstorage]=initCstorage(mpsB,mpoX,mpsA,N)
%This function essentially calculates all of the intermediary contractions when
%contracting 2 mps's, mpsB and mpsA, with an mpo mpoX in the middle. Used
%when calculating the minimum E of an mpo.
%INPUT: mps's to contract, mpsB and mpsA, mpo mpoX, and number of qubits /
%indices N
%OUTPUT: an N+1 sized structure Cstorage, with Cstorage{i} indicating the
%tensor created when the contraction mpsB-mpoX-mpsA is applied for all j,
%such that i<j<=N


%this function seems to work with initial R-canonical state, 
%and sweeps procedure going first from L to R
%because of the updateCright function, and not updateCleft


%Cright is a tensor with 3-indices because we are considering also mpo; 
% if it was only overlap of states it would have 2 indices;
Cstorage=cell(1,N+1); 
Cstorage{1}=1; 
Cstorage{N+1}=1;
for i=N:-1:2
    if isempty(mpoX), X=[]; else X=mpoX{i}; end
    Cstorage{i}=updateCright(Cstorage{i+1},mpsB{i},X,mpsA{i}); 
end