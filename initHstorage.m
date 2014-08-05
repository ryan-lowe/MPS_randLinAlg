function [Hstorage] = initHstorage(mps,mpo,d)
%This function essentially calculates all of the intermediary contractions when
%contracting an mps with itself, with an mpo in the middle. Used
%when calculating the minimum E of an mpo.
%INPUT: an mps and an mpo
%OUTPUT: an N+1 sized structure Cstorage, with Cstorage{i} indicating the
%tensor created when the contraction mps-mpo-mps is applied for all j,
%such that i<j<=N


%this function seems to work with initial R-canonical state, 
%and sweeps procedure going first from L to R
%because of the updateCright function, and not updateCleft


%Cright is a tensor with 3-indices because we are considering also mpo; 
% if it was only overlap of states it would have 2 indices;

N = size(mpo,2);

Hstorage = cell(1,N+1);
Hstorage{1} = 1; 
Hstorage{N+1} = 1;

% if isempty(mpo)
%     mpo=cell(1,N);
%     for i=1:N
%         mpo{i}=eye(d,d);
%         mpo{i}=reshape(mpo{i},[1,1,d,d]);
%     end
% end
for j = N:-1:2
    Hstorage{j} = updateCright(Hstorage{j+1},mps{j},mpo{j},mps{j});
end