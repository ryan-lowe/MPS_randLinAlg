function [Cleft]=updateCleft(Cleft,B,X,A)
% this function execute a series of contraction to evaluate 
% sequentially the overlap depicted in the right of fig. 28 pag 133 
% of Schollwock, in the case in which an mpo is also present.

% it is the efficient way of evaluating expectation values for mpo's.

if isempty(X), X=reshape(eye(size(B,3)),[1,1,2,2]); end

Cleft=contracttensors(A,3,1,Cleft,3,3); 
Cleft=contracttensors(X,4,[1,3],Cleft,4,[4,2]); %X [1,4]
Cleft=contracttensors(conj(B),3,[1,3],Cleft,4,[4,2]);

% Cleft=contracttensors(Cleft,3,3,A,3,1); 
% Cleft=contracttensors(Cleft,4,[2,4],Cleft,4,[1,3]); %X [1,4]
% Cleft=contracttensors(Cleft,4,[1,4],conj(B),4,[1,3]);