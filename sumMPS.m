function sumMPS=sumMPS(mps1,mps2)

%evaluates <MPS1|MPS2> [hence the norm is sqrt(overlap)] for open boundary conditions
%the complexity is O(N*D^3*d)
%!!!!REMEMBER TO PROVIDE INPUT STATE THAT ARE NORMALIZED!!!!

N1=length(mps1);N2=length(mps2);
d=size(mps1{1},3);
D=size(mps1{1},2);
if N1~=N2, error('the two MPSs do not have the same number of qubits'); end

sumMPS=cell(1,N1);

sumMPS{1}=[mps1{1}, zeros(1,D,d); zeros(1,D,d), mps2{1}];
% sumMPS{1}=reshape(sumMPS{1},[1,2*D,2*d]);
for i = 2:N1-1
    sumMPS{i}=[mps1{i}, zeros(D,D,d); zeros(D,D,d), mps2{i}];
%     sumMPS{i}=reshape(sumMPS{i},[2*D,2*D,2*d]);
end
sumMPS{N1}=[mps1{N1}, zeros(D,1,d); zeros(D,1,d), mps2{N1}];
% sumMPS{N1}=reshape(sumMPS{N1},[2*D,1,2*d]);