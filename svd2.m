function [U,S,V]=svd2(T)
%A slightly modified version of the built-in MATLAB svd function

[m,n]=size(T);
if m>=n, [U,S,V]=svd(T,0); else [V,S,U]=svd(T',0); end 
V=V';