% This is a modified version of matlab's building rref which calculates
% row-reduced echelon form in gf(2).  Useful for linear codes.
% Tolerance was removed because yolo, and because all values
% should only be 0 or 1.  @benathon

function [Ar, M, N, k] = g2rref(A)
%G2RREF   Reduced row echelon form in gf(2).
%   R = RREF(A) produces the reduced row echelon form of A in gf(2).
%
%   Class support for input A:
%      float: with values 0 or 1
%   Copyright 1984-2005 The MathWorks, Inc. 
%   $Revision: 5.9.4.3 $  $Date: 2006/01/18 21:58:54 $

% Modified to return the matrix M of row operations on A, i.e., Ar = M*A,
% and the matrix N of column operations which if applied to Ar results in
% a matrix of the form [I_k, 0; 0, 0] for the first m columns, 
% where k is the gf(2) rank of A, 
% i.e., (Ar*N)_{1:m,1:m} = (M*A*N)_{1:m,1:m} = [I_k, 0; 0, 0].
% For a square matrix A, Ar*N = M*A*N = [I_k, 0; 0, 0].

% By Narayanan Rengaswamy. Date: Feb. 22, 2018

[m,n] = size(A);
M = eye(m);
Ar = [A, M];
nr = size(Ar, 2);

% Loop over the entire matrix.
i = 1;
j = 1;

while (i <= m) && (j <= n)
   % Find value and index of largest element in the remainder of column j.
   k = find(Ar(i:m,j),1) + i - 1;

   % Swap i-th and k-th rows.
   Ar([i k],j:nr) = Ar([k i],j:nr);
   
   % Save the right hand side of the pivot row
   aijn = Ar(i,j:nr);
   
   % Column we're looking at
   col = Ar(1:m,j);
   
   % Never Xor the pivot row against itself
   col(i) = 0;
   
   % This builds an matrix of bits to flip
   flip = col*aijn;
   
   % Xor the right hand side of the pivot row with all the other rows
   Ar(1:m,j:nr) = xor( Ar(1:m,j:nr), flip );

   i = i + 1;
   j = j + 1;
end
M = Ar(1:m,(n+1):nr);
Ar = Ar(1:m,1:n);
N = [mod(eye(m,n) + Ar + [diag(diag(Ar)), zeros(m, n-m)], 2);
                zeros(n-m,m), eye(n-m)];
k = sum(diag(Ar));