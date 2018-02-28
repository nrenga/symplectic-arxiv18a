function [L, U, P] = gf2lu(A)
% Function to perform LU decomposition of A such that P*A = L*U,
% where A is an m x n binary matrix, P is a permutation matrix, 
% L is an m x m lower-triangular matrix and U is an m x n matrix in
% row-echelon form (it is upper triangular).

% See Algorithm 21.1 in "Numerical Linear Algebra" by Trefethen and Bau

% Author: Narayanan Rengaswamy, Date: Feb. 27, 2018

m = size(A,1);
U = A;
L = eye(m);
P = eye(m);

for k = 1:(m-1)
    i = find(U(k:m,k) == 1, 1, 'first') + k - 1;
    U([k i],k:m) = U([i k],k:m);
    L([k i],1:k-1) = L([i k],1:k-1);
    P([k i],:) = P([i k],:);
    for j = k+1:m
        L(j,k) = U(j,k);
        U(j,k:m) = mod(U(j,k:m) - L(j,k) * U(k,k:m), 2);
    end
end

end