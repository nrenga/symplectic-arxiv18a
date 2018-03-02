function inn = symp_inn_pdt(X, Y)
% Function to compute the symplectic inner product over GF(2) between the
% corresponding rows of X and Y.

% Author: Narayanan Rengaswamy, Date: Mar. 1, 2018

inn = mod(sum(X.*fftshift(Y,2), 2), 2);

end