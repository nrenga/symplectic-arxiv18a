function Ainv = gf2matinv(A)
% Function to calculate the inverse of a matrix over GF(2).

% Author: Narayanan Rengaswamy, Date: Mar. 1, 2018

if (gfrank(A) ~= size(A,1))
    Ainv = [];
else
    [~, Ainv, ~] = gf2rref(A);
end

end