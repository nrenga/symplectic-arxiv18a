function Ainv = g2matinv(A)

if (gfrank(A) ~= size(A,1))
    Ainv = [];
else
    [~, Ainv, ~] = g2rref(A);
end

end