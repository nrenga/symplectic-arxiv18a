function Ainv = gf2matinv(A)

if (gfrank(A) ~= size(A,1))
    Ainv = [];
else
    [~, Ainv, ~] = gf2rref(A);
end

end