function t = gf2trace(x, m)
% Compute Tr(x) = x + x^2 + x^(2^2) + ... + x^(2^(m-1)), x \in GF(2^m)
% x is in the exponential format for the field GF(2^m),
% i.e., x is an integer between -1 and 2^m-2

prim_poly = gfprimdf(m);
field = gftuple([-1:2^m-2]', prim_poly);

t = -Inf;
for i = 0:(m-1)
    t = gfadd(t, mod(x*(2^i), 2^m-1), field);
end

% Convert 't' back into the actual value in GF(2) since Tr(x) \in GF(2)
if (t == -Inf)
    t = 0;
elseif (t == 0)
    t = 1;
end

end