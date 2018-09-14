function [DG, field_A, P, F] = DGSet(m, r)
% Return the Delsarte-Goethals set DG(m,r) of symmetric m x m matrices P_z 
% in GF(2), z \in GF(2^m)^{r+1}, that satisfy rank(P_z) >= m-2r for z ~= 0

% Author: Narayanan Rengaswamy, Date: Sep. 14, 2018

prim_poly = gfprimdf(m);
field = gftuple([-1:2^m-2]', prim_poly);

% Construct a matrix A that generates the field GF(2^m). This matrix must
% represent right multiplication by the primitive element alpha in GF(2^m),
% i.e., if x = (x_0,x_1,...,x_{m-1}) then x*alpha = [x_0,...,x_{m-1}] * A.
% The matrix A is essentially the companion matrix for prim_poly but
% written in a form that is compatible with the above representation of x.

I = eye(m);
A = I(2:m, :);
A(m,:) = prim_poly(1,1:m);
field_A = cell(2^m,1);
field_A{1} = zeros(m);
field_A{2} = eye(m);
for i = 1:(2^m-2)
    field_A{i+2} = mod(A^i,2);
end

% Construct the symmetric matrix P representing the symmetric bilinear form
% beta(x,y) = Tr(xy), where Tr(x) = x + x^2 + x^(2^2) + ... + x^(2^(m-1)), 
% for x \in GF(2^m), i.e., beta(x,y) = Tr(xy) = x*P*y'; x,y \in GF(2)^m

P = zeros(m);
for i = 0:(m-1)
    for j = i:(m-1)
        P(i+1, j+1) = gf2trace(mod(i+j, 2^m-1), m);
        P(j+1, i+1) = P(i+1, j+1);
    end
end

% Construct the matrix F representing the linear transformation for
% squaring an element x = (x_0, x_1, ..., x_{m-1}) in GF(2^m) ~= GF(2)^m,
% i.e., x^2 = [x_0, x_1, ..., x_{m-1}] * F (mod 2).

F = zeros(m);
for i = 0:(m-1)
    % i represents the i^th basis element alpha^i, where alpha is primitive
    basis_i_squared = mod(2*i, 2^m-1);
    vec_rep = field(basis_i_squared + 2, :);
    F(i+1, :) = vec_rep;
end    

% Given a vector z = (z_0, z_1, ..., z_r) \in GF(2^m)^{r+1}, the
% corresponding Delsarte-Goethals matrix is defined as
% P_z = A_{z_0}*P + \sum_{i=1}^{r} [A_{z_i}*P*(F^i)' + (F^i)*P*A_{z_i}']

total_Ps = (2^m)^(r+1);
DG = cell(total_Ps, 2);
for i = 0:(total_Ps-1)
    z = de2bi(i, r+1, 2^m) - 1;  % Vector with entries in [-1:(2^m-2)]
    DG{i+1,1} = mod(field_A{z(1) + 2} * P, 2);
    for j = 1:r
        term = mod(field_A{z(j+1) + 2} * P * (F^j)', 2);
        term = term + mod(F^j * P * field_A{z(j+1) + 2}', 2);
        DG{i+1,1} = mod(DG{i+1} + term, 2);
    end
    DG{i+1,2} = bi2de(DG{i+1,1}(:)');  % just a hash for the matrix P_z
end


end