function [TG, PG] = kerdock_design(m, frobenius, circuits)
% Function to construct a unitary 2-design, using Kerdock sets, on m qubits

% The parameter 'frobenius' is a Boolean variable that indicates if the
% Frobenius automorphisms in the field GF(2^m) need to be included. If
% True (any non-zero value), then the size of the design is (N^2-1)*N*m,
% and if False (= 0), then the size of the design is (N^2-1)*N, N = 2^m.

% NOTE: To complete the unitary 2-design, each element in TG needs to be
% combined with each phaseless matrix in the Pauli group (PG) to generate N^2
% unitary elements. Therefore, the overall design size is (N^2-1)*N^3 or 
% (N^2-1)*N^3*m. The phaseless elements of the Pauli group can be obtained
% using the function pauli_group(m), which returns circuits, unitaries and
% binary representations.
% The output cell array PG is obtained as PG = pauli_group(m).

% The parameter 'circuits' is a Boolean variable that indicates if circuits
% need to be calculated for each element of the design.

% The output variable 'TG' is a cell array that consists of either 3 columns
% (if circuits = 0) or 4 columns (if circuits = 1). In the latter case, the
% last column give the quantum circuit for the element. The first column of
% 'TG' is a 5-tuple in GF(2^m) indexing the element, the second column is
% the symplectic matrix representation of the element, and the third column
% is simply a hash to compare different elements; it is just the decimal
% value of the vectorized form of the symplectic matrix.

% Author: Narayanan Rengaswamy, Date: Sep. 15, 2018

if (nargin == 1)
    frobenius = 0;
    circuits = 0;
elseif (nargin == 2)
    circuits = 0;
end

N = 2^m;
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

% Construct the matrix R representing the linear transformation for
% squaring an element x = (x_0, x_1, ..., x_{m-1}) in GF(2^m) ~= GF(2)^m,
% i.e., x^2 = [x_0, x_1, ..., x_{m-1}] * R (mod 2).

R = zeros(m);
for i = 0:(m-1)
    % i represents the i^th basis element alpha^i, where alpha is primitive
    basis_i_squared = mod(2*i, 2^m-1);
    vec_rep = field(basis_i_squared + 2, :);
    R(i+1, :) = vec_rep;
end    

PG = pauli_group(m);

if (frobenius)
    s_TG = (N + 1) * N * (N - 1) * m;  % including Frobenius automorphism
    mfro = m;
else
    s_TG = (N + 1) * N * (N - 1);
    mfro = 1;
end
if (circuits)
    TG = cell(s_TG,4);
else
    TG = cell(s_TG,3);
end    
ind = 1;

Z = zeros(m);

% Case 1: c ~= 0; ad + bc = 1
Rinv = gf2matinv(R);
for a = -1:(N-2)
    for c = 0:(N-2)
        for d = -1:(N-2)
            b = gfdiv(gfsub(gfmul(a,d,field), 0, field), c, field);
            if (b == -Inf)
                b = -1;
            end
            for i = 0:(mfro-1)
                R_i = mod(R^i, 2);
                Rinv_i = mod(Rinv^i, 2);
                Aat = mod(field_A{a+2,1}^2 * R_i, 2)';
                AbP = mod(Rinv_i * field_A{b+2,1}^2 * P, 2);
                PinvAc = mod(R_i' * gf2matinv(P) * field_A{c+2,1}^2, 2);
                Ad = mod(Rinv_i * field_A{d+2,1}^2, 2);
                TG{ind,1} = [a,b,c,d,i];
                TG{ind,2} = [Ad, AbP; PinvAc, Aat];
                TG{ind,3} = bi2de(TG{ind,2}(:)');
                if (circuits)
                    ckt = find_ckt(m,a,c,d,i,R_i,Rinv_i,field_A,P);
                    ckt2 = find_circuit(TG{ind,2});
                    if (size(ckt2,1) < size(ckt,1))
                        ckt = ckt2;
                    end
                    TG{ind,4} = ckt;
                end                
                ind = ind + 1;
            end
        end
    end
end

% Case 2: b ~= 0, c = 0
for a = 0:(N-2)
    for b = 0:(N-2)
        d = gfdiv(0, a, field);
        for i = 0:(mfro-1)
            R_i = mod(R^i, 2);
            Rinv_i = mod(Rinv^i, 2);
            Aat = mod(field_A{a+2,1}^2 * R_i, 2)';
            AbP = mod(Rinv_i * field_A{b+2,1}^2 * P, 2);
            Ad = mod(Rinv_i * field_A{d+2,1}^2, 2);
            TG{ind,1} = [a,b,-1,d,i];
            TG{ind,2} = [Ad, AbP; Z, Aat];
            TG{ind,3} = bi2de(TG{ind,2}(:)');
            if (circuits)
                e = 0;
                for ep = 0:(2^i - 1)
                    e = gfmul(e, gfmul(b,d,field), field);
                end
                Fe = [eye(m), mod(field_A{e+2,1}^2 * P,2); zeros(m), eye(m)];
                Fa = blkdiag(Ad, gf2matinv(Ad)');
                ckt = [find_circuit(Fe); find_circuit(Fa)];
                ckt2 = find_circuit(TG{ind,2});
                if (size(ckt2,1) < size(ckt,1))
                    ckt = ckt2;
                end
                TG{ind,4} = ckt;
            end
            ind = ind + 1;
        end
    end
end

% Case 3: b = 0, c = 0
for a = 0:(N-2)
    d = gfdiv(0, a, field);
    for i = 0:(mfro-1)
        R_i = mod(R^i, 2);
        Rinv_i = mod(Rinv^i, 2);
        Aat = mod(field_A{a+2,1}^2 * R_i, 2)';
        Ad = mod(Rinv_i * field_A{d+2,1}^2, 2);
        TG{ind,1} = [a,-1,-1,d,i];
        TG{ind,2} = [Ad, Z; Z, Aat];
        TG{ind,3} = bi2de(TG{ind,2}(:)');
        if (circuits)
            ckt = find_circuit(blkdiag(Ad, gf2matinv(Ad)'));
            ckt2 = find_circuit(TG{ind,2});
            if (size(ckt2,1) < size(ckt,1))
                ckt = ckt2;
            end
            TG{ind,4} = ckt;
        end        
        ind = ind + 1;
    end
end


    function ckt = find_ckt(m,a,c,d,i,R_i,Rinv_i,field_A,P)
        x = c;
        w = gfdiv(a,c,field);
        if (w == -Inf)
            w = -1;
        end
        y = 0;
        for lp = 0:(2^i - 1)
            y = gfmul(y, gfdiv(d,c,field), field);
        end
        if (y == -Inf)
            y = -1;
        end
        F1 = [eye(m), mod(field_A{y+2,1}^2 * P,2); zeros(m), eye(m)];
        F2 = blkdiag(mod(Rinv_i * gf2matinv(field_A{x+2,1})^2, 2),...
                     mod(R_i' * (field_A{x+2,1}')^2, 2));
        F3 = [zeros(m), P'; gf2matinv(P), zeros(m)];
        F4 = [eye(m), mod(field_A{w+2,1}^2 * P,2); zeros(m), eye(m)];
        ckt = [find_circuit(F1); find_circuit(F2); find_circuit(F3); ...
                   find_circuit(F4)];
    end

end