function F_all = qfind_all_symp_mat(U, V)
% Function fo find all symplectic matrices F that satisfy
% U(1:2*m-k,:)*F = V.
% This function is specifically for the application of finding logical
% Clifford operators for stabilizer codes.
% The general version of this function is find_all_symp_mat(U, V, I, J).

% Rows of U must form a symplectic basis for \mathbb{F}_2^{2m},
% i.e., U must satisfy U*Omega*U' = Omega, where Omega = [0 I_m; I_m 0].
% The rows of U must be [Xbar; S; Zbar; Spair], where Spair are the vectors
% that complete the symplectic basis.
% Number of rows of V must be equal to (2m - k).

% Author: Narayanan Rengaswamy, Date: Feb. 20, 2018

m = size(U,2)/2;
k = 2*m - size(V,1);
tot = 2^(k*(k+1)/2);
F_all = cell(tot,1);

% Find one solution using symplectic transvections
F0 = find_symp_mat(U(1:(2*m-k), :), V);

A = mod(U * F0, 2);
Ainv = gf2dec(inv(gf(A)), 1, 3);
Basis = A([(m - k + 1):m, (2*m - k + 1):(2*m)],:);
Subspace = mod(de2bi((0:2^(2*k)-1)',2*k) * Basis, 2);
StabF0 = A((m - k + 1):m, :);
Choices = cell(k,1);

% Calculate all choices for each vector in Spair using just the conditions
% imposed by the (modified) stabilizer generators, i.e., the rows of S*F0
for i = 1:k
    % Impose symplectic inner product of 1 with the "fixed" stabilizer
    h = [zeros(1,i-1) 1 zeros(1,k-i)];
    Innpdts = mod(Subspace * fftshift(StabF0, 2)', 2);
    Choices{i,1} = Subspace(bi2de(Innpdts) == bi2de(h), :);
end

% First free vector, i.e. row 1 of Spair has 2^k choices, 
% second free vector has 2^(k-1) choices and so on...
for l = 0:(tot - 1)
    Bl = A;
    V = zeros(k,2*m);
    lbin = de2bi(l,k*(k+1)/2,'left-msb');
    v1_ind = bi2de(lbin(1,1:k),'left-msb') + 1;
    V(1,:) = Choices{1,1}(v1_ind,:);
    for i = 2:k
        vi_ind = bi2de(lbin(1,sum(k:-1:k-(i-2)) + (1:(k-(i-1)))),'left-msb') + 1;
        Innprods = mod(Choices{i,1} * fftshift(V,2)', 2);
        
        % Impose symplectic inner product of 0 with chosen free vectors
        Ch_i = Choices{i,1}(bi2de(Innprods) == 0, :);
        V(i,:) = Ch_i(vi_ind,:);
    end
    Bl((2*m - k + 1):(2*m), :) = V;
    F = mod(Ainv * Bl, 2);
    F_all{l+1,1} = mod(F0 * F, 2);
end    

end
