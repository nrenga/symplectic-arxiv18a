function F_all = find_all_symp_mat(U, V, I, J)
% Function fo find all symplectic matrices F that satisfy 
% U([I,m+J],:)*F = V.
% Rows of U must form a symplectic basis for \mathbb{F}_2^{2m},
% i.e., U must satisfy U*Omega*U' = Omega, where Omega = [0 I_m; I_m 0].
% Number of rows of V must be equal to (length(I) + length(J)).

% Author: Narayanan Rengaswamy, Date: Feb. 20, 2018

m = size(U,2)/2;
Omega = [zeros(m), eye(m); eye(m), zeros(m)];
if (~all(all(mod(U*Omega*U',2) == Omega)))
    F_all = cell(1,1);
    fprintf('\nInvalid matrix U in function find_all_symp_mat.m\n');
    return;
end

I = I(:)';
J = J(:)';
Ibar = setdiff(1:m,I);
Jbar = setdiff(1:m,J);
alpha = length(Ibar) + length(Jbar); 
tot = 2^(alpha*(alpha+1)/2);
F_all = cell(tot,1);

% Find one solution using symplectic transvections
F0 = find_symp_mat(U([I, m+J], :), V);

A = mod(U * F0, 2);
Ainv = gf2matinv(A);  % can replace by any fn. that inverts A
IbJb = union(Ibar,Jbar);
Basis = A([IbJb, m+IbJb],:);
Subspace = mod(de2bi((0:2^(2*length(IbJb))-1)',2*length(IbJb)) * Basis, 2);

% Collect indices of free vectors in the top and bottom halves of Basis
% Note: these are now row indices of Basis, not row indices of A!!
[~, Basis_fixed_I, ~] = intersect(IbJb,I);  % = intersect(I,Jbar)
[~, Basis_fixed_J, ~] = intersect(IbJb,J);  % = intersect(Ibar,J)
Basis_fixed = [Basis_fixed_I, length(IbJb) + Basis_fixed_J];
Basis_free = setdiff(1:2*length(IbJb), Basis_fixed);

% length(Basis) = 2*length(IbJb)
% length(Basis_fixed) = 2*length(IbJb) - alpha
% length(Basis_free) = alpha

Choices = cell(alpha,1);

% Calculate all choices for each free vector using just conditions imposed
% by the fixed vectors in Basis (or equivalently in A)
for i = 1:alpha
    ind = Basis_free(i);
    h = zeros(1,length(Basis_fixed));
    
    % Impose symplectic inner product of 1 with the "fixed" symplectic pair
    if (i <= length(Ibar))
        h(Basis_fixed == length(IbJb) + ind) = 1;
    else
        h(Basis_fixed == ind - length(IbJb)) = 1;
    end    
    
    Innpdts = mod(Subspace * fftshift(Basis(Basis_fixed,:), 2)', 2);
    Choices{i,1} = Subspace(bi2de(Innpdts) == bi2de(h), :);
end

% First free vector has 2^(alpha) choices, second has 2^(alpha-1) and so on
for l = 0:(tot - 1)
    Bl = A;
    W = zeros(alpha,2*m);   % Rows are choices made for free vectors
                            % W(i,:) corresponds to Basis(Basis_free(i),:)
    lbin = de2bi(l,alpha*(alpha+1)/2,'left-msb');
    v1_ind = bi2de(lbin(1,1:alpha),'left-msb') + 1;
    W(1,:) = Choices{1,1}(v1_ind,:);
    for i = 2:alpha
        vi_ind = bi2de(lbin(1,sum(alpha:-1:alpha-(i-2)) + (1:(alpha-(i-1)))),'left-msb') + 1;
        Innprods = mod(Choices{i,1} * fftshift(W,2)', 2);
        
        % Impose symplectic inner product of 0 with chosen free vectors
        h = zeros(1,alpha);
        % Handle case when Basis contains a symplectic pair of free vectors
        if (i > length(Ibar))
            h(Basis_free == Basis_free(i) - length(IbJb)) = 1;
        end
        
        Ch_i = Choices{i,1}(bi2de(Innprods) == bi2de(h), :);
        W(i,:) = Ch_i(vi_ind,:);
    end
    Bl([Ibar, m+Jbar], :) = W;
    F = mod(Ainv * Bl, 2);
    F_all{l+1,1} = mod(F0 * F, 2);
end    

end
