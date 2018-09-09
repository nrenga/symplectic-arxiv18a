function circuit = find_circuit(F)
% Function to find a circuit for the given symplectic transformation F
% Uses Trung Can's algorithm to decompose F into elementary forms

% Author: Narayanan Rengaswamy, Date: Mar. 1, 2018

m = size(F,1)/2;
I = eye(m);
O = zeros(m);
Omega = [O, I; I, O];

if (~all(all(mod(F * Omega * F', 2) == Omega)))
    fprintf('\nInvalid symplectic matrix!\n');
    circuit = [];
    return;
end

if (all(all(F == eye(2*m))))
    circuit = [];
    return;
end

Decomp = symp_mat_decompose(F);   % Trung Can's decomposition
circuit = cell(1,2);
ckt_ind = 1;

for i = 1:size(Decomp,1)
    if (all(all(Decomp{i} == eye(2*m))))
        continue;
    elseif (all(all(Decomp{i} == Omega)))
        % Transversal Hadamard
        circuit{ckt_ind,1} = 'H';
        circuit{ckt_ind,2} = 1:m;
        ckt_ind = ckt_ind + 1;
        continue;
    end
    A = Decomp{i}(1:m,1:m);
    B = Decomp{i}(1:m,m+(1:m));
    C = Decomp{i}(m+(1:m),1:m);
    D = Decomp{i}(m+(1:m),m+(1:m));
    
    if (all(A(:) == I(:)) && all(C(:) == O(:)) && all(D(:) == I(:)))
        % CZs and Phase
        P_ind = find(diag(B) == 1)';
        if (~isempty(P_ind))
            circuit{ckt_ind,1} = 'P';
            circuit{ckt_ind,2} = P_ind;
            ckt_ind = ckt_ind + 1;
        end
        
        % Clear diagonal entries, extract upper triangular part as B = B'
        B = triu(mod(B + diag(diag(B)), 2));
        for j = 1:m
            CZ_ind = find(B(j,:) == 1);
            for k = 1:length(CZ_ind)
                circuit{ckt_ind,1} = 'CZ';
                circuit{ckt_ind,2} = [j CZ_ind(k)];
                ckt_ind = ckt_ind + 1;
            end
        end
    elseif (all(B(:) == O(:)) && all(C(:) == O(:)))
        % CNOTs and Permutations
        % CAUTION: For qubits that act as both control and target, first
        %          implement the CNOTs where they act as control!
        %          To ensure this, we use LU decomposition over GF(2).
        [L, U, P] = gf2lu(A);   % P' * L * U = A (P is a permutation matrix)
        if (~all(P(:) == I(:)))
            circuit{ckt_ind,1} = 'Permute';
            circuit{ckt_ind,2} = (1:m)*P';
            ckt_ind = ckt_ind + 1;
        end
        for j = 1:m
            inds = setdiff(find(L(j,:) == 1), j);
            for k = 1:length(inds)
                circuit{ckt_ind,1} = 'CNOT';
                circuit{ckt_ind,2} = [j inds(k)];  % CNOT_{j->inds(k)}
                ckt_ind = ckt_ind + 1;
            end
        end    
        for j = m:-1:1
            inds = setdiff(find(U(j,:) == 1), j);
            for k = 1:length(inds)
                circuit{ckt_ind,1} = 'CNOT';
                circuit{ckt_ind,2} = [j inds(k)];  % CNOT_{j->inds(k)}
                ckt_ind = ckt_ind + 1;
            end
        end    
    else
        % Partial Hadamards
        k = m - sum(diag(A));
        Uk = blkdiag(eye(k), zeros(m-k));
        Lmk = blkdiag(zeros(k), eye(m-k));
        if (all(A(:) == Lmk(:)) && all(B(:) == Uk(:)) && ...
            all(C(:) == Uk(:)) && all(D(:) == Lmk(:)))
            circuit{ckt_ind,1} = 'H';
            circuit{ckt_ind,2} = 1:k;
            ckt_ind = ckt_ind + 1;
        else
            fprintf('\nUnknown elementary symplectic form!\n');
            circuit = [];
            break;
        end
    end
end    
            
end