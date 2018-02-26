function circuit = find_circuit(F)
% Function to find a circuit for the given symplectic transformation F
% Uses Trung Can's algorithm to decompose F into elementary forms

% Author: Narayanan Rengaswamy, Date: Feb. 22, 2018

% CAUTION: In any consecutive sequence of CNOTs produced by this function,
%          for qubits that act as both control and target, first
%          implement the CNOTs where they act as control!

m = size(F,1)/2;
I = eye(m);
Z = zeros(m);
Omega = [Z, I; I, Z];

if (~all(all(mod(F * Omega * F', 2) == Omega)))
    fprintf('\nInvalid symplectic matrix!\n');
    circuit = [];
    return;
end

Decomp = symp_mat_decompose(F);
circuit = cell(1,2);
ckt_ind = 1;

for i = 1:length(Decomp)
    if (all(all(Decomp{i} == eye(2*m))))
        continue;
    elseif (all(all(Decomp{i} == Omega)))
        circuit{ckt_ind,1} = 'H';
        circuit{ckt_ind,2} = 1:m;
        ckt_ind = ckt_ind + 1;
        continue;
    end
    A = Decomp{i}(1:m,1:m);
    B = Decomp{i}(1:m,m+(1:m));
    C = Decomp{i}(m+(1:m),1:m);
    D = Decomp{i}(m+(1:m),m+(1:m));
    
    if (all(A(:) == I(:)) && all(C(:) == Z(:)) && all(D(:) == I(:)))
        S_ind = find(diag(B) == 1);
        for j = 1:length(S_ind)
            circuit{ckt_ind,1} = 'S';
            circuit{ckt_ind,2} = S_ind(j);
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
    elseif (all(B(:) == Z(:)) && all(C(:) == Z(:)))
        % CAUTION: For qubits that act as both control and target, first
        %          implement the CNOTs where they act as control!
        AP = A;
        if (~all(diag(A) == 1))
            valid_perm = find_valid_perm(A);
            P = I(:,valid_perm);
            AP = mod(A*P,2);
            Pinv = g2matinv(P);
        end
        for j = 1:m
            inds = setdiff(find(AP(j,:) == 1), j);
            for k = 1:length(inds)
                circuit{ckt_ind,1} = 'CNOT';
                circuit{ckt_ind,2} = [j inds(k)];  % CNOT_{j->inds(k)}
                ckt_ind = ckt_ind + 1;
            end
        end    
        if (~all(diag(A) == 1))
            circuit{ckt_ind,1} = 'Permute';
            circuit{ckt_ind,2} = (1:m)*Pinv;
            ckt_ind = ckt_ind + 1;
        end
        
        % Another possibility: naively implement the transform given by A
%         ancilla = m;
%         for j = 1:m
%             inds = setdiff(find(A(:,j) == 1), j);
%             if (A(j,j) == 1)
%                 for k = 1:length(inds)
%                     circuit{ckt_ind,1} = 'CNOT';
%                     circuit{ckt_ind,2} = [inds(k) j];  % CNOT_{inds(k)->j}
%                     ckt_ind = ckt_ind + 1;
%                 end
%             else
%                 ancilla = ancilla + 1;
%                 circuit{ckt_ind,1} = 'Ancilla in |0>';
%                 circuit{ckt_ind,2} = ancilla;
%                 ckt_ind = ckt_ind + 1;
%                 circuit{ckt_ind,1} = 'CNOTs';
%                 circuit{ckt_ind,2} = [inds', ancilla];
%                 ckt_ind = ckt_ind + 1;
%             end
%         end
%         if (ancilla > m)
%             diagzs = find(diag(A) == 0);
%             circuit{ckt_ind,1} = 'Swap with ancilla';
%             circuit{ckt_ind,2} = strcat(mat2str(diagzs),' = ',mat2str((m+1):ancilla));
%             ckt_ind = ckt_ind + 1;
%         end
    else
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
            
    function valid_perm = find_valid_perm(A)
        Indices = cell(size(A,1),1);
        for iter = 1:size(A,1)
            Indices{iter,1} = find(A(iter,:) == 1);
        end
        Set = cartprod(Indices{:,1});
        validP = find(sum(abs(bsxfun(@minus,sort(Set,2),1:size(A,1))), 2) == 0, 1, 'first');
        valid_perm = Set(validP,:);
    end

end