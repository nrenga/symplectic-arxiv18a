function F_all = find_logical_cliff(S, Xbar, Zbar, gate, qubits, Snorm, no_of_solns)
% Function to find all symplectic matrices for a logical Clifford operator.

% 'Snorm' specifies if one wants to only normalize the stabilizer S,
% in which case the rows of S (i.e., the stabilizer generators) are mapped
% to the rows of Snorm.
% If Snorm = S, then it specifies that the stabilizer must be centralized.
% Default is S.

% 'no_of_solns' indicates the desired number of solutions: 
% no_of_solns = 1 produces only one solution, which need not be the
% cheapest (in terms of circuit depth), and any other value produces all
% solutions. Default is 1.

% Author: Narayanan Rengaswamy, Date: Feb. 28, 2018

%            Logical Gate              |   'gate'   | 'qubits'
% -----------------------------------------------------------------
% Phase on qubit 1                     |     'S'    |   [1]
% Targeted Hadamard on qubit 4         |     'H'    |   [4]
% Controlled-Z on qubits 3,6           |    'CZ'    |   [3 6]
% Controlled-NOT: qubit 2 controls 1   |   'CNOT'   |   [2 1]
% -----------------------------------------------------------------

if (nargin <= 5)
    Snorm = S;   % Assume physical operator for 'gate' must centralize S
    no_of_solns = 1;  % Produce only one solution
end
if (nargin <= 6)
    no_of_solns = 1;
end
if (isempty(Snorm))
    Snorm = S;
end

[k, m] = size(S);
m = m/2;
tot = 2^(k*(k+1)/2);

G = [Xbar; Snorm; Zbar];

if (strcmpi(gate, 'S'))
    if (isempty(qubits) || length(qubits) ~= 1)
        fprintf('\nNeed to specify one qubit!\n');
        return;
    else
        H = G;
        H(qubits,:) = mod(Xbar(qubits,:) + Zbar(qubits,:), 2);
    end
elseif (strcmpi(gate, 'CZ'))
    if (isempty(qubits) || length(qubits) ~= 2)
        fprintf('\nNeed to specify two qubits!\n');
        return;
    else
        H = G;
        H(qubits(1),:) = mod(Xbar(qubits(1),:) + Zbar(qubits(2),:), 2);
        H(qubits(2),:) = mod(Xbar(qubits(2),:) + Zbar(qubits(1),:), 2);
    end
elseif (strcmpi(gate, 'CNOT'))
    if (isempty(qubits) || length(qubits) ~= 2)
        fprintf('\nNeed to specify two qubits!\n');
        return;
    else
        H = G;
        H(qubits(1),:) = mod(Xbar(qubits(1),:) + Xbar(qubits(2),:), 2);
        H(m+qubits(2),:) = mod(Zbar(qubits(2),:) + Zbar(qubits(1),:), 2);
    end
elseif (strcmpi(gate, 'H'))
    if (isempty(qubits) || length(qubits) ~= 1)
        fprintf('\nNeed to specify one qubit!\n');
        return;
    else
        H = G;
        H(qubits,:) = Zbar(qubits,:);
        H(m+qubits,:) = Xbar(qubits,:);
    end
else
    fprintf('\nUnknown gate!\n');
    return;
end

% Need to find symplectic matrices that map the first 2m-k rows of H to U.
% This corresponds to the action g E(a,b) g^{\dagger} = E[(a,b) F_g].
% However our elementaty transformations only satisfy the equation
% g^{\dagger} E(a,b) g = E[(a,b) F_g]. Hence we will find solutions F that
% map rows of U to H and then take the inverse of each such solution.

% First we need to complete a symplectic basis for \mathbb{F}_2^{2m}
U = H;
H = [Xbar; S; Zbar];        
for i = 1:k
    h = zeros(2*m-k+(i-1),1);
    h(m-k+i) = 1;
    U(2*m-k+i,:) = gflineq(fftshift(U,2),h)';
end

if (no_of_solns == 1)
    F_all = cell(1,3);
    F_all{1,1} = find_symp_mat(U(1:(2*m-k), :), H);
else
    F_all = cell(tot,3);
    F_all(:,1) = qfind_all_symp_mat(U, H);

    % Can also use the general algorithm as below, 
    % but need to specify I = 1:m and J = 1:m-k
    % F_all(:,1) = find_all_symp_mat(U, H, 1:m, 1:m-k);
end


for i = 1:size(F_all,1)
    F_all{i,1} = gf2matinv(F_all{i,1});
    F_all{i,2} = find_circuit(F_all{i,1});
    F_all{i,3} = size(F_all{i,2},1);   % Circuit depth
    
    % Check signs on conjugation with stabilizer generators
    Oper = find_unitary(m, F_all{i,2});
    v = zeros(2*m-k, 1);
    for j = 1:(2*m-k)
        hx = H(j, 1:m);
        hz = H(j, m+(1:m));
        iota = sqrt(-1);
        h_u = find_unitary(m, {'X', find(hx + iota*hz == 1); ...
                                  'Z', find(hx + iota*hz == iota); ...
                                     'Y', find(hx + iota*hz == 1 + iota)});
        hnewx = U(j, 1:m);
        hnewz = U(j, m+(1:m));
        hnew_u = find_unitary(m, {'X', find(hnewx + iota*hnewz == 1); ...
                                    'Z', find(hnewx + iota*hnewz == iota); ...
                                     'Y', find(hnewx + iota*hnewz == 1 + iota)});
        h_u_new = Oper * h_u * Oper';
        if (norm(h_u_new(:) - (-hnew_u(:)), 'fro') < 1e-10)
            v(j) = 1;
        else
            if (~(norm(h_u_new(:) - hnew_u(:), 'fro') < 1e-10))
                fprintf('\nSomething is wrong for stabilizer %d!!\n', j);
            end
        end
    end
    if (any(v == 1))
        choices = fftshift(gflineq_all(H, v)',2);
        choices = choices(:,1:m) + sqrt(-1)*choices(:,m+(1:m));
        [~, cheap_ind] = min(sum(choices ~= 0, 2));
        x = choices(cheap_ind, :);
        ckt_ind = F_all{i,3} + 1;
        if (any(x == 1))
            F_all{i,2}(ckt_ind, :) = {'X', find(x == 1)};
            ckt_ind = ckt_ind + 1;
        end
        if (any(x == sqrt(-1)))
            F_all{i,2}(ckt_ind, :) = {'Z', find(x == sqrt(-1))};
            ckt_ind = ckt_ind + 1;
        end
        if (any(x == 1 + sqrt(-1)))
            F_all{i,2}(ckt_ind, :) = {'Y', find(x == 1+sqrt(-1))};
        end
        F_all{i,3} = F_all{i,3} + 1;  % As 1-qubit gates add only depth 1
    end
end

end
