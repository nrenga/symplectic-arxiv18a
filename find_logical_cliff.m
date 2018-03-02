function F_all = find_logical_cliff(S, Xbar, Zbar, circuit, Snorm, no_of_solns)
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

% Each row in the cell array 'circuit' consists of two columns. The first
% column is the logical gate and the second column is the vector of logical
% qubits the gate is applied on. Here are the formats.

%            Logical Gate              |  Column 1  |  Column 2
% -----------------------------------------------------------------
% Phase on (logical) qubits 1,3,5      |     'S'    |   [1 3 5]
% Targeted Hadamard on qubits 2,4      |     'H'    |   [2 4]
% Controlled-Z on qubits 3,6           |    'CZ'    |   [3 6]
% Controlled-NOT: qubit 2 controls 1   |   'CNOT'   |   [2 1]
% -----------------------------------------------------------------

% Example Circuit (m-k = 4 logical qubits): 
% U = CZ_{24} * H1 * CNOT_{12} * H2 * CNOT_{13} * H3 * CZ_{14}
% In circuit diagram the last CZ_{14} will appear first. This is because the
% operator acts on state |v> as U|v>, and so |v> goes through the last
% CZ_{14} first. Hence this is the required order for this function too.

% In this case, our specification for this function will be:
% circuit = {'CZ', [1 4]; 'H', 3; 'CNOT', [1 3]; 'H', 2; 'CNOT', [1 2]; ...
%                         'H', 1; 'CZ', [2 4]};

% Author: Narayanan Rengaswamy, Date: Mar. 1, 2018

if (nargin <= 4)
    Snorm = S;   % Assume physical operator for 'gate' must centralize S
    no_of_solns = 1;  % Produce only one solution
end
if (nargin <= 5)
    no_of_solns = 1;
end
if (isempty(Snorm))
    Snorm = S;
end

[k, m] = size(S);
m = m/2;
tot = 2^(k*(k+1)/2);
F_all = cell(1,3);

H = [Xbar; Snorm; Zbar];

for i = 1:size(circuit,1)
    gate = circuit{i,1};
    qubits = circuit{i,2};
    if (any(qubits > m-k))
        fprintf('\nLogical qubit index exceeds %d...!\n', m-k);
        return;
    end
    if (strcmpi(gate, 'S'))
        if (isempty(qubits))
            fprintf('\nPhase Gate: Need to specify at least one qubit!\n');
            return;
        else
            H(qubits,:) = mod(H(qubits,:) + H(m+qubits,:), 2);
        end
    elseif (strcmpi(gate, 'CZ'))
        if (isempty(qubits) || length(qubits) ~= 2)
            fprintf('\nCZ Gate: Need to specify two qubits!\n');
            return;
        else
            H(qubits(1),:) = mod(H(qubits(1),:) + H(m+qubits(2),:), 2);
            H(qubits(2),:) = mod(H(qubits(2),:) + H(m+qubits(1),:), 2);
        end
    elseif (strcmpi(gate, 'CNOT'))
        if (isempty(qubits) || length(qubits) ~= 2)
            fprintf('\nCNOT Gate: Need to specify two qubits!\n');
            return;
        else
            H(qubits(1),:) = mod(H(qubits(1),:) + H(qubits(2),:), 2);
            H(m+qubits(2),:) = mod(H(m+qubits(2),:) + H(m+qubits(1),:), 2);
        end
    elseif (strcmpi(gate, 'H'))
        if (isempty(qubits))
            fprintf('\nHadamard Gate: Need to specify atleast one qubit!\n');
            return;
        else
            H([qubits, m+qubits],:) = H([m+qubits, qubits],:);
        end
    else
        fprintf('\nUnknown gate!\n');
        return;
    end
end

% Need to find symplectic matrices that map the first 2m-k rows of H to U.
% This corresponds to the action g*E(a,b)*g^{\dagger} = E([a,b]*F_g).

% First we need to complete a symplectic basis for \mathbb{F}_2^{2m}
U = [Xbar; S; Zbar];
for i = 1:k
    h = zeros(2*m-k+(i-1),1);
    h(m-k+i) = 1;
    U(2*m-k+i,:) = gflineq(fftshift(U,2),h)';
end

if (no_of_solns == 1)
    F_all = cell(1,3);
    F_all{1,1} = find_symp_mat(U(1:(2*m-k), :), H);
else
    fprintf('\nCalculating all %d solutions! Please be patient...\n', tot);
    F_all = cell(tot,3);
    F_all(:,1) = qfind_all_symp_mat(U, H);

    % Can also use the general algorithm as below, 
    % but need to specify I = 1:m and J = 1:m-k
    % F_all(:,1) = find_all_symp_mat(U, H, 1:m, 1:m-k);
end


for i = 1:size(F_all,1)
    F_all{i,2} = find_circuit(F_all{i,1});
    F_all{i,3} = size(F_all{i,2},1);   % Circuit depth
    
    % Check signs on conjugation with stabilizer generators, logical Paulis
    Oper = find_unitary(m, F_all{i,2});
    v = zeros(2*m-k, 1);
    for j = 1:(2*m-k)
        iota = sqrt(-1);
        h = U(j, 1:m) + iota * U(j, m+(1:m));
        hx = find(h == 1);
        hz = find(h == iota);
        hy = find(h == 1 + iota);
        h_ckt = cell(1,2);
        ind = 1;
        if (~isempty(hx))
            h_ckt(ind,:) = {'X', hx};
            ind = ind + 1;
        end
        if (~isempty(hz))
            h_ckt(ind,:) = {'Z', hz};
            ind = ind + 1;
        end
        if (~isempty(hy))
            h_ckt(ind,:) = {'Y', hy};
        end
        h_u = find_unitary(m, h_ckt);
        
        hnew = H(j, 1:m) + iota * H(j, m+(1:m));
        hnewx = find(hnew == 1);
        hnewz = find(hnew == iota);
        hnewy = find(hnew == 1 + iota);
        hnew_ckt = cell(1,2);
        ind = 1;
        if (~isempty(hnewx))
            hnew_ckt(ind,:) = {'X', hnewx};
            ind = ind + 1;
        end
        if (~isempty(hnewz))
            hnew_ckt(ind,:) = {'Z', hnewz};
            ind = ind + 1;
        end
        if (~isempty(hnewy))
            hnew_ckt(ind,:) = {'Y', hnewy};
        end
        hnew_u = find_unitary(m, hnew_ckt);
        
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
        choices = fftshift(gflineq_all(U(1:2*m-k,:), v)',2);
        choices = choices(:,1:m) + sqrt(-1)*choices(:,m+(1:m));
        [~, cheap_ind] = min(sum(choices ~= 0, 2));
        x = choices(cheap_ind, :);
        ckt_ind = F_all{i,3} + 1;
        if (any(x == 1))
            F_all{i,2}(ckt_ind, :) = {'X', find(x == 1)};
            ckt_ind = ckt_ind + 1;
        end
        if (any(x == iota))
            F_all{i,2}(ckt_ind, :) = {'Z', find(x == sqrt(-1))};
            ckt_ind = ckt_ind + 1;
        end
        if (any(x == 1 + iota))
            F_all{i,2}(ckt_ind, :) = {'Y', find(x == 1+sqrt(-1))};
        end
        F_all{i,3} = F_all{i,3} + 1;  % As 1-qubit gates add only depth 1
    end
end

end
