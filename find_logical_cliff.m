function F_all = find_logical_cliff(S, Xbar, Zbar, gate, qubits, Snorm)
% Function to find all symplectic matrices for a logical Clifford operator.
% Snorm specifies if one wants to only normalize the stabilizer S,
% in which case the rows of S (i.e., the stabilizer generators) are mapped
% to the rows of Snorm.
% If Snorm = S, then it specifies that the stabilizer must be centralized.

% Author: Narayanan Rengaswamy, Date: Feb. 20, 2018

%            Logical Gate              |   'gate'   | 'qubits'
% -----------------------------------------------------------------
% Phase on qubit 1                     |     'S'    |   [1]
% Targeted Hadamard on qubit 4         |     'H'    |   [4]
% Controlled-Z on qubits 3,6           |    'CZ'    |   [3 6]
% Controlled-NOT: qubit 2 controls 1   |   'CNOT'   |   [2 1]
% -----------------------------------------------------------------

if (nargin == 5)
    Snorm = S;   % Assume physical operator for 'gate' must centralize S
end

[k, m] = size(S);
m = m/2;
tot = 2^(k*(k+1)/2);
F_all = cell(tot,2);

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

% Need to complete a symplectic basis for \mathbb{F}_2^{2m}
U = [Xbar; S; Zbar];        
for i = 1:k
    h = zeros(2*m-k+(i-1),1);
    h(m-k+i) = 1;
    U(2*m-k+i,:) = gflineq(fftshift(U,2),h)';
end

F_all(:,1) = qfind_all_symp_mat(U, H);

% Can also use the general algorithm as below, 
% but need to specify I = 1:m and J = 1:m-k
% F_all(:,1) = find_all_symp_mat(U, H, 1:m, 1:m-k);

for i = 1:tot
    F_all{i,2} = find_circuit(F_all{i,1});
end

end
