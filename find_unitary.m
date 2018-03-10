function U = find_unitary(m, circuit)
% 'm' is the number of qubits in the system
% 'circuit' is an 'd x 2' cell array with each row defining a gate as 
% specified below; 'd' is the depth of the circuit
% 'U' is the overall unitary matrix for the given circuit

% Each row of the cell array must be one of the following Clifford gates:
% Gate is the desired gate; columns 1 and 2 are specifications for the gate
% that form a row of the 'circuit' cell array.

%                Gate                  |  Column 1  |   Column 2
% -----------------------------------------------------------------
% Pauli X on qubits 2,4                |     'X'    |   [2 4]
% Pauli Z on qubits 1,5                |     'Z'    |   [1 5]
% Pauli Y on qubits 1,2,5              |     'Y'    |   [1 2 5]
% Phase on qubits 1,3                  |     'P'    |   [1 3]
% Hadamard on qubits 2,4,5             |     'H'    |   [2 4 5]
% Controlled-Z on qubits 3,6           |    'CZ'    |   [3 6]
% Controlled-NOT: qubit 2 controls 1   |   'CNOT'   |   [2 1]
% Permutation (m=3): [1 2 3] -> [2 3 1]|  'Permute' |   [2 3 1]
% -----------------------------------------------------------------

% Example Circuit (m = 6 qubits): 
% U = CZ_{26} * H1 * CNOT_{12} * H2 * CNOT_{24} * H3 * CZ_{14}
% In circuit diagram the last CZ_{14} will appear first. This is because the
% operator acts on state |v> as U|v>, and so |v> goes through the last
% CZ_{14} first. Hence this is the required order for this function too.

% In this case, our specification for this function will be:
% circuit = {'CZ', [1 4]; 'H', 3; 'CNOT', [2 4]; 'H', 2; 'CNOT', [1 2]; ...
%                         'H', 1; 'CZ', [2 6]};

% Author: Narayanan Rengaswamy, Date: Mar. 1, 2018


I = eye(2);
X = [0 1; 1 0];
Z = [1 0; 0 -1];
Y = sqrt(-1) * X * Z;
P = [1 0; 0 sqrt(-1)];
H = 1/sqrt(2) * (X + Z);
e0 = [1; 0];
e1 = [0; 1];
E00 = e0 * e0';
E11 = e1 * e1';
E01 = e0 * e1';
E10 = e1 * e0';

U = eye(2^m);
for i = 1:size(circuit,1)
    gate = circuit{i,1};
    qubits = circuit{i,2}(:)';
    if (strcmpi(gate, 'X'))
        if (isempty(qubits))
            fprintf('\nPauli X Gate: Need to specify atleast one qubit!\n');
            U = [];
            return;
        end
        UX = 1;
        for j = 1:m
            if (~isempty(intersect(qubits,j)))
                UX = kron(UX, X);
            else
                UX = kron(UX, I);
            end
        end
        U = UX * U;
    elseif (strcmpi(gate, 'Z'))
        if (isempty(qubits))
            fprintf('\nPauli Z Gate: Need to specify atleast one qubit!\n');
            U = [];
            return;
        end
        UZ = 1;
        for j = 1:m
            if (~isempty(intersect(qubits,j)))
                UZ = kron(UZ, Z);
            else
                UZ = kron(UZ, I);
            end
        end
        U = UZ * U;
    elseif (strcmpi(gate, 'Y'))
        if (isempty(qubits))
            fprintf('\nPauli Y Gate: Need to specify atleast one qubit!\n');
            U = [];
            return;
        end
        UY = 1;
        for j = 1:m
            if (~isempty(intersect(qubits,j)))
                UY = kron(UY, Y);
            else
                UY = kron(UY, I);
            end
        end
        U = UY * U;
    elseif (strcmpi(gate, 'P'))
        if (isempty(qubits))
            fprintf('\nPhase Gate: Need to specify atleast one qubit!\n');
            U = [];
            return;
        end
        UP = 1;
        for j = 1:m
            if (~isempty(intersect(qubits,j)))
                UP = kron(UP, P);
            else
                UP = kron(UP, I);
            end
        end
        U = UP * U;
    elseif (strcmpi(gate, 'H'))
        if (isempty(qubits))
            fprintf('\nHadamard Gate: Need to specify atleast one qubit!\n');
            U = [];
            return;
        end
        UH = 1;
        for j = 1:m
            if (~isempty(intersect(qubits,j)))
                UH = kron(UH, H);
            else
                UH = kron(UH, I);
            end
        end
        U = UH * U;
    elseif (strcmpi(gate, 'CZ'))
        if (isempty(qubits) || length(qubits) ~= 2)
            fprintf('\nCZ Gate: Need to specify two qubits!\n');
            U = [];
            return;
        end
        UCZ1 = 1;
        UCZ2 = 1;
        for j = 1:m
            if (j == qubits(1))
                UCZ1 = kron(UCZ1, E00);
                UCZ2 = kron(UCZ2, E11);
            elseif (j == qubits(2))
                UCZ1 = kron(UCZ1, I);
                UCZ2 = kron(UCZ2, Z);
            else
                UCZ1 = kron(UCZ1, I);
                UCZ2 = kron(UCZ2, I);
            end
        end
        U = (UCZ1 + UCZ2) * U;
    elseif (strcmpi(gate, 'CNOT'))
        if (isempty(qubits) || length(qubits) ~= 2)
            fprintf('\nCNOT Gate: Need to specify two qubits!\n');
            U = [];
            return;
        end
        UCNOT1 = 1;
        UCNOT2 = 1;
        for j = 1:m
            if (j == qubits(1))
                UCNOT1 = kron(UCNOT1, E00);
                UCNOT2 = kron(UCNOT2, E11);
            elseif (j == qubits(2))
                UCNOT1 = kron(UCNOT1, I);
                UCNOT2 = kron(UCNOT2, X);
            else
                UCNOT1 = kron(UCNOT1, I);
                UCNOT2 = kron(UCNOT2, I);
            end
        end
        U = (UCNOT1 + UCNOT2) * U;
    elseif (strcmpi(gate, 'Permute'))
        desired_order = circuit{i,2}(:)';
        if (isempty(desired_order) || length(desired_order) ~= m)
            fprintf('\nPermutation: Need to specify %d qubits!\n', m);
            U = [];
            return;
        end
        current_order = 1:m;
        for j = 1:m
            % Swap j with the position of desired_order(j) in current_order
            k = find(current_order == desired_order(j));
            if (k ~= j)
                UPerm1 = 1;
                UPerm2 = 1;
                UPerm3 = 1;
                UPerm4 = 1;
                for l = 1:m
                    if (l == j)
                        UPerm1 = kron(UPerm1, E00);
                        UPerm2 = kron(UPerm2, E11);
                        UPerm3 = kron(UPerm3, E01);
                        UPerm4 = kron(UPerm4, E10);
                    elseif (l == k)
                        UPerm1 = kron(UPerm1, E00);
                        UPerm2 = kron(UPerm2, E11);
                        UPerm3 = kron(UPerm3, E10);
                        UPerm4 = kron(UPerm4, E01);
                    else
                        UPerm1 = kron(UPerm1, I);
                        UPerm2 = kron(UPerm2, I);
                        UPerm3 = kron(UPerm3, I);
                        UPerm4 = kron(UPerm4, I);
                    end
                end
                U = (UPerm1 + UPerm2 + UPerm3 + UPerm4) * U;
                current_order(1, [j k]) = current_order(1, [k j]);
            end
        end
    else  % handles all unrecognized gates
        fprintf('\nfind_unitary: Unrecognized gate encountered!\n');
        U = [];
        return;
    end
end

end