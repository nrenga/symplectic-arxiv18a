function op_out = calculate_conj(m, op_in, circuit)
% Function to compute the effect of a Clifford circuit on an input Pauli
% operator under conjugation, i.e. technically this computes
% op_out = circuit * op_in * circuit'.

% 'op_in' is the input Pauli operator specified as a 1 x 2 cell array. The
% first column is a string with characters X,Y and Z, and the second column
% is a vector listing their corresponding qubit indices. 'op_out' will be
% in the same format but the first character in its first column could be
% '-' or 'i' or 'j' representing a negative sign or sqrt(-1) or -sqrt(-1)
% obtained as a result of conjugation, respectively. In that case the first
% entry in its 2nd column will be -1 or sqrt(-1) or -sqrt(-1), respectively.

% Each row of the cell array 'circuit' must be one of the following gates:
% Gate is the desired gate; columns 1 and 2 are specifications for the gate
% that form a row of the 'circuit' cell array.

%                Gate                  |  Column 1  |   Column 2
% -----------------------------------------------------------------
% Pauli X on qubits 2,4                |     'X'    |   [2 4]
% Pauli Z on qubits 1,5                |     'Z'    |   [1 5]
% Pauli Y on qubits 1,2,5              |     'Y'    |   [1 2 5]
% Phase on qubits 1,3                  |     'S'    |   [1 3]
% Hadamard on qubits 2,4,5             |     'H'    |   [2 4 5]
% Controlled-Z on qubits 3,6           |    'CZ'    |   [3 6]
% Controlled-NOT: qubit 2 controls 1   |   'CNOT'   |   [2 1]
% Permutation (m=3): [1 2 3] -> [2 3 1]|  'Permute' |   [2 3 1]
% -----------------------------------------------------------------

% Example Circuit (m = 6 qubits): 
% Unitary U = CZ_{26} * H1 * CNOT_{12} * H2 * CNOT_{24} * H3 * CZ_{14}
% In circuit diagram the last CZ_{14} will appear first. This is because the
% operator acts on state |v> as U|v>, and so |v> goes through the last
% CZ_{14} first. Hence this is the required order for this function too.

% In this case, our specification for this function will be:
% circuit = {'CZ', [1 4]; 'H', 3; 'CNOT', [2 4]; 'H', 2; 'CNOT', [1 2]; ...
%                         'H', 1; 'CZ', [2 6]};
% If we want to compute its effect on the Pauli operator X_1*Y_2*Z_4*X_6,
% then op_in = {'XYZX', [1, 2, 4, 6]}.

% Author: Narayanan Rengaswamy, Date: Mar. 2, 2018

I = eye(2);
X = [0 1; 1 0];
Z = [1 0; 0 -1];
Y = sqrt(-1) * X * Z;
S = [1 0; 0 sqrt(-1)];
H = 1/sqrt(2) * (X + Z);
e0 = [1; 0];
e1 = [0; 1];
E00 = e0 * e0';
E11 = e1 * e1';
CNOT = kron(E00, I) + kron(E11, X);
CZ = kron(E00, I) + kron(E11, Z);

Pauli_in = cell(m,1);
for i = 1:m
    ind = find(op_in{1,2} == i, 1);
    if (~isempty(ind))
        P = op_in{1,1}(ind);
    else
        P = 'I';
    end
    if (strcmpi(P, 'X'))
        Pauli_in{i,1} = X;
    elseif (strcmpi(P, 'Z'))
        Pauli_in{i,1} = Z;
    elseif (strcmpi(P, 'Y'))
        Pauli_in{i,1} = Y;
    else
        Pauli_in{i,1} = I;
    end
end

Pauli_out = Pauli_in;
for i = 1:size(circuit,1)
    gate = circuit{i,1};
    qubits = circuit{i,2}(:)';
    if (strcmpi(gate, 'X'))
        if (isempty(qubits))
            fprintf('\nPauli X Gate: Need to specify atleast one qubit!\n');
            op_out = [];
            return;
        end
        for j = 1:length(qubits)
            Pauli_out{qubits(j),1} = X * Pauli_out{qubits(j),1} * X';
        end
    elseif (strcmpi(gate, 'Z'))
        if (isempty(qubits))
            fprintf('\nPauli Z Gate: Need to specify atleast one qubit!\n');
            op_out = [];
            return;
        end
        for j = 1:length(qubits)
            Pauli_out{qubits(j),1} = Z * Pauli_out{qubits(j),1} * Z';
        end
    elseif (strcmpi(gate, 'Y'))
        if (isempty(qubits))
            fprintf('\nPauli Y Gate: Need to specify atleast one qubit!\n');
            op_out = [];
            return;
        end
        for j = 1:length(qubits)
            Pauli_out{qubits(j),1} = Y * Pauli_out{qubits(j),1} * Y';
        end
    elseif (strcmpi(gate, 'S'))
        if (isempty(qubits))
            fprintf('\nPhase Gate: Need to specify atleast one qubit!\n');
            op_out = [];
            return;
        end
        for j = 1:length(qubits)
            Pauli_out{qubits(j),1} = S * Pauli_out{qubits(j),1} * S';
        end
    elseif (strcmpi(gate, 'H'))
        if (isempty(qubits))
            fprintf('\nHadamard Gate: Need to specify atleast one qubit!\n');
            op_out = [];
            return;
        end
        for j = 1:length(qubits)
            Pauli_out{qubits(j),1} = H * Pauli_out{qubits(j),1} * H';
        end
    elseif (strcmpi(gate, 'CNOT'))
        if (isempty(qubits) || length(qubits) ~= 2)
            fprintf('\nCNOT Gate: Need to specify two qubits!\n');
            op_out = [];
            return;
        end
        out12 = find_other(CNOT, Pauli_out{qubits(1),1}, 1);
        out21 = find_other(CNOT, Pauli_out{qubits(2),1}, 2);
        Pauli_out{qubits(1),1} = out21 * Pauli_out{qubits(1),1};
        Pauli_out{qubits(2),1} = Pauli_out{qubits(2),1} * out12;
    elseif (strcmpi(gate, 'CZ'))
        if (isempty(qubits) || length(qubits) ~= 2)
            fprintf('\nCZ Gate: Need to specify two qubits!\n');
            op_out = [];
            return;
        end
        out12 = find_other(CZ, Pauli_out{qubits(1),1}, 1);
        out21 = find_other(CZ, Pauli_out{qubits(2),1}, 2);
        Pauli_out{qubits(1),1} = out21 * Pauli_out{qubits(1),1};
        Pauli_out{qubits(2),1} = Pauli_out{qubits(2),1} * out12;
    elseif (strcmpi(gate, 'Permute'))
        desired_order = circuit{i,2}(:)';
        if (isempty(desired_order) || length(desired_order) ~= m)
            fprintf('\nPermutation: Need to specify %d qubits!\n', m);
            op_out = [];
            return;
        end
        Pauli_out(1:m, 1) = Pauli_out(desired_order, 1);
    end
end

out_sign = 1;
op_out = {'', []};
for i = 1:m
    Pauli_out{i,1} = round(Pauli_out{i,1});
    if (all(all(Pauli_out{i,1} == X)))
        op_out{1,1} = strcat(op_out{1,1}, 'X');
        op_out{1,2} = [op_out{1,2}, i];
    elseif (all(all(Pauli_out{i,1} == -X)))
        op_out{1,1} = strcat(op_out{1,1}, 'X');
        op_out{1,2} = [op_out{1,2}, i];
        out_sign = out_sign * (-1);
    elseif (all(all(Pauli_out{i,1} == sqrt(-1)*X)))
        op_out{1,1} = strcat(op_out{1,1}, 'X');
        op_out{1,2} = [op_out{1,2}, i];
        out_sign = out_sign * sqrt(-1);
    elseif (all(all(Pauli_out{i,1} == -sqrt(-1)*X)))
        op_out{1,1} = strcat(op_out{1,1}, 'X');
        op_out{1,2} = [op_out{1,2}, i];
        out_sign = out_sign * (-1) * sqrt(-1);
    elseif (all(all(Pauli_out{i,1} == Z)))
        op_out{1,1} = strcat(op_out{1,1}, 'Z');
        op_out{1,2} = [op_out{1,2}, i];
    elseif (all(all(Pauli_out{i,1} == -Z)))
        op_out{1,1} = strcat(op_out{1,1}, 'Z');
        op_out{1,2} = [op_out{1,2}, i];
        out_sign = out_sign * (-1);
    elseif (all(all(Pauli_out{i,1} == sqrt(-1)*Z)))
        op_out{1,1} = strcat(op_out{1,1}, 'Z');
        op_out{1,2} = [op_out{1,2}, i];
        out_sign = out_sign * sqrt(-1);
    elseif (all(all(Pauli_out{i,1} == -sqrt(-1)*Z)))
        op_out{1,1} = strcat(op_out{1,1}, 'Z');
        op_out{1,2} = [op_out{1,2}, i];
        out_sign = out_sign * (-1) * sqrt(-1);
    elseif (all(all(Pauli_out{i,1} == Y)))
        op_out{1,1} = strcat(op_out{1,1}, 'Y');
        op_out{1,2} = [op_out{1,2}, i];
    elseif (all(all(Pauli_out{i,1} == -Y)))
        op_out{1,1} = strcat(op_out{1,1}, 'Y');
        op_out{1,2} = [op_out{1,2}, i];
        out_sign = out_sign * (-1);
    elseif (all(all(Pauli_out{i,1} == sqrt(-1)*Y)))
        op_out{1,1} = strcat(op_out{1,1}, 'Y');
        op_out{1,2} = [op_out{1,2}, i];
        out_sign = out_sign * sqrt(-1);
    elseif (all(all(Pauli_out{i,1} == -sqrt(-1)*Y)))
        op_out{1,1} = strcat(op_out{1,1}, 'Y');
        op_out{1,2} = [op_out{1,2}, i];
        out_sign = out_sign * (-1) * sqrt(-1);
    elseif (all(all(Pauli_out{i,1} == -I)))
        out_sign = out_sign * (-1);
    elseif (all(all(Pauli_out{i,1} == sqrt(-1)*I)))
        out_sign = out_sign * sqrt(-1);
    elseif (all(all(Pauli_out{i,1} == -sqrt(-1)*I)))
        out_sign = out_sign * (-1) * sqrt(-1);
    elseif (all(all(Pauli_out{i,1} == I)))
        continue;
    else
        fprintf('\ncalculate_conj: Unknown gate encountered...\n');
        op_out = [];
    end
end
if (out_sign == -1)
    op_out{1,1} = strcat('-', op_out{1,1});
    op_out{1,2} = [-1, op_out{1,2}];
elseif (out_sign == sqrt(-1))
    op_out{1,1} = strcat('i', op_out{1,1});
    op_out{1,2} = [sqrt(-1), op_out{1,2}];
elseif (out_sign == -sqrt(-1))
    op_out{1,1} = strcat('j', op_out{1,1});  % j represents -i
    op_out{1,2} = [-sqrt(-1), op_out{1,2}];    
end
    

    function other = find_other(gate, inp, id)
        % 'gate' is a 2-qubit Clifford gate - a 4 x 4 matrix
        % 'inp' is a 1-qubit Pauli gate - a 2 x 2 matrix
        % 'id' is the terminal of the gate where 'inp' is input - 1 or 2
        if (id == 1)
            Pin = kron(inp, I);
            PinX = kron(inp, X);
            PinZ = kron(inp, Z);
            PinY = kron(inp, Y);
        else
            Pin = kron(I, inp);
            PinX = kron(X, inp);
            PinZ = kron(Z, inp);
            PinY = kron(Y, inp);
        end
        Pout = gate * Pin * gate';
        if (all(all(Pout == PinX)))
            other = X;
        elseif (all(all(Pout == PinZ)))
            other = Z;
        elseif (all(all(Pout == PinY)))
            other = Y;
        elseif (all(all(Pout == Pin)))
            other = I;
        else
            fprintf('\ncalculate_conj: Something wrong in find_other...');
            other = [];
        end
    end

end