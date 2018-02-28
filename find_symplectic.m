function [F, A, B, C, D] = find_symplectic(n, circuit)
% 'n' is the number of qubits in the system
% 'circuit' is an 'm x 2' cell array with each row defining a gate as 
% specified below; 'm' is the number of gates (stages) in the circuit
% 'F' is the overall '2n x 2n' symplectic matrix for the circuit, which is
% of the form F = [ A, B; C, D]

% Each row of the cell array must be one of the following Clifford gates:
% Gate is the desired gate; columns 1 and 2 are specifications for the gate
% that form a row of the 'circuit' cell array.
%                Gate                  |  Column 1  |   Column 2
% -----------------------------------------------------------------
% Phase on qubit 1                     |     'S'    |   [1]
% Phase on qubits 1,3                  |     'S'    |   [1 3]
% Hadamard on qubit 4                  |     'H'    |   [4]
% Hadamard on qubits 2,4,5             |     'H'    |   [2 4 5]
% Controlled-Z on qubits 3,6           |    'CZ'    |   [3 6]
% Controlled-NOT: qubit 2 controls 1   |   'CNOT'   |   [2 1]
% -----------------------------------------------------------------
% Note: The symplectic matrix for all Pauli gates is the identity matrix.
%       So you could ignore those gates in the 'circuit' cell array.

% Example Circuit (n = 6): 
% U = CZ_{26} * H1 * CNOT_{12} * H2 * CNOT_{24} * H3 * CZ_{14}
% In circuit diagram the last CZ_{14} will appear first. This is because the
% operator acts on state |v> as U|v>, and so |v> goes through the last
% CZ_{14} first. Hence this is the required order for this function too.

% In this case, our specification for this function will be:
% circuit = {'CZ', [1 4]; 'H', 3; 'CNOT', [2 4]; 'H', 2; 'CNOT', [1 2]; ...
%                         'H', 1; 'CZ', [2 6]};


F = eye(2*n);
m = size(circuit,1);
I = eye(n);
O = zeros(n);

for i = 1:m
    gate = circuit{i,1};
    if (strcmpi(gate, 'S'))
        qubits = circuit{i,2}(:)';
        if (isempty(qubits))
            continue;
        end
        B = diag(mod(sum(I(qubits,:),1),2));
        FS = [I, B; O, I];
        F = mod(F * FS, 2);
    elseif (strcmpi(gate, 'CZ'))
        qubits = circuit{i,2}(:)';
        if (isempty(qubits) || length(qubits) ~= 2)
            continue;
        end
        B = zeros(n);
        B(qubits(1), qubits(2)) = 1; 
        B(qubits(2), qubits(1)) = 1;
        FCZ = [I, B; O, I];
        F = mod(F * FCZ, 2);
    elseif (strcmpi(gate, 'CNOT'))
        qubits = circuit{i,2}(:)';
        if (isempty(qubits) || length(qubits) ~= 2)
            continue;
        end
        M = I;
        M(qubits(1), qubits(2)) = 1; 
        FCNOT = blkdiag(M, mod(inv(M),2)');
        F = mod(F * FCNOT, 2);
    elseif (strcmpi(gate, 'H'))
        qubits = circuit{i,2}(:)';
        if (isempty(qubits))
            continue;
        end
        FH = eye(2*n);
        FH([qubits, n + qubits],:) = fftshift(FH([qubits, n + qubits],:), 2);
        F = mod(F * FH, 2);
    else  % handles all Pauli gates and unrecognized gates
        continue;
    end
end
        
A = F(1:n, 1:n);
B = F(1:n, n+(1:n));
C = F(n+(1:n), 1:n);
D = F(n+(1:n), n+(1:n));

end