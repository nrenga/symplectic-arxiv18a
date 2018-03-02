function [F, A, B, C, D] = find_symplectic(m, circuit)
% 'm' is the number of qubits in the system
% 'circuit' is a 'd x 2' cell array with each row defining a gate as 
% specified below; 'd' is the depth of the circuit
% 'F' is the overall '2m x 2m' symplectic matrix for the circuit, which is
% of the form F = [ A, B; C, D]. It satisfies F*Omega*F' = Omega, where
% Omega = [0 I_m; I_m 0].

% Each row of the cell array must be one of the following Clifford gates:
% Gate is the desired gate; columns 1 and 2 are specifications for the gate
% that form a row of the 'circuit' cell array.

%                Gate                  |  Column 1  |   Column 2
% -----------------------------------------------------------------
% Phase on qubit 1                     |     'S'    |   [1]
% Phase on qubits 1,3                  |     'S'    |   [1 3]
% Hadamard on qubits 2,4,5             |     'H'    |   [2 4 5]
% Controlled-Z on qubits 3,6           |    'CZ'    |   [3 6]
% Controlled-NOT: qubit 2 controls 1   |   'CNOT'   |   [2 1]
% Permutation: (m=3) [1 2 3] -> [2 3 1]|  'Permute' |   [2 3 1]
% -----------------------------------------------------------------
% Note: The symplectic matrix for all Pauli gates is the identity matrix.
%       So you could ignore those gates in the 'circuit' cell array.

% Example Circuit (m = 6): 
% U = CZ_{26} * H1 * CNOT_{12} * H2 * CNOT_{24} * H3 * CZ_{14}
% In circuit diagram the last CZ_{14} will appear first. This is because the
% operator acts on state |v> as U|v>, and so |v> goes through the last
% CZ_{14} first. Hence this is the required order for this function too.

% In this case, our specification for this function will be:
% circuit = {'CZ', [1 4]; 'H', 3; 'CNOT', [2 4]; 'H', 2; 'CNOT', [1 2]; ...
%                         'H', 1; 'CZ', [2 6]};

% Author: Narayanan Rengaswamy, Date: Mar. 1, 2018


F = eye(2*m);
d = size(circuit,1);
I = eye(m);
O = zeros(m);

for i = 1:d
    gate = circuit{i,1};
    qubits = circuit{i,2}(:)';
    if (strcmpi(gate, 'S'))
        if (isempty(qubits))
            fprintf('\nPhase Gate: Need to specify atleast one qubit!\n');
            F = [];
            return;
        end
        B = diag(mod(sum(I(qubits,:),1),2));
        FS = [I, B; O, I];
        F = mod(F * FS, 2);
    elseif (strcmpi(gate, 'H'))
        if (isempty(qubits))
            fprintf('\nHadamard Gate: Need to specify atleast one qubit!\n');
            F = [];
            return;
        end
        FH = eye(2*m);
        FH([qubits, m + qubits],:) = fftshift(FH([qubits, m + qubits],:), 2);
        F = mod(F * FH, 2);
    elseif (strcmpi(gate, 'CZ'))
        if (isempty(qubits) || length(qubits) ~= 2)
            fprintf('\nCZ Gate: Need to specify two qubits!\n');
            F = [];
            return;
        end
        B = zeros(m);
        B(qubits(1), qubits(2)) = 1; 
        B(qubits(2), qubits(1)) = 1;
        FCZ = [I, B; O, I];
        F = mod(F * FCZ, 2);
    elseif (strcmpi(gate, 'CNOT'))
        if (isempty(qubits) || length(qubits) ~= 2)
            fprintf('\nCNOT Gate: Need to specify two qubits!\n');
            F = [];
            return;
        end
        M = I;
        M(qubits(1), qubits(2)) = 1; 
        FCNOT = blkdiag(M, gf2matinv(M)');
        F = mod(F * FCNOT, 2);
    elseif (strcmpi(gate, 'Permute'))
        if (isempty(qubits) || length(qubits) ~= m)
            fprintf('\nPermutation: Need to specify %d qubits!\n', m);
            F = [];
            return;
        end
        M = I;
        M = M(:,qubits); 
        FPermute = blkdiag(M, gf2matinv(M)');
        F = mod(F * FPermute, 2);
    else  % handles all Pauli gates and unrecognized gates
        fprintf('\nfind_symplectic: Unrecognized gate encountered!\n');
        F = [];
        return;
    end
end
        
A = F(1:m, 1:m);
B = F(1:m, m+(1:m));
C = F(m+(1:m), 1:m);
D = F(m+(1:m), m+(1:m));

end