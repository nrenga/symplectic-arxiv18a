function PG_logical = pauli_group_logical(S, Xbar, Zbar, PG)
% Function to return m qubit Pauli elements as gates and unitaries realizing 
% logical Paulis on the (m-k) protected qubits of a given [[m,m-k]] 
% stabilizer code.

% The inputs 'S', 'Xbar', 'Zbar' are k x 2m, (m-k) x 2m, and (m-k) x 2m
% binary matrices, respectively, specifying the stabilizer generators,
% logical Pauli X generators, and logical Pauli Z generators in their
% binary vector representation.

% The input 'PG' is a binary matrix whose rows represent the Paulis on
% (m-k) qubits for which physical realizations are desired. Note that this
% matrix needs to have 2*(m-k) columns. If this input is not given, then by
% default PG will be taken to be all 2^(2m) Pauli elements on (m-k) qubits.

% 'PG_logical' is a cell array with three columns; the first column contains
% circuit representation of the Pauli element, the second column
% contains the unitary matrix representation, and the third column contains
% the binary vector representation.

% Author: Narayanan Rengaswamy, Date: Sep. 15, 2018

[k, m] = size(S);
m = m/2;
XZbar = [Xbar; Zbar];

if (nargin == 3)
    PG = pauli_group(m-k);
    PG = cell2mat(PG(:,3));    % Just need the binary representation of the Paulis
end

PG_logical = cell(size(PG,1),3);
for i = 1:size(PG,1)
    bi = mod(PG(i,:) * XZbar, 2);
    bi2 = bi(1:m) + sqrt(-1) * bi(m+(1:m));
    Xinds = find(bi2 == 1);
    Yinds = find(bi2 == 1 + sqrt(-1));
    Zinds = find(bi2 == sqrt(-1));
    gates = [];
    if (~isempty(Xinds))
        gates = [gates; {'X', Xinds}];
    end
    if (~isempty(Yinds))
        gates = [gates; {'Y', Yinds}];
    end
    if (~isempty(Zinds))
        gates = [gates; {'Z', Zinds}];
    end
    PG_logical{i,1} = gates;
    PG_logical{i,2} = find_unitary(m,gates);
    PG_logical{i,3} = bi;
end

end