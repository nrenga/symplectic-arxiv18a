function PG = pauli_group(m, scalars)
% Function to return Pauli elements as gates and unitaries on 'm' qubits

% The parameter 'scalars' is a Boolean variable indicating whether the
% returning set 'PG' needs to be a set (if scalars = 0) or a group (if
% scalars = 1). If scalars = 0, then 'PG' contains all phaseless
% elements of the Pauli group on m qubits.

% 'PG' is a cell array with three columns; the first column contains
% circuit representation of the Pauli element, the second column
% contains the unitary matrix representation, and the third column contains
% the binary vector representation.

% Author: Narayanan Rengaswamy, Date: Sep. 15, 2018

if (nargin == 1)
    scalars = 0;
end

if (scalars)
    f = 4;
else
    f = 1;
end

PG = cell(f*2^(2*m),3);
for i = 0:(2^(2*m)-1)
    bi = de2bi(i,2*m);
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
    PG{f*i+1,1} = gates;
    PG{f*i+1,2} = find_unitary(m,gates);
    PG{f*i+1,3} = bi;
    if (scalars)
        PG(f*i+2,1:3) = {gates, 1i*find_unitary(m,gates), bi};
        PG(f*i+3,1:3) = {gates, (-1)*find_unitary(m,gates), bi};
        PG(f*i+4,1:3) = {gates, (-1i)*find_unitary(m,gates), bi};
    end
end

end