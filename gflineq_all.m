function X = gflineq_all(a, b)

[m_a, n_a] = size(a);
% [m_b, n_b] = size(b);

% The following Gaussian elimination code is borrowed from gflineq.m
% -----------------------------------------------------------------------

% Construct an AA = [A, B] composite matrix and assign initial values.
aa = [a b];
[m_aa, n_aa] = size(aa);
p = 2;

row_idx = 1;
column_idx = 1;
row_store = [];
column_store = [];

% Find the multiplicative inverse of the field elements.
% This will be used for setting major elements in the matrix to one.
[field_inv ignored] = find( mod( (1:(p-1)).' * (1:(p-1)) , p ) == 1 );  %#ok

% Search for major elements, trying to form 'identity' matrix.
while (row_idx <= m_aa) && (column_idx < n_aa)
    
    % Look for a major element in the current column.
    while (aa(row_idx,column_idx) == 0) && (column_idx < n_aa)
        
        % In the current column, search below all the rows that already
        % have major elements.
        idx = find( aa(row_idx:m_aa, column_idx) ~= 0 );
        
        if isempty(idx)
            % There are no major elements in this column.
            % Move to the next column.
            column_idx = column_idx + 1;
            
        else
            % There is a major element in this column.
            % See if any are already equal to one.
            idx = [ find(aa(row_idx:m_aa, column_idx) == 1); idx ];
            idx = idx(1);
            
            % Swap the current row with a row containing a major element.
            temp_row = aa(row_idx,:);
            aa(row_idx,:) = aa(row_idx+idx-1,:);
            aa(row_idx+idx-1,:) = temp_row;
            
        end
    end
    
    
    % Clear all non-zero elements in the column except the major element,
    % and set the major element to one.
    if ( ( aa(row_idx,column_idx) ~= 0 ) && ( column_idx < n_aa ) )
        
        % The current element is a major element.
        row_store = [row_store row_idx];
        column_store = [column_store column_idx];
        
        % If the major element is not already one, set it to one.
        if (aa(row_idx,column_idx) ~= 1)
            aa(row_idx,:) = mod( field_inv( aa(row_idx,column_idx) ) * aa(row_idx,:), p );
        end
        
        % Find the other elements in the column that must be cleared,
        idx = find(aa(:,column_idx)~=0)';
        % and set those elements to zero.
        for i = idx
            if i ~= row_idx
                aa(i,:) = mod( aa(i,:) + aa(row_idx,:) * (p - aa(i,column_idx)), p );
            end
        end
        
        column_idx = column_idx + 1;
        
    end
    
    row_idx = row_idx + 1;
    
end
% -----------------------------------------------------------------------

non_pivot = setdiff(1:n_a, column_store);
n_np = length(non_pivot);
X = zeros(n_a,1);
no_of_Xs = 0;

for i = 0:(2^n_np-1)
    x = zeros(n_a, 1);
    x(non_pivot) = de2bi(i, n_np, 'left-msb')';
    bnew = mod(aa(:,1:n_a) * x + aa(:,n_aa), 2);
    x(column_store) = bnew(1:length(column_store));
    no_of_Xs = no_of_Xs + 1;
    X(:,no_of_Xs) = x;
end

end
