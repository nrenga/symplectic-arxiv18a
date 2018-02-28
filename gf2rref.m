% This is a modified version of matlab's building rref which calculates
% row-reduced echelon form in gf(2).  Useful for linear codes.
% Tolerance was removed because yolo, and because all values
% should only be 0 or 1.  @benathon

function [Arref, M, N, rnk] = gf2rref(A)
%G2RREF   Reduced row echelon form in gf(2).
%   R = RREF(A) produces the reduced row echelon form of A in gf(2).
%
%   Class support for input A:
%      float: with values 0 or 1
%   Copyright 1984-2005 The MathWorks, Inc. 
%   $Revision: 5.9.4.3 $  $Date: 2006/01/18 21:58:54 $

% Modified to return the matrix M of row operations on A, i.e., Arref = M*A,
% and the matrix N of column operations which if applied to Arref results in
% a matrix of the form [I_rnk, 0; 0, 0] for the first m columns, 
% where rnk is the gf(2) rank of A, 
% i.e., (Arref*N)_{1:m,1:m} = (M*A*N)_{1:m,1:m} = [I_rnk, 0; 0, 0].
% For a square matrix A, Arref*N = M*A*N = [I_rnk, 0; 0, 0].

% By Narayanan Rengaswamy. Date: Feb. 28, 2018

[Arref, M] = gf2redref(A);
[Ardiag, Nt] = gf2redref(Arref');
N = Nt';
rnk = sum(diag(Ardiag));

    function [Arref, Row_ops] = gf2redref(A)
        [m,n] = size(A);
        Row_ops = eye(m);
        Arref = [A, Row_ops];
        nr = size(Arref, 2);
        
        % Loop over the entire matrix.
        i = 1;
        j = 1;
        
        while (i <= m) && (j <= n)
            while (Arref(i,j) == 0) && (j <= n)
                % Find value and index of largest element in the remainder of column j.
                k = find(Arref(i:m,j),1) + i - 1;
                
                if (isempty(k))
                    j = j + 1;
                    continue;
                else
                    % Swap i-th and k-th rows.
                    Arref([i k],j:nr) = Arref([k i],j:nr);
                end
            end
            
            if (Arref(i,j) == 1) && (j <= n)
                % Save the right hand side of the pivot row
                aijn = Arref(i,j:nr);
                
                % Column we're looking at
                col = Arref(1:m,j);
                
                % Never Xor the pivot row against itself
                col(i) = 0;
                
                % This builds an matrix of bits to flip
                flip = col*aijn;
                
                % Xor the right hand side of the pivot row with all the other rows
                Arref(1:m,j:nr) = xor( Arref(1:m,j:nr), flip );
                
                j = j + 1;
            end
            
            i = i + 1;
        end
        Row_ops = Arref(1:m,(n+1):nr);
        Arref = Arref(1:m,1:n);
    end

end