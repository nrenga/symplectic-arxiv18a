function [F, Transvecs] = find_symp_mat_transvecs(X, Y)
% Find a binary symplectic matrix F that satisfies X*F = Y

% Author: Narayanan Rengaswamy, Date: Mar. 1, 2018

m = size(X,1);
n = size(X,2)/2;
F = eye(2*n);
Z_h = @(h) ( mod(eye(2*n) + mod(fftshift(h)' * h, 2), 2) );  % h is a row vector
Transvecs = cell(1,2);
ind = 1;

for i = 1:m
    x_i = X(i,:);
    y_i = Y(i,:);
    x_it = mod(x_i * F, 2);
    if (all(x_it == y_i))
        continue;
    end
    if (symp_inn_pdt(x_it, y_i) == 1)
        h_i = mod(x_it + y_i, 2);
        F = mod(F * Z_h(h_i), 2);
        Transvecs{ind,1} = Z_h(h_i);
        Transvecs{ind,2} = h_i;
        ind = ind + 1;
    else
        % Need to pick w_i s.t. <x_it, w_i> = <w_i, y_i> = 1 and
        %                       <y_j, w_i> = <y_j, y_i> for all j < i
        w_i = find_w(x_it, y_i, Y(1:i-1,:));
        h_i1 = mod(w_i + y_i, 2);
        h_i2 = mod(x_it + w_i, 2);
        F = mod(mod(F * Z_h(h_i1), 2) * Z_h(h_i2), 2);
        Transvecs{ind,1} = Z_h(h_i1);
        Transvecs{ind,2} = h_i1;
        ind = ind + 1;
        Transvecs{ind,1} = Z_h(h_i2);
        Transvecs{ind,2} = h_i2;
        ind = ind + 1;
    end
end

    function w = find_w(x,y,Ys)
        A = fftshift([x; y; Ys],2);
        b = [1; 1; symp_inn_pdt(Ys,repmat(y,size(Ys,1),1))];
        w = gflineq(A,b)';
    end

end