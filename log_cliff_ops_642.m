%% Solve symplectic linear system to get logical Clifford operators 
% Example: The [[6,4,2]] CSS code. 
% This script is standalone unless the commented lines suggesting to use
% the functions qfind_all_symp_mat or find_all_symp_mat are uncommented.

% This script reproduces the results published in Appendix II of the
% paper "Synthesis of Logical Clifford Operators via Symplectic Geometry",
% available at https://arxiv.org/abs/

% Author: Narayanan Rengaswamy, Date: Feb. 20, 2018

clc
clear
close all

% [[m, m-k, d]] = [[6, 4, 2]]
m = 6;
k = 2;

% Stabilizers
S = [ 1 1 1 1 1 1, 0 0 0 0 0 0 ;
      0 0 0 0 0 0, 1 1 1 1 1 1 ];

% Logical Paulis
Xbar = [ 1 1 0 0 0 0, 0 0 0 0 0 0 ;
         1 0 1 0 0 0, 0 0 0 0 0 0 ;
         1 0 0 1 0 0, 0 0 0 0 0 0 ;
         1 0 0 0 1 0, 0 0 0 0 0 0 ];

Zbar = [ 0 0 0 0 0 0, 0 1 0 0 0 1 ;
         0 0 0 0 0 0, 0 0 1 0 0 1 ;
         0 0 0 0 0 0, 0 0 0 1 0 1 ;
         0 0 0 0 0 0, 0 0 0 0 1 1 ];

Ybar = mod(Xbar + Zbar, 2);     

% A symplectic basis for \mathbb{F}_2^{12} using S, Xbar, Zbar
U = [ Xbar; S; Zbar; [zeros(1,11) 1]; [1 zeros(1,11)]];

%% Solve for logical Phase gate S1bar

G_S1 = [ Xbar(1,:); Xbar(2:end,:); S; Zbar];
H_S1 = [ Ybar(1,:); Xbar(2:end,:); S; Zbar];

% Find one solution using transvections
F_S1 = find_one_symp(G_S1,H_S1);

% Find all solutions using dual subspace idea
F_all_S1 = qfind_all_symp(k, U, F_S1);    

% Can also do one of the following:
% F_all_S1 = qfind_all_symp_mat(U, H_S1);    
% F_all_S1 = find_all_symp_mat(U, H_S1, 1:m, 1:m-k);    

%% Solve for logical CZ gate CZ12bar

G_CZ12 = [ Xbar(1,:); Xbar(2,:); Xbar(3:end,:); S; Zbar];
H_CZ12 = [ mod(Xbar(1,:) + Zbar(2,:), 2); ...
           mod(Xbar(2,:) + Zbar(1,:), 2); Xbar(3:end,:); S; Zbar];

% Find one solution using transvections
F_CZ12 = find_one_symp(G_CZ12,H_CZ12);

% Find all solutions using dual subspace idea
F_all_CZ12 = qfind_all_symp(k, U, F_CZ12);

% Can also do one of the following:
% F_all_CZ12 = qfind_all_symp_mat(U, H_CZ12);    
% F_all_CZ12 = find_all_symp_mat(U, H_CZ12, 1:m, 1:m-k);    

%% Solve for logical CNOT gate CNOT21bar

G_CNOT21 = [ Xbar(1,:); Xbar(2,:); Xbar([3 4],:); S; Zbar(1,:); Zbar(2:end,:)];
H_CNOT21 = [ Xbar(1,:); mod(Xbar(1,:) + Xbar(2,:), 2); Xbar([3 4],:); S; ...
             mod(Zbar(1,:) + Zbar(2,:), 2); Zbar(2:end,:)];

% Find one solution using transvections
F_CNOT21 = find_one_symp(G_CNOT21,H_CNOT21);

% Find all solutions using dual subspace idea
F_all_CNOT21 = qfind_all_symp(k, U, F_CNOT21);

% Can also do one of the following:
% F_all_CNOT21 = qfind_all_symp_mat(U, H_CNOT21);    
% F_all_CNOT21 = find_all_symp_mat(U, H_CNOT21, 1:m, 1:m-k);    

%% Solve for logical targered Hadamard H1bar

G_H1 = [ Xbar(1,:); Xbar(2:end,:); S; Zbar(1,:); Zbar(2:end,:)];
H_H1 = [ Zbar(1,:); Xbar(2:end,:); S; Xbar(1,:); Zbar(2:end,:)];

% Find one solution using transvections
F_H1 = find_one_symp(G_H1,H_H1);

% Find all solutions using dual subspace idea
F_all_H1 = qfind_all_symp(k, U, F_H1);

% Can also do one of the following:
% F_all_H1 = qfind_all_symp_mat(U, H_H1);    
% F_all_H1 = find_all_symp_mat(U, H_H1, 1:m, 1:m-k);    

%% Function to calculate one solution for each logical Clifford operator

function F = find_one_symp(X, Y)
% Find a binary symplectic matrix F that satisfies X*F = Y

m = size(X,1);
n = size(X,2)/2;
F = eye(2*n);
Z_h = @(h) ( mod(eye(2*n) + mod(fftshift(h)' * h, 2), 2) );  % h is a row vector

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
    else
        % Need to pick w_i s.t. <x_it, w_i> = <w_i, y_i> = 1 and
        %                       <y_j, w_i> = <y_j, y_i> for all j < i
        w_i = find_w(x_it, y_i, Y(1:i-1,:));
        h_i1 = mod(w_i + y_i, 2);
        h_i2 = mod(x_it + w_i, 2);
        F = mod(mod(F * Z_h(h_i1), 2) * Z_h(h_i2), 2);
    end
end

    function w = find_w(x,y,Ys)
        A = fftshift([x; y; Ys],2);
        b = [1; 1; symp_inn_pdt(Ys,repmat(y,size(Ys,1),1))];
        w = gflineq(A,b)';
    end

end

%% Function to calculate all solutions for each logical Clifford operator

function F_all = qfind_all_symp(k, U, F0)
% U = [Xbar; S; Zbar; Spair]
% Spair completes the symplectic basis for \mathbb{F}_2^{2m}

m = size(U,2)/2;
tot = 2^(k*(k+1)/2);
F_all = cell(tot,1);

A = mod(U * F0, 2);
Ainv = g2matinv(A);
Basis = A([(m - k + 1):m, (2*m - k + 1):(2*m)],:);  % Basis = [S; Spair]
Subspace = mod(de2bi((0:2^(2*k)-1)',2*k) * Basis, 2);
StabF0 = A((m - k + 1):m, :);  % StabF0 = S*F0
Choices = cell(k,1);

% Calculate all choices for each vector in Spair using just the conditions
% imposed by the (modified) stabilizer generators, i.e., the rows of S*F0
for i = 1:k
    % Impose symplectic inner product of 1 with the "fixed" stabilizer
    h = [zeros(1,i-1) 1 zeros(1,k-i)];
    Innpdts = mod(Subspace * fftshift(StabF0, 2)', 2);
    Choices{i,1} = Subspace(bi2de(Innpdts) == bi2de(h), :);
end

% First free vector, i.e. row 1 of Spair has 2^k choices, 
% second free vector has 2^(k-1) choices and so on...
for l = 0:(tot - 1)
    Bl = A;
    V = zeros(k,2*m);
    lbin = de2bi(l,k*(k+1)/2,'left-msb');
    v1_ind = bi2de(lbin(1,1:k),'left-msb') + 1;
    V(1,:) = Choices{1,1}(v1_ind,:);
    for i = 2:k
        vi_ind = bi2de(lbin(1,sum(k:-1:k-(i-2)) + (1:(k-(i-1)))),'left-msb') + 1;
        Innprods = mod(Choices{i,1} * fftshift(V,2)', 2);
        
        % Impose symplectic inner product of 0 with chosen free vectors
        Ch_i = Choices{i,1}(bi2de(Innprods) == 0, :);
        V(i,:) = Ch_i(vi_ind,:);
    end
    Bl((2*m - k + 1):(2*m), :) = V;
    F = mod(Ainv * Bl, 2);
    F_all{l+1,1} = mod(F0 * F, 2);
end    

end
