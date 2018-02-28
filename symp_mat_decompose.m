function Decomp = symp_mat_decompose(F)
% Function to decompose an arbitrary 2m x 2m binary symplectic matrix F
% into a product of elementary symplectic transformations from below.
% Based on Trung Can's algorithm

% Author: Narayanan Rengaswamy, Date: Feb. 28, 2018

m = size(F,1)/2;
I = eye(m);
O = zeros(m);
U = @(k) (blkdiag(eye(k), zeros(m-k)));
L = @(k) (blkdiag(zeros(m-k), eye(k)));

% Elementary Symplectic Transformations
Omega = [O, I; I, O];     % Transversal Hadamards
Elem1 = @(Q) (blkdiag(Q, gf2matinv(Q)'));  % CNOTs, Permutations
Elem2 = @(P) ([I, P; O, I]);    % CZs, Phase gates
Elem3 = @(k) ([L(m-k), U(k); U(k), L(m-k)]);   % Partial Hadamards

% Step 1
A = F(1:m, 1:m);
[~, M_A, N_A, k] = gf2rref(A);
Qleft1 = Elem1(M_A);
Qright = Elem1(N_A);
Fcp = mod(Qleft1 * F * Qright, 2);
% if (k == m)
%     Pright = Elem2(Fcp(1:m, m+(1:m)));
%     Fcp = mod(Fcp * Pright, 2);  % Will result in [I_m, Z; Cp, I_m]
%     P = Fcp(m+(1:m), 1:m);
%     Decomp = {gf2matinv(Qleft1); Omega; Elem2(P); Omega; gf2matinv(Pright); gf2matinv(Qright)};
%     return;
% end

% Step 2
Bmk = Fcp((k+1):m, m+((k+1):m));
[~, M_Bmk1, ~, ~] = gf2rref(Bmk);
M_Bmk = blkdiag(eye(k), M_Bmk1);
Qleft2 = Elem1(M_Bmk);
Fcp = mod(Qleft2 * Fcp, 2);

% Step 3
R = Fcp(1:k, m+((k+1):m));
M_R = [eye(k), R; zeros(m-k,k), eye(m-k)];
Qleft3 = Elem1(M_R);
Fcp = mod(Qleft3 * Fcp, 2);

% Step 4
S = Fcp(1:k, m+(1:k));
Pright = Elem2(blkdiag(S,zeros(m-k)));
Fcp = mod(Fcp * Pright, 2);

% Step 5
Fright = mod(Omega * Elem3(k), 2);
Fcp = mod(Fcp * Fright, 2);    % Will be of the form [I, 0; P, I]
P = Fcp(m+(1:m), 1:m);

Q = mod(Qleft3 * Qleft2 * Qleft1, 2);
Decomp = {gf2matinv(Q); Omega; Elem2(P); Elem3(k); gf2matinv(Pright); gf2matinv(Qright)}; 


end