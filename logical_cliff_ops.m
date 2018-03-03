%% Code to get logical Clifford operators for any stabilizer code
% Examples: The [[6,4,2]] CSS code and the [[5,1,3]] perfect code.

% In each cell array produced as output, the first column will contain a 
% symplectic solution, the corresponding second column will contain a
% circuit for that solution, and the third column gives the circuit depth.

% For details, please see the paper 
% "Synthesis of Logical Clifford Operators via Symplectic Geometry", 
% available at https://arxiv.org/abs/

% Author: Narayanan Rengaswamy, Date: Mar. 3, 2018

clc
clear
close all
tic;
m = 6;
k = 2;  % implies 2^(k*(k+1)/2) solutions for a logical Clifford circuit!

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

% Phase gate on logical qubit 1
F_all_S1 = find_logical_cliff(S, Xbar, Zbar, {'S', 1}, [], 'all');
% Find cheapest circuit in terms of circuit depth
[S1depth, S1cheap_ind] = min(cellfun(@(a) a(1,1), F_all_S1(:,3)));
F_S1 = F_all_S1(S1cheap_ind, :);

% Controlled-Z gate on logical qubits (1,2)
F_all_CZ12 = find_logical_cliff(S, Xbar, Zbar, {'CZ', [1 2]}, [], 'all');
% Find cheapest circuit
[CZ12depth, CZ12cheap_ind] = min(cellfun(@(a) a(1,1), F_all_CZ12(:,3)));
F_CZ12 = F_all_CZ12(CZ12cheap_ind, :);

% CNOT gate where logical qubit 2 controls 1
F_all_CNOT21 = find_logical_cliff(S, Xbar, Zbar, {'CNOT', [2 1]}, [], 'all');
% Find cheapest circuit
[CNOT21depth, CNOT21cheap_ind] = min(cellfun(@(a) a(1,1), F_all_CNOT21(:,3)));
F_CNOT21 = F_all_CNOT21(CNOT21cheap_ind, :);

% Targeted Hadamard gate on logical qubit 1
% Snorm = [S(1,:) + S(2,:); S(2,:)];
Snorm = S;
F_all_H1 = find_logical_cliff(S, Xbar, Zbar, {'H', 1}, Snorm, 'all');
% Find cheapest circuit
[H1depth, H1cheap_ind] = min(cellfun(@(a) a(1,1), F_all_H1(:,3)));
F_H1 = F_all_H1(H1cheap_ind, :);

% Transversal Hadamard gate
Snorm = [S(2,:); S(1,:)];
F_all_Htrans = find_logical_cliff(S, Xbar, Zbar, {'H', 1:4}, Snorm, 'all');
% Find cheapest circuit
[Htransdepth, Htranscheap_ind] = min(cellfun(@(a) a(1,1), F_all_Htrans(:,3)));
F_Htran = F_all_Htrans(Htranscheap_ind, :);

% Logical Clifford circuit
Snorm = S;
% Snorm = [S(1,:) + S(2,:); S(2,:)];
% Snorm = [S(2,:); S(1,:)];
ckt1 = {'S', 2; 'CNOT', [1 3]; 'CZ', [3 4]};
F_all_ckt1 = find_logical_cliff(S, Xbar, Zbar, ckt1, Snorm, 'all');
% Find cheapest circuit
[ckt1depth, ckt1cheap_ind] = min(cellfun(@(a) a(1,1), F_all_ckt1(:,3)));
F_ckt1 = F_all_Htrans(ckt1cheap_ind, :);

toc;

%% Code to get logical Clifford operators for the [[5,1,3]] perfect code
% clc
% clear
% close all
% tic;
% m = 5;
% k = 4;  % implies 2^(k*(k+1)/2) solutions for a logical Clifford circuit!
% 
% % Stabilizers
% S = [ 1 0 0 1 0 , 0 1 1 0 0 ;
%       0 1 0 0 1 , 0 0 1 1 0 ;
%       1 0 1 0 0 , 0 0 0 1 1 ;
%       0 1 0 1 0 , 1 0 0 0 1 ];
% 
% % Logical Paulis
% Zbar = [ 0 0 0 0 0 , 1 1 1 1 1 ];
% Xbar = [ 1 1 1 1 1 , 0 0 0 0 0 ];
% 
% % Phase gate on the logical qubit
% F_all_S = find_logical_cliff(S, Xbar, Zbar, {'S', 1}, [], 'all');
% % Find cheapest circuit
% [Sdepth, Scheap_ind] = min(cellfun(@(a) a(1,1), F_all_S(:,3)));
% F_S = F_all_S(Scheap_ind, :);
% 
% % Hadamard gate on the logical qubit
% F_all_H = find_logical_cliff(S, Xbar, Zbar, {'H', 1}, [], 'all');
% % Find cheapest circuit
% [Hdepth, Hcheap_ind] = min(cellfun(@(a) a(1,1), F_all_H(:,3)));
% F_H = F_all_H(Hcheap_ind, :);
% 
% toc;

%% The [[15,7,3]] Hamming code

% clc
% clear
% close all
% 
% m = 15;
% k = 8;  % implies 2^(k*(k+1)/2) solutions for a logical Clifford circuit!
% 
% % Parity-check matrix for the [15,11,3] binary Hamming code
% H = de2bi((1:15)',4,'left-msb')';
% S = [ H, zeros(4,m);
%       zeros(4,m), H ];
% 
% L = [ 1 1 0 1 0 0 0 1 0 0 0 0 0 0 1 ;...
%       1 1 0 0 1 0 0 0 0 1 0 1 0 0 0 ;...
%       1 1 0 0 0 1 0 0 0 0 1 0 0 1 0 ;...
%       1 1 0 0 0 0 1 0 1 0 0 0 1 0 0 ;...
%       1 0 0 1 0 1 0 0 1 1 0 0 0 0 0 ;...
%       1 0 0 1 0 0 1 0 0 0 0 1 0 1 0 ;...
%       1 0 0 0 0 0 0 1 0 1 0 0 1 1 0 ];    % Logical Xs and Zs for [[15,7,3]]
% 
% Xbar = [ L, zeros(7,m) ];
% Zbar = [ zeros(7,m), L ];
% 
% Snorm = S;
% ckt1 = {'Permute', [2 3 1 4:7]};
% F_all_ckt1 = find_logical_cliff(S, Xbar, Zbar, ckt1, Snorm, 1);
% % Find cheapest circuit
% % [ckt1depth, ckt1cheap_ind] = min(cellfun(@(a) a(1,1), F_all_ckt1(:,3)));
% % F_ckt1 = F_all_Htrans(ckt1cheap_ind, :);

