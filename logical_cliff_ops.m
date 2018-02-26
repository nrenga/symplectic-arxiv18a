%% Code to get logical Clifford operators for any stabilizer code
% Examples: The [[6,4,2]] CSS code and the [[5,1,3]] perfect code.

% In each cell array produced as output, the first column will contain a 
% symplectic solution and the corresponding second column will contain a
% circuit for that solution.

% For details, please see the paper 
% "Synthesis of Logical Clifford Operators via Symplectic Geometry", 
% available at https://arxiv.org/abs/

% Author: Narayanan Rengaswamy, Date: Feb. 20, 2018

clc
clear
close all
tic;
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

% Phase gate on logical qubit 1
F_all_S1 = find_logical_cliff(S, Xbar, Zbar, 'S', 1);

% Controlled-Z gate on logical qubits (1,2)
F_all_CZ12 = find_logical_cliff(S, Xbar, Zbar, 'CZ', [1 2]);

% CNOT gate where logical qubit 2 controls 1
F_all_CNOT21 = find_logical_cliff(S, Xbar, Zbar, 'CNOT', [2 1]);

% Targeted Hadamard gate on logical qubit 1
F_all_H1 = find_logical_cliff(S, Xbar, Zbar, 'H', 1);

toc;

%% Code to get logical Clifford operators for the [[5,1,3]] perfect code
% clc
% clear
% close all
% tic;
% m = 5;
% k = 4;

% % Stabilizers
% S = [ 1 0 0 1 0 , 0 1 1 0 0 ;
      % 0 1 0 0 1 , 0 0 1 1 0 ;
      % 1 0 1 0 0 , 0 0 0 1 1 ;
      % 0 1 0 1 0 , 1 0 0 0 1 ];

% % Logical Paulis
% Zbar = [ 0 0 0 0 0 , 1 1 1 1 1 ];
% Xbar = [ 1 1 1 1 1 , 0 0 0 0 0 ];

% % Phase gate on the logical qubit
% F_all_S = find_logical_cliff(S, Xbar, Zbar, 'S', 1);

% % Hadamard gate on the logical qubit
% F_all_H = find_logical_cliff(S, Xbar, Zbar, 'H', 1);

% toc;
