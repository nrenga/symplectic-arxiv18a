% Script to translate a given Kerdock unitary 2-design, synthesized via the
% function kerdock_design(m, frobenius, circuits), into a logical unitary
% 2-design on the protected qubits of a given stabilizer code.

% The outputs are cell arrays TG_logical and PG_logical. The full logical
% unitary 2-design is obtained by multiplying each unitary element of
% TG_logical with all unitary elements of PG_logical. So the size of the
% design is the product of the sizes of TG_logical and PG_logical, which is
% same as the sizes of TG and PG, respectively.

% Author: Narayanan Rengaswamy, Date: Sep. 15, 2018

clc
clear
close all
tic;

% Setup the parameters, stabilizers, logical Paulis for the chosen
% stabilizer code
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

code_name = '[[6,4,2]]-CSS';

% Replace the filename below with any file containing cell arrays 'TG',  
% 'PG' and an integer 'm', generated using the 
% kerdock_design(m, frobenius, circuits) function. This value of 'm' 
% corresponds to the number of protected qubits in the stabilizer code used
% in this script, i.e., (m-k) here.

% [TG, PG] = kerdock_design(m-k, 0, 1);
load('TG4_no_frob.mat','TG','PG');
PG_logical = [PG, pauli_group_logical(S, Xbar, Zbar, cell2mat(PG(:,3)))];
TG_logical = TG;

% For codes with large redundancy k (stabilizer dimension), it can take 
% very long for the program to end since for each design element there are 
% 2^{k(k+1)/2} physical realizations and the program chooses the best of 
% them in terms of smallest depth for the circuit.
% If only one physical realization is desired, not necessarily the best,
% then replace the parameter 'all' below for find_logical_cliff with 1.

t1 = clock;
for i = 1:size(TG,1)
    if (mod(i,1000) == 0)
        disp(i);
    end
    F_all = find_logical_cliff(S, Xbar, Zbar, TG{i,4}, [], 'all');
    % Find cheapest circuit in terms of circuit depth
    [depth, cheap_ind] = min(cell2mat(F_all(:,3)));
    F = F_all(cheap_ind, :);
    TG_logical{i,5} = F{1,1};  % symplectic matrix for physical circuit
    TG_logical{i,6} = F{1,2};  % physical circuit
    TG_logical{i,7} = F{1,3};  % depth of circuit
end
t2 = clock;
fprintf('\nElapsed time: %1.4f minutes\n',etime(t2,t1)/60);

% save('TG4_642_no_frob.mat','m','k','S','Xbar','Zbar','TG_logical','PG_logical','t1','t2','code_name');

