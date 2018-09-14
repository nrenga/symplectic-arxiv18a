% Script to translate a given Kerdock unitary 2-design, synthesized via the
% function KerdockDesign(m, frobenius, circuits), into a logical unitary
% 2-design on the protected qubits of a given stabilizer code.

% Author: Narayanan Rengaswamy, Date: Sep. 14, 2018

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

% Replace the filename below with any file containing a cell array 'TG' and 
% an integer 'm', generated using the KerdockDesign(m, frobenius, circuits)
% function. This value of 'm' corresponds to the number of protected qubits
% in the stabilizer code used in this script, i.e., (m-k) here.
load('TG4_with_frob.mat');
TG_logical = TG;
t1 = clock;
for i = 1:size(TG,1)
    if (mod(i,1000) == 0)
        disp(i);
    end
    opt = 5;
    if (size(TG{i,4},1) < size(TG{i,5},1))
        opt = 4;
    end
    F_all = find_logical_cliff(S, Xbar, Zbar, TG{i,opt}, [], 'all');
    % Find cheapest circuit in terms of circuit depth
    [depth, cheap_ind] = min(cell2mat(F_all(:,3)));
    F = F_all(cheap_ind, :);
    TG_logical{i,6} = F{1,1};
    TG_logical{i,7} = F{1,2};
    TG_logical{i,8} = F{1,3};
end
t2 = clock;
fprintf('\nElapsed time: %1.4f minutes\n',etime(t2,t1)/60);
