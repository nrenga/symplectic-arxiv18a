% Script to reproduce figures in the paper showing gate complexity

clc
clear 
close all

m_max = 10;  % beyond 10 qubits will take quite some time
ms = 1:m_max;

% Elementary form L_{A_c^{-1}} (L_{A_{\beta}} in the paper)
Alens = zeros(m_max,1);
% Alens = [0; 2; 6; 11; 14; 24; 28; 42; 51; 64; 85; 93; 117; 127; 161; 177];

% Elementary form L_{P^{-1}} (L_{W^{-1}} in the paper)
Plens = zeros(m_max,1);
% Plens = [0; 2; 1; 3; 14; 3; 5; 31; 26; 15; 47; 63; 77; 33; 1; 53];

% Elementary form T_{A_c * P} (T_{A_{\alpha} W} in the paper)
APlens = zeros(m_max,1);
% APlens = [1; 2; 4; 6; 9; 14; 18; 25; 29; 37; 46; 57; 67; 78; 90; 103];

for m = 1:m_max
    disp(m);
    N = 2^m;
    prim_poly = gfprimdf(m);
    field = gftuple([-1:2^m-2]', prim_poly);
    
    p = zeros(1,2*m-1);
    for i = 0:(2*m-2)
%         disp(i);
        p(i+1) = gf2trace(i, m, field);
    end
    P = hankel(p(1:m), p(m:(2*m-1)));
    CP = find_circuit(blkdiag(gf2matinv(P),P'));
    Plens(m) = size(CP,1);
    
    I = eye(m);
    A = I(2:m, :);
    A(m,:) = prim_poly(1,1:m);
    CA = find_circuit(blkdiag(gf2matinv(A), A'));
    Alens(m) = size(CA,1);
    CAP = find_circuit([eye(m), mod(A*P,2); zeros(m), eye(m)]);
    APlens(m) = size(CAP,1);
    fprintf('\n0');
    for i = 1:(2^m-2)
        for j = 1:length(num2str(i-1))
            fprintf('\b');
        end
        fprintf('%d',i);
        Ai = eye(m);
        for j = 1:i
            Ai = mod(Ai*A,2);
        end
        CAi = find_circuit(blkdiag(gf2matinv(Ai), Ai'));
        Alens(m) = max(Alens(m), size(CAi,1));
        CAiP = find_circuit([eye(m), mod(Ai*P,2); zeros(m), eye(m)]);
        APlens(m) = max(APlens(m), size(CAiP,1));
    end
end

fprintf('\n');

figure;
stem(ms,Alens); 
grid on; 
hold on; 
stem(ms,ms.*log2(ms).*log2(log2(ms)))
stem(ms,1.5*ms.*log2(ms).*log2(log2(ms)));
legend({'$$L_{A_{\beta}}$$ worst-case over all $$\beta$$', ...
            '$$m \log_2(m) \log_2(\log_2(m))$$', ...
             '$$1.5 m \log_2(m) \log_2(\log_2(m))$$'}, ...
             'Interpreter', 'latex');

figure;
stem(ms,Plens); 
grid on; 
hold on; 
stem(ms,ms.*log2(ms).*log2(log2(ms)));
legend({'$$L_{W^{-1}}$$', ...
            '$$m \log_2(m) \log_2(\log_2(m))$$'}, ...
             'Interpreter', 'latex');
         
figure;
stem(ms,APlens); 
grid on; 
hold on; 
stem(ms,ms.*log2(ms).*log2(log2(ms)));
legend({'$$T_{A_{\alpha} W}$$ worst-case over all $$\alpha$$', ...
            '$$m \log_2(m) \log_2(\log_2(m))$$'}, ...
             'Interpreter', 'latex');
         
         