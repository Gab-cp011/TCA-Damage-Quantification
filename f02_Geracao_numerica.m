% 02_Geracao_Numerica.m
% Módulo 2: Geração do Domínio de Origem (Modelo Numérico 4-DOF)
clc; clear; close all;

fprintf('Inicializando a modelagem do sistema dinâmico...\n');

% 1. Parâmetros geométricos e propriedades nominais
bc = 0.025; 
hc_normal = 0.006; 
L = 0.177 - 0.025;
I_normal = (bc * hc_normal^3) / 12;

% Massas dos componentes (kg)
mf = 6.26; mc = 0.06; mbumper = 0.26; msp = 0.20; mcf = 0.04;

% Matriz de massa nominal (sem adições operacionais)
m1_nom = mf + mc*4/2 + mcf*4; 
m2_nom = mf + mc*4 + mcf*4; 
m3_nom = mf + mc*4 + mcf*4 + msp; 
m4_nom = mf + mc*4/2 + mcf*4 + mbumper; 

% 2. Mapeamento estocástico do Módulo de Elasticidade
mu_E = 65e9;
sigma_E = 0.15e9;
A_gamma = (mu_E / sigma_E)^2; 
B_gamma = sigma_E^2 / mu_E;

% Semente fixa para garantir a reprodutibilidade da incerteza
rng(42);
amostras_totais = 150 * 9; % 150 iterações para 9 condições
E_samples = gamrnd(A_gamma, B_gamma, [amostras_totais, 1]);

% 3. Configuração da progressão do dano
hc_array = 0.0045:-0.0003:0.0015;
N_rodadas = length(hc_array);

% Estruturas de armazenamento para as 11 rodadas
Xs_todas = cell(N_rodadas, 1);
Labels_todas = cell(N_rodadas, 1);

fprintf('Calculando frequências teóricas para %d rodadas de degradação...\n', N_rodadas);

% 4. Loop principal de simulação
for rodada = 1:N_rodadas
    hc_dano = hc_array(rodada);
    I_dano = (bc * hc_dano^3) / 12;
    
    Xs_atual = zeros(amostras_totais, 3);
    Labels_atual = zeros(amostras_totais, 1);
    amostra_idx = 1;
    
    for condicao = 1:9
        % A. Configuração das variações de massa operacionais
        m1 = m1_nom;
        m2 = m2_nom;
        
        if condicao == 2
            m1 = m1_nom + 1.2; % Massa de 1.2 kg adicionada na base
        elseif condicao == 3
            m2 = m2_nom + 1.2; % Massa de 1.2 kg adicionada no 1º andar
        end
        M = diag([m1, m2, m3_nom, m4_nom]);
        
        % B. Configuração das perdas de rigidez
        mult1 = 1; mult2 = 1; mult3 = 1;
        fator_1col = (3 * I_normal + 1 * I_dano) / (4 * I_normal);
        fator_2col = (2 * I_normal + 2 * I_dano) / (4 * I_normal);
        
        if condicao == 4, mult1 = fator_1col; end
        if condicao == 5, mult2 = fator_1col; end
        if condicao == 6, mult3 = fator_1col; end
        if condicao == 7, mult1 = fator_2col; end
        if condicao == 8, mult2 = fator_2col; end
        if condicao == 9, mult3 = fator_2col; end
        
        % C. Resolução do problema de autovalores
        for iter = 1:150
            E_atual = E_samples(amostra_idx);
            k_base = 4 * 12 * E_atual * I_normal / L^3;
            ka = 10; % Rigidez residual dos trilhos de deslizamento
            
            k1 = k_base * 0.5 * mult1;
            k2 = k_base * mult2;
            k3 = k_base * mult3;
            
            K = [k1+ka, -k1,      0,        0; 
                -k1,    k1+k2,   -k2,       0; 
                 0,    -k2,       k2+k3,   -k3; 
                 0,     0,       -k3,       k3];
            
            [~, D] = eig(K, M);
            freq = sqrt(diag(D)) / (2*pi);
            
            % D. Seleção dos Modos de Flexão
            % Os índices 2 a 4 ignoram a translação de corpo rígido (~0.1 Hz)
            Xs_atual(amostra_idx, :) = freq(2:4)';
            Labels_atual(amostra_idx) = condicao; 
            
            amostra_idx = amostra_idx + 1;
        end
    end
    
    Xs_todas{rodada} = Xs_atual;
    Labels_todas{rodada} = Labels_atual;
end

% 5. Salvamento da base de dados teórica
nome_arquivo_origem = 'Dados_Origem.mat';
save(nome_arquivo_origem, 'Xs_todas', 'Labels_todas', 'hc_array');
fprintf('\nFase 2 concluída. Matrizes de origem consolidadas em: %s\n', nome_arquivo_origem);