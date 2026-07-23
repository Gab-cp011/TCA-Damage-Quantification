% 02_Geracao_Numerica.m
% Módulo 2: Geração do Domínio de Origem (Modelo Numérico 4-DOF)
clc; clear; close all;

fprintf('Inicializando a modelagem do sistema dinâmico...\n');

% 1. Parâmetros geométricos e propriedades nominais
bc = 0.025; % largura da coluna
hc_normal = 0.006; % espessura da coluna íntegra
L = 0.152; % comprimento da coluna
I_normal = (bc * hc_normal^3) / 12; 

% Massas dos componentes (kg)
mf = 6.26; % massa da placa de alumínio de cada andar
mc = 0.06; % massa da coluna de alumínio
mbumper = 0.20; % massa do bumper
msp = 0.26; % massa da coluna suspensa do mecanismo de dano 
mcf = 0.04; % massa do bloco de parafuso da coluna

% Matriz de massa nominal (sem adições operacionais)
m1_nom = mf + mc*4/2 + mcf*4; % massa da base
m2_nom = mf + mc*4 + mcf*4; % massa do primeiro andar
m3_nom = mf + mc*4 + mcf*4 + mbumper; % massa do segundo andar
m4_nom = mf + mc*4/2 + mcf*4 + msp; % massa do terceiro andar

% Razões de amortecimento extraídas do modelo experimental (dcv)
dcv = [0.0001; 0.063; 0.020; 0.0097];

% 2. Mapeamento estocástico do Módulo de Elasticidade
mu_E = 65e9; % média assumida do material
sigma_E = 0.15e9; % desvio padrão

% parâmetros da distribuição Gamma
A_gamma = (mu_E / sigma_E)^2;  
B_gamma = sigma_E^2 / mu_E;

% Semente fixa para garantir a reprodutibilidade da incerteza
rng(42);
amostras_totais = 150 * 9; % 150 iterações para 9 condições
% sorteia os valores de rigidez do material em cada simulação
E_samples = gamrnd(A_gamma, B_gamma, [amostras_totais, 1]);

% 3. Configuração da progressão do dano
hc_array = 0.0045:-0.0003:0.0015;
N_rodadas = length(hc_array);

% Estruturas de armazenamento para as 11 rodadas, representando a redução
% de espessura gradativa
Xs_todas = cell(N_rodadas, 1);
Labels_todas = cell(N_rodadas, 1);

fprintf('Calculando frequências teóricas para %d rodadas de degradação...\n', N_rodadas);

% 4. Loop principal de simulação
for rodada = 1:N_rodadas
    hc_dano = hc_array(rodada);
    I_dano = (bc * hc_dano^3) / 12; %  inércia de coluna com dano
    
    % Fatores de redução de rigidez (percentual entre a rigidez reduzida e
    % a rigidez "inicial")
    fator_1col = (3 * I_normal + 1 * I_dano) / (4 * I_normal);
    fator_2col = (2 * I_normal + 2 * I_dano) / (4 * I_normal);
    
    Xs_atual = zeros(amostras_totais, 3);
    Labels_atual = zeros(amostras_totais, 1);
    amostra_idx = 1;
    
    for condicao = 1:9
        % A. Vetores de propriedades base (Condição 1: Saudável)
        m_andares = [m1_nom, m2_nom, m3_nom, m4_nom];
        mult_rigidez = [1.0, 1.0, 1.0]; % [1º andar, 2º andar, 3º andar]
        
        % B. Mapeamento dos Cenários Experimentais
        switch condicao
            % Variações Operacionais (Massa adicionada)
            case 2
                m_andares(1) = m_andares(1) + 1.2; % +1.2 kg na base
            case 3
                m_andares(2) = m_andares(2) + 1.2; % +1.2 kg no 1º andar
                
            % Dano Nível 1 (Troca de 1 coluna)
            case 4
                mult_rigidez(1) = fator_1col; % Dano 1 no 1º andar
            case 5
                mult_rigidez(2) = fator_1col; % Dano 1 no 2º andar
            case 6
                mult_rigidez(3) = fator_1col; % Dano 1 no 3º andar
                
            % Dano Nível 2 (Troca de 2 colunas)
            case 7
                mult_rigidez(1) = fator_2col; % Dano 2 no 1º andar
            case 8
                mult_rigidez(2) = fator_2col; % Dano 2 no 2º andar
            case 9
                mult_rigidez(3) = fator_2col; % Dano 2 no 3º andar
        end
        
        % Consolidação da matriz de massa do cenário atual
        M = diag(m_andares);
        
        % C. Resolução do problema de autovalores
        for iter = 1:150
            E_atual = E_samples(amostra_idx);
            
            % Rigidez base usando a amostra estocástica de E
            k_base = 4 * 12 * E_atual * I_normal / L^3;
            ka = 10; % Atrito residual dos trilhos
            
            % Aplicação direta dos multiplicadores de rigidez
            k1 = k_base * mult_rigidez(1);
            k2 = k_base * mult_rigidez(2);
            k3 = k_base * mult_rigidez(3);
            
            K = [k1+ka, -k1,      0,        0; 
                -k1,    k1+k2,   -k2,       0; 
                 0,    -k2,       k2+k3,   -k3; 
                 0,     0,       -k3,       k3];
            
            % 1. Problema de autovalores não amortecido (intermediário)
            [V, D_und] = eig(K, M);
            wn = sqrt(diag(D_und));
            
            % 2. Cálculo da massa modal
            Mn = diag(V' * M * V);
            
            % 3. Matriz de amortecimento modal e reconstrução da matriz C física
            A_modal = diag(2 .* dcv .* wn .* Mn);
            C= (V'\A_modal)/V;
            
            % 4. Montagem da Matriz de Estado global (SSA)
            SSA = [zeros(4), eye(4); -M\K, -M\C];
            
            % 5. Solução do problema de autovalores amortecidos
            autovalores = eig(SSA);
            
            % 6. Filtragem das raízes (apenas conjugados com parte imaginária positiva)
            lambda = autovalores(imag(autovalores) > 0);
            
            % A frequência registrada é a parte imaginária do autovalor
            % (frequências naturais amortecidas)
            freq_amortecida = sort(imag(lambda) / (2*pi));
            
            % D. Seleção dos Modos de Flexão
            % Os índices 2 a 4 ignoram a translação de corpo rígido (~0.1 Hz)
            Xs_atual(amostra_idx, :) = freq_amortecida(2:4)';
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
