%% 01_Extracao_Experimental.m
% Módulo 1: Extração Experimental (PSD única, PSD multicanal ou FRF Fitting)
clc; clear; close all;

% =========================================================================
% 1. Parâmetros Iniciais e Diretórios
% =========================================================================
data_dir = 'Data';
estados_validos = {'state#13', 'state#02', 'state#01', 'state#17', ...
                   'state#24', 'state#21', 'state#18', 'state#23', 'state#22'};

read_aux = readmatrix('Data/time.txt');
t = read_aux(:, 3);
fs = 1 / (t(2) - t(1)); 

% =========================================================================
% 2. Parâmetros de Processamento de Sinal
% =========================================================================
% Configurações de janela temporal (Redução de ruído via médias)
tam_janela = 1024 / 4; 
window = hanning(tam_janela);
noverlap = tam_janela * 3 / 4; % Sobreposição de 75%
ch = 4; % Canal de aceleração alvo

% Configuração de Resolução Espectral (Zero-Padding)
nfft = tam_janela; 

Xt = []; 
Labels_t = [];

% Estimativas iniciais para o algoritmo de Fitting
f0 = [30; 55; 70];
maxiter = 20;

fprintf('Iniciando extração de características para o canal %d\n', ch);

% =========================================================================
% 3. Loop de Extração de Parâmetros Modais
% =========================================================================
for i = 1:length(estados_validos)
    state_name = estados_validos{i};
    arquivos = dir(fullfile(data_dir, state_name, 'data*.txt'));
    
    if isempty(arquivos)
        error('Arquivos não encontrados no diretório %s.', state_name);
    end
    
    for j = 1:length(arquivos)
        raw_data = readmatrix(fullfile(arquivos(j).folder, arquivos(j).name));
        
        % Separação dos canais
        force = raw_data(:, 1);
        accels = raw_data(:, 2:5);

        % Estimativa da Função de Resposta em Frequência (FRF)
        [G, f] = tfestimate(force, accels(:, ch), window, noverlap, nfft, fs);
    
        % Ajuste de Curva para extração dos parâmetros modais
        [q, err] = Fitting(G, f, f0, maxiter);
    
        u = q(1);
        v = q(2);
        r = q(3:5);
        s = q(6:8);
        
        % Reconstrução da FRF analítica 
        Hs = zeros(size(G));
        W = 2 * pi * f;
        for k = 1:length(r)
            Hs = Hs + r(k) ./ (1i * W - s(k));
        end
        Hs = Hs + u + 1i * (W - mean(W)) * v;
    
        % Extração da Frequência Natural Amortecida (polos do sistema)
        wn = abs(s) / 2 / pi;
        damp = -real(s) ./ abs(s) * 100;

        % Concatenação das features (Matriz Alvo)
        Xt = [Xt wn];
        Labels_t = [Labels_t; i];
    end
end

% =========================================================================
% 4. Salvamento da Base de Dados
% =========================================================================
hc_array = 0.0045:-0.0003:0.0015;
nome_arquivo_origem = 'Dados_Alvo.mat';

save(nome_arquivo_origem, 'Xt', 'Labels_t', 'hc_array');
fprintf('\nFase 1 concluída. Matrizes de alvo consolidadas em: %s\n', nome_arquivo_origem);
