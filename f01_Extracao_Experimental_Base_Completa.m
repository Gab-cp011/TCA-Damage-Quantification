%% 01_Extracao_Experimental_Array3D_Validos.m
% Módulo 1: Extração Experimental via Array 3D (Apenas Estados Válidos)
clc; clear; close all;

% =========================================================================
% 1. Definição da Ordem dos Dados e Estados Original
% =========================================================================
% Lista de TODAS as pastas na exata ordem do diretório (e do array 3D)
estados_todos = {'state#01', 'state#02', 'state#03', 'state#04', ...
                 'state#05', 'state#06', 'state#07', 'state#08', ...
                 'state#09', 'state#10', 'state#11', 'state#12', ...
                 'state#13', 'state#14', 'state#15', 'state#16', ...
                 'state#17'};

% Lista dos estados alvos para extração (ordem definida)
estados_validos = {'state#01', 'state#02', 'state#03', 'state#04', 'state#05', ...
                   'state#06', 'state#07', 'state#08', 'state#09'};

ensaios_por_estado = 50;

% =========================================================================
% 2. Carregamento dos Dados e Parâmetros Iniciais
% =========================================================================

% O array 3D possui dimensões [8192 x 5 x 850]:
%  - 8192: Amostras no tempo por ensaio (linhas da matriz 2D)
%  - 5: Canais (1 coluna de força + 4 colunas de aceleração)
%  - 850: Total de ensaios empilhados (fatias / 17 estados * 50 ensaios)

load('data3SS2009.mat', 'dataset'); 
Data3D = dataset; 

% Frequência de amostragem
read_aux = readmatrix('Data/time.txt'); 
t = read_aux(:, 3);  
fs = 1 / (t(2) - t(1)); 

% =========================================================================
% 3. Parâmetros de Processamento de Sinal
% =========================================================================
tam_janela = 1024; 
window = hanning(tam_janela);
noverlap = tam_janela * 3 / 4; 
ch = 4; % Canal de aceleração alvo
nfft = tam_janela; 

Xt = []; 
Labels_t = []; 

f0 = [30; 55; 70];
maxiter = 20;

fprintf('Iniciando extração para %d estados válidos...\n', length(estados_validos));

% =========================================================================
% 4. Loop de Extração Mapeada
% =========================================================================
% Loop 1: Itera apenas sobre os estados que queremos processar (válidos)
for i = 1:length(estados_validos)
    state_name = estados_validos{i};
    
    % Encontra a posição deste estado na matriz 3D original
    idx_no_array = find(strcmp(estados_todos, state_name));
    
    if isempty(idx_no_array)
        error('Estado %s não encontrado na lista de estados totais.', state_name);
    end
    
    % Calcula o range de índices (fatias) correspondentes a este estado,
    % basicamente onde nos 850 dados está localizado o estado atual
    idx_inicio = (idx_no_array - 1) * ensaios_por_estado + 1;
    idx_fim = idx_no_array * ensaios_por_estado;
    
    fprintf('Processando %s (Bloco %d a %d do Array 3D)...\n', state_name, idx_inicio, idx_fim);
    
    % Loop 2: Itera sobre as 50 medições específicas deste estado
    for j = idx_inicio:idx_fim
        raw_data = Data3D(:, :, j);
        
        force = raw_data(:, 1);       
        accels = raw_data(:, 2:5);    
        
        % Estimativa da FRF
        [G, f] = tfestimate(force, accels(:, ch), window, noverlap, nfft, fs);
        
        % Ajuste de Curva 
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
        wn = abs(s) / 2 / pi; % vetor com as  frequências naturais, que vai para o espaço de características
        damp = -real(s) ./ abs(s) * 100; 

        % Concatenação das features (Matriz Alvo)
        Xt = [Xt wn];
        Labels_t = [Labels_t; i];
    end
end

% =========================================================================
% 5. Salvamento da Base de Dados
% =========================================================================
hc_array = 0.0045:-0.0003:0.0015;
nome_arquivo_origem = 'Dados_Alvo_Completos.mat';
save(nome_arquivo_origem, 'Xt', 'Labels_t', 'hc_array');

fprintf('\nExtração concluída com sucesso.\n');
fprintf('Matriz Xt: [%d x %d] (%d features, %d amostras totais)\n', size(Xt,1), size(Xt,2), size(Xt,1), size(Xt,2));

semilogy(f, abs(G), f, abs(Hs))
