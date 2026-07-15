%% 01_Extracao_Experimental.m
% Módulo 1: Extração Experimental via Densidade Espetral de Potência (Welch)
clc; clear; close all;

% 1. Definição dos parâmetros iniciais e diretórios
data_dir = 'Data';
estados_validos = {'state#13', 'state#02', 'state#01', 'state#17', ...
    'state#24', 'state#21', 'state#18', 'state#23', 'state#22'};
fs = 322.58; % Frequência de amostragem definida pelo equipamento (Hz)

% =========================================================================
% CONFIGURAÇÕES DE EXTRAÇÃO E PROCESSAMENTO
% =========================================================================
% nfft ajustado para 4096 para equilibrar resolução e média temporal
nfft = 4096; 
window = hanning(nfft);
noverlap = nfft / 2; % Sobreposição de 50% entre os blocos

% Escolha o método de extração: 'unico' (artigo original) ou 'multicanal'
metodo_extracao = 'multicanal'; 

% Se for multicanal, escolha como fundir as PSDs: 'media' ou 'soma'
metodo_fusao = 'media'; 
% =========================================================================

% Inicialização das matrizes de armazenamento finais
Xt = [];
Labels_t = [];

% Fronteiras de contenção para os três modos de flexão
bounds = [27.0, 33.0;  
    51.0, 59.0;  
    67.0, 73.0];

fprintf('Iniciando processamento dos sinais experimentais...\n');

% 3. Processo iterativo de leitura e extração
for i = 1:length(estados_validos)
    state_name = estados_validos{i};
    state_path = fullfile(data_dir, state_name);

    % Busca todos os arquivos de texto dentro da pasta do estado atual
    arquivos = dir(fullfile(state_path, 'data*.txt'));

    if isempty(arquivos)
        error('Arquivos não encontrados para %s. Verifique o caminho.', state_name);
    end

    num_testes = length(arquivos);

    for j = 1:num_testes
        file_path = fullfile(arquivos(j).folder, arquivos(j).name);
        raw_data = readmatrix(file_path);

        % --- BLOCO ALTERNÁVEL: EXTRAÇÃO DA PSD ---
        if strcmp(metodo_extracao, 'unico')
            % Método 1: Apenas aceleração do topo (Canal 4 = coluna 5)
            sinal_topo = raw_data(:, 5);
            
            % Cálculo da Densidade Espetral de Potência (PSD)
            [Pxx, f_psd] = pwelch(sinal_topo, window, noverlap, nfft, fs);
            
        elseif strcmp(metodo_extracao, 'multicanal')
            % Método 2: Fusão dos 4 acelerômetros (Canais 1 a 4 = colunas 2 a 5)
            sinais_andares = raw_data(:, 2:5);
            
            % Inicializa a matriz para guardar as PSDs individuais dos 4 andares
            % O tamanho de saída do f_psd no MATLAB para sinais reais é (nfft/2)+1
            Pxx_all = zeros(nfft/2 + 1, 4);
            
            for canal = 1:4
                [Pxx_all(:, canal), f_psd] = pwelch(sinais_andares(:, canal), window, noverlap, nfft, fs);
            end
            
            % Fusão espacial da energia espectral
            if strcmp(metodo_fusao, 'soma')
                Pxx = sum(Pxx_all, 2);
            elseif strcmp(metodo_fusao, 'media')
                Pxx = mean(Pxx_all, 2);
            else
                error('Método de fusão inválido. Escolha "soma" ou "media".');
            end
        else
            error('Método de extração inválido. Escolha "unico" ou "multicanal".');
        end
        % -----------------------------------------

        fn = zeros(1, 3);
        for k = 1:3
            % Isolamento da região de busca para o modo específico
            idx_range = f_psd >= bounds(k,1) & f_psd <= bounds(k,2);
            f_range = f_psd(idx_range);
            Pxx_range = Pxx(idx_range);

            % Identificação da frequência com maior concentração de energia
            [~, max_idx] = max(Pxx_range);
            fn(k) = f_range(max_idx);
        end

        % Acumulação das três frequências e do rótulo numérico da condição
        Xt = [Xt; fn(1), fn(2), fn(3)];
        Labels_t = [Labels_t; i];
    end
    fprintf('Estado processado com sucesso: %s\n', state_name);
end

% 4. Salvamento seguro da base de dados física
nome_arquivo_alvo = 'Dados_Alvo.mat';
save(nome_arquivo_alvo, 'Xt', 'Labels_t');
fprintf('\nFase 1 concluída. Dados experimentais consolidados em: %s\n', nome_arquivo_alvo);