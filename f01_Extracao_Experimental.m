%% 01_Extracao_Experimental.m
% Módulo 1: Extração Experimental (PSD única, PSD multicanal ou Circle Fit)
clc; clear; close all;

% 1. Parâmetros iniciais e diretórios
data_dir = 'Data';
estados_validos = {'state#13', 'state#02', 'state#01', 'state#17', ...
    'state#24', 'state#21', 'state#18', 'state#23', 'state#22'};
fs = 322.58; 

% =========================================================================
% CHAVEAMENTO METODOLÓGICO
% =========================================================================
nfft = 4096; 
window = hanning(nfft);
noverlap = nfft / 2;

% Opções: 'unico', 'multicanal', 'circle_fit'
metodo_extracao = 'circle_fit'; 
% Opções: 'soma', 'media' (usado apenas se metodo_extracao = 'multicanal')
metodo_fusao = 'media'; 
% =========================================================================

Xt = []; Labels_t = [];
bounds = [27.0, 33.0; 51.0, 59.0; 67.0, 73.0];

fprintf('Iniciando extração via método: %s\n', metodo_extracao);

for i = 1:length(estados_validos)
    state_name = estados_validos{i};
    arquivos = dir(fullfile(data_dir, state_name, 'data*.txt'));
    
    if isempty(arquivos)
        error('Arquivos não encontrados no diretório %s.', state_name);
    end
    
    for j = 1:length(arquivos)
        raw_data = readmatrix(fullfile(arquivos(j).folder, arquivos(j).name));
        
        % Separação das matrizes de entrada (força) e saída (aceleração)
        force = raw_data(:, 1);
        accels = raw_data(:, 2:5);

        % Lógica de roteamento de extração
        if strcmp(metodo_extracao, 'unico')
            [Pxx, f_vec] = pwelch(accels(:, 4), window, noverlap, nfft, fs);
            
        elseif strcmp(metodo_extracao, 'multicanal')
            Pxx_all = zeros(nfft/2 + 1, 4);
            for canal = 1:4
                [Pxx_all(:, canal), f_vec] = pwelch(accels(:, canal), window, noverlap, nfft, fs);
            end
            
            % Estrutura de decisão validada para fusão
            if strcmp(metodo_fusao, 'soma')
                Pxx = sum(Pxx_all, 2);
            else
                Pxx = mean(Pxx_all, 2);
            end
            
        elseif strcmp(metodo_extracao, 'circle_fit')
            H_all = zeros(nfft/2 + 1, 4);
            for canal = 1:4
                % H1 Estimator para isolar polos mecânicos reais
                [H_all(:, canal), f_vec] = tfestimate(force, accels(:, canal), window, noverlap, nfft, fs);
            end
            H_complex = mean(H_all, 2);
        else
            error('Método de extração inválido configurado.');
        end

        % Extração paramétrica por banda
        fn = zeros(1, 3);
        for k = 1:3
            idx_range = f_vec >= bounds(k,1) & f_vec <= bounds(k,2);
            
            if strcmp(metodo_extracao, 'circle_fit')
                f_band_indices = find(idx_range);
                
                % Identificação do pico preliminar em magnitude para ancorar a busca de Nyquist
                [~, max_local_idx] = max(abs(H_complex(idx_range)));
                peak_f = f_vec(f_band_indices(max_local_idx));
                
                % Acionamento da extração de fase
                fn(k) = circle_fit_modal(f_vec, H_complex, peak_f);
            else
                f_range = f_vec(idx_range);
                [~, max_idx] = max(Pxx(idx_range));
                fn(k) = f_range(max_idx);
            end
        end
        
        Xt = [Xt; fn];
        Labels_t = [Labels_t; i];
    end
end

save('Dados_Alvo.mat', 'Xt', 'Labels_t');
fprintf('Extração física finalizada. Matriz salva para etapa de Adaptação de Domínio.\n');
