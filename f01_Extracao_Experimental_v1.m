%% 01_Extracao_Experimental_Atualizado.m
% Módulo 1: Extração Experimental via PSD ou Circle Fit Modal
clc; clear; close all;

% 1. Definição dos parâmetros iniciais
data_dir = 'Data';
estados_validos = {'state#13', 'state#02', 'state#01', 'state#17', ...
    'state#24', 'state#21', 'state#18', 'state#23', 'state#22'};
fs = 322.58; 

% =========================================================================
% CONFIGURAÇÕES
% =========================================================================
% nfft sugerido: 4096 para manter boa resolução sem ruído excessivo
nfft = 4096; 
window = hanning(nfft);
noverlap = nfft / 2;

% Escolha: 'unico', 'multicanal' (PSD) ou 'circle_fit' (FRF)
metodo_extracao = 'circle_fit'; 
metodo_fusao = 'media'; % Aplicável se metodo_extracao for 'multicanal'
% =========================================================================

Xt = []; Labels_t = [];
bounds = [27.0, 33.0; 51.0, 59.0; 67.0, 73.0];

fprintf('Iniciando processamento com método: %s...\n', metodo_extracao);

for i = 1:length(estados_validos)
    state_name = estados_validos{i};
    arquivos = dir(fullfile(data_dir, state_name, 'data*.txt'));

    for j = 1:length(arquivos)
        raw_data = readmatrix(fullfile(arquivos(j).folder, arquivos(j).name));

        % Força (Input): Coluna 1 | Acelerômetros (Output): Colunas 2 a 5
        force = raw_data(:, 1);
        accels = raw_data(:, 2:5);

        % --- BLOCO ALTERNÁVEL: MÉTODO DE EXTRAÇÃO ---
        if strcmp(metodo_extracao, 'unico')
            [Pxx, f_psd] = pwelch(accels(:, 4), window, noverlap, nfft, fs);

        elseif strcmp(metodo_extracao, 'multicanal')
            Pxx_all = zeros(nfft/2 + 1, 4);
            for canal = 1:4
                [Pxx_all(:, canal), f_psd] = pwelch(accels(:, canal), window, noverlap, nfft, fs);
            end
            Pxx = (strcmp(metodo_fusao, 'soma')) ? sum(Pxx_all, 2) : mean(Pxx_all, 2);

        elseif strcmp(metodo_extracao, 'circle_fit')
            % Cálculo da FRF Complexa (H1 Estimator)
            H_all = zeros(nfft/2 + 1, 4);
            for canal = 1:4
                [H_all(:, canal), f_frf] = tfestimate(force, accels(:, canal), window, noverlap, nfft, fs);
            end
            % Média complexa das FRFs (Preserva fase e magnitude)
            H_complex = mean(H_all, 2);
        end
        % --------------------------------------------

        fn = zeros(1, 3);
        for k = 1:3
            if strcmp(metodo_extracao, 'circle_fit')
                % Encontra o pico aproximado dentro do range para o circle fit
                band_idx = (f_frf >= bounds(k,1)) & (f_frf <= bounds(k,2));
                [~, max_idx] = max(abs(H_complex(band_idx)));
                f_band_indices = find(band_idx);
                peak_f = f_frf(f_band_indices(max_idx));

                % Executa o Circle Fit
                fn(k) = circle_fit_modal(f_frf, H_complex, peak_f);
            else
                % Método clássico PSD (unico ou multicanal)
                idx_range = f_psd >= bounds(k,1) & f_psd <= bounds(k,2);
                [~, max_idx] = max(Pxx(idx_range));
                f_range = f_psd(idx_range);
                fn(k) = f_range(max_idx);
            end
        end
        Xt = [Xt; fn];
        Labels_t = [Labels_t; i];
    end
end

save('Dados_Alvo.mat', 'Xt', 'Labels_t');
fprintf('Extração concluída.\n');