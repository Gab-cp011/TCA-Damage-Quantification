% 03_Adaptacao_Classificacao.m
% Módulo 3: Adaptação de Domínio (TCA) e Classificação (Mahalanobis)
clc; clear; close all;

% 1. Carregamento unificado dos dados
arquivo_alvo = 'Dados_Alvo.mat';
arquivo_origem = 'Dados_Origem.mat';

if isfile(arquivo_alvo) && isfile(arquivo_origem)
    load(arquivo_alvo, 'Xt', 'Labels_t');
    load(arquivo_origem, 'Xs_todas', 'Labels_todas', 'hc_array');
else
    error('Arquivos não encontrados. Execute as Fases 1 e 2 previamente.');
end

Xt=Xt'; 

N_rodadas = length(hc_array);
historico_macro_f1 = zeros(N_rodadas, 1);

% Binarização dos rótulos reais: condições 1 a 3 são saudáveis, 4 a 9 representam dano
Labels_t_binario = Labels_t > 3; 

fprintf('Iniciando alinhamento e classificação para %d rodadas...\n\n', N_rodadas);

for i = 1:N_rodadas
    Xs = Xs_todas{i};
    Labels_s = Labels_todas{i};

    % 2. Normalização Z-score ancorada no domínio de origem
    % Preserva a separabilidade não supervisionada do domínio alvo
    mu_X = mean(Xs);
    sigma_X = std(Xs);
    
    Xs_norm = (Xs - mu_X) ./ sigma_X;
    Xt_norm = (Xt - mu_X) ./ sigma_X;

    %Xs_norm = zscore(Xs);
    %Xt_norm = zscore(Xt); 

    % 3. Análise de Componentes de Transferência
    % Configuração replicada do artigo: 1 componente principal, regularização 0.1
    [Zs, Zt] = TCA_linear(Xs_norm, Xt_norm, 1, 0.1);

    % 4. Classificação de Mahalanobis
    % A fronteira aprende exclusivamente com as amostras intactas (condições 1 a 3)
    Zs_saudavel = Zs(Labels_s <= 3, :);
    mu_s = mean(Zs_saudavel);
    Sigma_s = cov(Zs_saudavel);

    % Estabilização da matriz de covariância
    Sigma_s = Sigma_s + 1e-8 * eye(size(Sigma_s));

    % Limiar percentílico (99%) para blindagem contra outliers estocásticos
    dist_saudavel = pdist2(Zs_saudavel, mu_s, 'mahalanobis', Sigma_s).^2;

    %limiar = prctile(dist_saudavel, 99);
    limiar = max(dist_saudavel);

    % Avaliação da integridade do domínio físico
    dist_target = pdist2(Zt, mu_s, 'mahalanobis', Sigma_s).^2;
    previsoes_target = dist_target > limiar;

    % 5. Avaliação de Performance (Macro F1)
    tp = sum(previsoes_target == 1 & Labels_t_binario == 1);
    fp = sum(previsoes_target == 1 & Labels_t_binario == 0);
    fn = sum(previsoes_target == 0 & Labels_t_binario == 1);
    tn = sum(previsoes_target == 0 & Labels_t_binario == 0);

    f1_danificado = 0; 
    f1_saudavel = 0;

    if (2*tp + fn + fp) > 0
        f1_danificado = 2*tp / (2*tp + fn + fp);
    end
    if (2*tn + fn + fp) > 0
        f1_saudavel = 2*tn / (2*tn + fn + fp);
    end

    historico_macro_f1(i) = (f1_danificado + f1_saudavel) / 2;

    fprintf('Rodada %02d (hc = %.4f m) | Macro F1: %.4f\n', i, hc_array(i), historico_macro_f1(i));
end

nome_arquivo_origem = 'Dados_Espaco_Latente.mat';
save(nome_arquivo_origem, 'Zs', 'Zt');
fprintf('\nFase 3 concluída.');
