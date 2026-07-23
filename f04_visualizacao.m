% =========================================================================
% MÓDULO DE VISUALIZAÇÃO - CORRIGIDO
% =========================================================================

load('Dados_Alvo.mat');    % Carrega Xt e Labels_t
load('Dados_Origem.mat');  % Carrega Xs e Labels_s

% 1. Tratamento Prévio dos Dados (Correção de Bugs)
% Transpor Xt para o formato N x 3 (amostras x características)
if size(Xt, 1) == 3
    Xt = Xt'; 
end

% Extrair a primeira rodada de degradação do modelo numérico para plotagem
% (Você pode alterar o índice {1} para visualizar outras rodadas)
if iscell(Xs_todas)
    Xs_mat = Xs_todas{1}; 
    Labels_s_mat = Labels_todas{1};
else
    Xs_mat = Xs_todas;
    Labels_s_mat = Labels_todas;
end

% Definição da paleta de cores
mapa_cores = [0 0 1; 1 0 0; 0 0 0; 1 0 1; 0 0.6 0; 0.5 0.5 0.5; 0.8 0.4 0; 0 0.7 1; 1 0.8 0];
cores_dano = [0 0 1; 1 0 0; 0 0.6 0]; 

% -------------------------------------------------------------------------
% FIGURA 4: Espaço 3D Original (Frequências Naturais)
% -------------------------------------------------------------------------
fig4 = figure('Name', 'Figura 4 - Espaço de Características 3D', 'Color', 'w');
hold on; grid on;
for c = 1:9
    idx_s = (Labels_s_mat == c);
    idx_t = (Labels_t == c);

    % Origem (Triângulos) e Alvo (Cruzes) usando as matrizes corrigidas
    scatter3(Xs_mat(idx_s, 1), Xs_mat(idx_s, 2), Xs_mat(idx_s, 3), 40, mapa_cores(c,:), '^', 'filled');
    scatter3(Xt(idx_t, 1), Xt(idx_t, 2), Xt(idx_t, 3), 40, mapa_cores(c,:), 'x', 'LineWidth', 1.2);
end
xlabel('f_1 [Hz]'); ylabel('f_2 [Hz]'); zlabel('f_3 [Hz]');
view(-45, 20);
set(gca, 'FontSize', 12);

% -------------------------------------------------------------------------
% FIGURA 5: Espaço 2D Original (Projeções)
% -------------------------------------------------------------------------
fig5 = figure('Name', 'Figura 5 - Espaço de Características 2D', 'Color', 'w', 'Position', [100, 100, 800, 600]);

subplot(2,1,1); hold on; grid on;
for c = 1:9
    idx_s = (Labels_s_mat == c);
    idx_t = (Labels_t == c);
    scatter(Xs_mat(idx_s, 1), Xs_mat(idx_s, 2), 30, mapa_cores(c,:), '^', 'filled');
    scatter(Xt(idx_t, 1), Xt(idx_t, 2), 30, mapa_cores(c,:), 'x', 'LineWidth', 1.2);
end
ylabel('f_2 [Hz]'); xlabel('f_1 [Hz]');
set(gca, 'FontSize', 11);

subplot(2,1,2); hold on; grid on;
for c = 1:9
    idx_s = (Labels_s_mat == c);
    idx_t = (Labels_t == c);
    scatter(Xs_mat(idx_s, 1), Xs_mat(idx_s, 3), 30, mapa_cores(c,:), '^', 'filled');
    scatter(Xt(idx_t, 1), Xt(idx_t, 3), 30, mapa_cores(c,:), 'x', 'LineWidth', 1.2);
end
ylabel('f_3 [Hz]'); xlabel('f_1 [Hz]');
set(gca, 'FontSize', 11);

% -------------------------------------------------------------------------
% FIGURA 7: Espaço Latente após TCA
% -------------------------------------------------------------------------
% Verificação de existência das variáveis latentes para evitar crash
if exist('Zs', 'var') && exist('Zt', 'var')
    fig7 = figure('Name', 'Figura 7 - Espaço Latente (TCA)', 'Color', 'w', 'Position', [150, 150, 900, 600]);

    % Subplot 1
    subplot(2,2,1); hold on; grid on;
    for c = 1:9
        if c <= 3, cor = cores_dano(1,:); elseif c <= 6, cor = cores_dano(2,:); else, cor = cores_dano(3,:); end
        idx_s = (Labels_s_mat == c);
        testes_s = find(idx_s); 
        scatter(testes_s, Zs(idx_s), 25, cor, '^', 'filled');
    end
    ylabel('$\mathcal{F}_1$', 'Interpreter', 'latex'); xlabel('Tests'); title('Latent features of source domain');
    set(gca, 'FontSize', 10);

    % Subplot 2
    subplot(2,2,2); hold on; grid on;
    for c = 1:9
        if c <= 3, cor = cores_dano(1,:); elseif c <= 6, cor = cores_dano(2,:); else, cor = cores_dano(3,:); end
        idx_t = (Labels_t == c);
        testes_t = find(idx_t);
        scatter(testes_t, Zt(idx_t), 25, cor, 'x', 'LineWidth', 1.2);
    end
    ylabel('$\mathcal{F}_1$', 'Interpreter', 'latex'); xlabel('Tests'); title('Latent features of target domain');
    set(gca, 'FontSize', 10);

    % Subplot 3
    subplot(2,2,[3,4]); hold on;
    for lvl = 1:3
        conds = (lvl-1)*3 + 1 : lvl*3;
        idx_s_lvl = ismember(Labels_s_mat, conds);
        idx_t_lvl = ismember(Labels_t, conds);

        % Histograma do Domínio de Origem (Cor sólida e preenchida)
        histogram(Zs(idx_s_lvl), ...
            'FaceColor', cores_dano(lvl,:), ...
            'EdgeColor', cores_dano(lvl,:), ...
            'FaceAlpha', 1);

        % Histograma do Domínio Alvo (Barras vazadas com contorno colorido)
        histogram(Zt(idx_t_lvl), ...
            'DisplayStyle', 'bar', ...
            'FaceColor', 'none', ...
            'EdgeColor', cores_dano(lvl,:), ...
            'LineWidth', 1.2);
    end
    xlabel('$\mathcal{F}_1$', 'Interpreter', 'latex'); 
    ylabel('Amplitude');
    title('Histograms of source (filled) and target (empty) features');
    set(gca, 'FontSize', 11);


else
    warning('As variáveis latentes Zs e Zt não foram encontradas no Workspace. A Figura 7 (TCA) não será gerada.');
end
