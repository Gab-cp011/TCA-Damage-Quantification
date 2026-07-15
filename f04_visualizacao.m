% =========================================================================
% MÓDULO DE VISUALIZAÇÃO - REPRODUÇÃO DAS FIGURAS 4, 5 E 7
% =========================================================================

% Definição da paleta de cores oficial das 9 condições (RGB)
% 1:Azul, 2:Vermelho, 3:Preto, 4:Magenta, 5:Verde, 6:Cinza, 7:Laranja, 8:Ciano, 9:Amarelo
mapa_cores = [0 0 1; 1 0 0; 0 0 0; 1 0 1; 0 0.6 0; 0.5 0.5 0.5; 0.8 0.4 0; 0 0.7 1; 1 0.8 0];

% Cores agrupadas por Nível de Dano (Azul: Saudável, Vermelho: Dano 1, Verde: Dano 2)
cores_dano = [0 0 1; 1 0 0; 0 0.6 0]; 

%% FIGURA 4: Espaço 3D Original (Frequências Naturais)
fig4 = figure('Name', 'Figura 4 - Espaço de Características 3D', 'Color', 'w');
hold on; grid on;
for c = 1:9
    idx_s = (Labels_s == c);
    idx_t = (Labels_t == c);

    % Origem (Triângulos preenchidos) e Alvo (Cruzes)
    scatter3(Xs(idx_s, 1), Xs(idx_s, 2), Xs(idx_s, 3), 40, mapa_cores(c,:), '^', 'filled');
    scatter3(Xt(idx_t, 1), Xt(idx_t, 2), Xt(idx_t, 3), 40, mapa_cores(c,:), 'x', 'LineWidth', 1.2);
end
xlabel('f_1 [Hz]'); ylabel('f_2 [Hz]'); zlabel('f_3 [Hz]');
view(-45, 20); % Ajuste de câmera para perspectiva similar ao artigo
set(gca, 'FontSize', 12);

%% FIGURA 5: Espaço 2D Original (Projeções)
fig5 = figure('Name', 'Figura 5 - Espaço de Características 2D', 'Color', 'w', 'Position', [100, 100, 800, 600]);

% Subplot Superior: f1 vs f2
subplot(2,1,1); hold on; grid on;
for c = 1:9
    idx_s = (Labels_s == c);
    idx_t = (Labels_t == c);
    scatter(Xs(idx_s, 1), Xs(idx_s, 2), 30, mapa_cores(c,:), '^', 'filled');
    scatter(Xt(idx_t, 1), Xt(idx_t, 2), 30, mapa_cores(c,:), 'x', 'LineWidth', 1.2);
end
ylabel('f_2 [Hz]'); xlabel('f_1 [Hz]');
set(gca, 'FontSize', 11);

% Subplot Inferior: f1 vs f3
subplot(2,1,2); hold on; grid on;
for c = 1:9
    idx_s = (Labels_s == c);
    idx_t = (Labels_t == c);
    scatter(Xs(idx_s, 1), Xs(idx_s, 3), 30, mapa_cores(c,:), '^', 'filled');
    scatter(Xt(idx_t, 1), Xt(idx_t, 3), 30, mapa_cores(c,:), 'x', 'LineWidth', 1.2);
end
ylabel('f_3 [Hz]'); xlabel('f_1 [Hz]');
set(gca, 'FontSize', 11);

%% FIGURA 7: Espaço Latente após TCA (1ª Componente Principal)
fig7 = figure('Name', 'Figura 7 - Espaço Latente (TCA)', 'Color', 'w', 'Position', [150, 150, 900, 600]);

% Subplot 1: Dispersão da Origem (Esquerda)
subplot(2,2,1); hold on; grid on;
for c = 1:9
    if c <= 3, cor = cores_dano(1,:); elseif c <= 6, cor = cores_dano(2,:); else, cor = cores_dano(3,:); end
    idx_s = (Labels_s == c);
    testes_s = find(idx_s); 
    scatter(testes_s, Zs(idx_s), 25, cor, '^', 'filled');
end
% CORREÇÃO AQUI
ylabel('$\mathcal{F}_1$', 'Interpreter', 'latex'); xlabel('Tests'); title('Latent features of source domain');
set(gca, 'FontSize', 10);

% Subplot 2: Dispersão do Alvo (Direita)
subplot(2,2,2); hold on; grid on;
for c = 1:9
    if c <= 3, cor = cores_dano(1,:); elseif c <= 6, cor = cores_dano(2,:); else, cor = cores_dano(3,:); end
    idx_t = (Labels_t == c);
    testes_t = find(idx_t);
    scatter(testes_t, Zt(idx_t), 25, cor, 'x', 'LineWidth', 1.2);
end
% CORREÇÃO AQUI
ylabel('$\mathcal{F}_1$', 'Interpreter', 'latex'); xlabel('Tests'); title('Latent features of target domain');
set(gca, 'FontSize', 10);

% Subplot 3: Histogramas Sobrepostos (Base)
subplot(2,2,[3,4]); hold on;
for lvl = 1:3
    conds = (lvl-1)*3 + 1 : lvl*3;
    idx_s_lvl = ismember(Labels_s, conds);
    idx_t_lvl = ismember(Labels_t, conds);
    
    histogram(Zs(idx_s_lvl), 'BinWidth', 10, 'FaceColor', cores_dano(lvl,:), 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    histogram(Zt(idx_t_lvl), 'BinWidth', 10, 'DisplayStyle', 'stairs', 'EdgeColor', cores_dano(lvl,:), 'LineWidth', 1.5);
end
% CORREÇÃO AQUI
xlabel('$\mathcal{F}_1$', 'Interpreter', 'latex'); ylabel('Amplitude');
title('Histograms of source (filled) and target (empty) features');
set(gca, 'FontSize', 11);