% 02_Extrair_Dados.m
clc; clear;

fprintf('Importando dados numéricos brutos...\n');
% 1. Carrega o arquivo original que tem a variável wn_num (1x11 cell)
load('NatFreq_numerical.mat', 'wn_num');

N_rodadas = 11;
Xs_todas = wn_num;
Labels_todas = cell(1, N_rodadas);

% 2. Criação das labels em formato de célula (compartimentado por rodada)
labels_base = repelem((1:9)', 150); % 1350 amostras (150 de cada uma das 9 condições)

for rodada = 1:N_rodadas
    Labels_todas{rodada} = labels_base;
end

% 3. Recriação do vetor hc_array (que o seu script 03 usa nos prints)
hc_array = 0.0045:-0.0003:0.0015;

% 4. Salva tudo no arquivo 'Dados_Origem.mat'
save('Dados_Origem.mat', 'Xs_todas', 'Labels_todas', 'hc_array');
fprintf('Pronto! Arquivo "Dados_Origem.mat" gerado com sucesso.\n');