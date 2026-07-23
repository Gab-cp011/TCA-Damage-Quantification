function [Zs, Zt] = TCA(Xs, Xt, m, mu)
    % Obtém o número de amostras para os domínios fonte e alvo.
    ns = size(Xs, 1);
    nt = size(Xt, 1);
    n = ns + nt;
    
    % Concatena os domínios para construir as matrizes globais
    X = [Xs; Xt];
    
    % 1. Construção da Matriz de Kernel Linear
    % Armazena em cada elemento o produto interno entre duas amostras,
    % ilustrando a "similaridade"
    K = X * X'; 
    
    % 2. Construção da Matriz L (Maximum Mean Discrepancy)
    L = zeros(n, n);
    L(1:ns, 1:ns) = 1 / (ns * ns);
    L(ns+1:end, ns+1:end) = 1 / (nt * nt);
    L(1:ns, ns+1:end) = -1 / (ns * nt);
    L(ns+1:end, 1:ns) = -1 / (ns * nt);
    
    % 3. Construção da Matriz de Centralização H
    H = eye(n) - (1/n) * ones(n, n);
    
    % 4. Formulação do Problema Generalizado de Autovalores
    % (K*L*K + mu*I) * W = lambda * (K*H*K) * W
    I = eye(n);
    A = K * L * K + mu * I;
    B = K * H * K;
    
    % Adiciona um valor infinitesimal à diagonal de B para garantir estabilidade numérica
    B = B + 1e-8 * eye(n);
    
    % Resolve extraindo os autovalores e autovetores
    [W, D] = eig(A, B);
    autovalores = diag(D);
    
    % Ordena os autovalores em ordem crescente para minimizar a discrepância
    [~, indices] = sort(autovalores, 'ascend');
    W_ordenado = W(:, indices);
    
    % 5. Seleciona as 'm' primeiras dimensões para definir o espaço latente
    W_reduzido = W_ordenado(:, 1:m);
    
    % 6. Projeta os dados originais no novo espaço
    Z = K * W_reduzido;

    % Separa os dados de volta para Origem e Destino
    Zs = Z(1:ns, :);
    Zt = Z(ns+1:end, :);
end
