function [Zs, Zt, Yt_pseudo] = JDA(Xs, Xt, Ys, m, mu, max_iter, classifier_type)
    % JDA - Joint Distribution Adaptation
    % classifier_type: 'knn' (padrão), 'svm_linear', ou 'svm_rbf'

    % Define 'knn' como padrão caso o usuário não passe o 7º argumento
    if nargin < 7
        classifier_type = 'knn'; 
    end

    ns = size(Xs, 1);
    nt = size(Xt, 1);
    n = ns + nt;
    X = [Xs; Xt];
    
    % 1. Construção do Kernel Gaussiano (RBF)
    D2 = pdist2(X, X, 'euclidean').^2;
    sigma = median(sqrt(D2(:))); 
    K = exp(-D2 ./ (2 * sigma^2));
    
    % 2. Matriz de Centralização (Preservação da variância)
    H = eye(n) - (1/n) * ones(n, n);
    
    % Mapeamento das classes do domínio de origem
    classes = unique(Ys);
    C = length(classes);
    
    % 3. Inicialização: Obtém os primeiros pseudo-rótulos
    Yt_pseudo = gerar_pseudo_rotulos(Xs, Ys, Xt, classifier_type);
    
    % Loop Iterativo do JDA
    for iter = 1:max_iter
        Yt_pseudo_old = Yt_pseudo;
    
        % 4. Matriz de MMD Marginal (M0)
        M0 = zeros(n, n);
        M0(1:ns, 1:ns) = 1 / (ns^2);
        M0(ns+1:end, ns+1:end) = 1 / (nt^2);
        M0(1:ns, ns+1:end) = -1 / (ns * nt);
        M0(ns+1:end, 1:ns) = -1 / (ns * nt);
    
        % 5. Matriz de MMD Condicional (Mc)
        Mc = zeros(n, n);
        for c = 1:C
            classe_atual = classes(c);
            idx_s = find(Ys == classe_atual);
            idx_t = find(Yt_pseudo == classe_atual);
    
            n_sc = length(idx_s);
            n_tc = length(idx_t);
    
            if n_sc > 0 && n_tc > 0
                Mc_temp = zeros(n, n);
                Mc_temp(idx_s, idx_s) = 1 / (n_sc^2);
                Mc_temp(ns + idx_t, ns + idx_t) = 1 / (n_tc^2);
                Mc_temp(idx_s, ns + idx_t) = -1 / (n_sc * n_tc);
                Mc_temp(ns + idx_t, idx_s) = -1 / (n_sc * n_tc);
    
                Mc = Mc + Mc_temp;
            end
        end
    
        % Matriz de Discrepância Conjunta
        M = M0 + Mc;
    
        % 6. Problema Generalizado de Autovalores
        I = eye(n);
        A = K * M * K + mu * I;
        B = K * H * K + 1e-8 * I; % Piso numérico
        
        A = (A + A') / 2;
        B = (B + B') / 2;
        
        [W, D_eig] = eig(A, B);
        
        autovalores = diag(D_eig);
        [~, indices] = sort(autovalores, 'ascend');
        W_reduzido = W(:, indices(1:m));
    
        % 7. Extração do subespaço latente
        Z = K * W_reduzido;
        Zs = Z(1:ns, :);
        Zt = Z(ns+1:end, :);
    
        % 8. Atualização iterativa dos pseudo-rótulos no espaço latente
        Yt_pseudo = gerar_pseudo_rotulos(Zs, Ys, Zt, classifier_type);
        
        % 9. Critério de Parada Antecipada
        mudancas_rotulos = sum(Yt_pseudo ~= Yt_pseudo_old);
        if mudancas_rotulos == 0
            fprintf('JDA convergiu na iteração %d usando %s.\n', iter, classifier_type);
            break;
        end
    end

    % =====================================================================
    % SUBFUNÇÃO: Lógica de Classificação Flexível
    % =====================================================================
    function Y_pred = gerar_pseudo_rotulos(X_treino, Y_treino, X_teste, tipo)
        switch lower(tipo)
            case 'knn'
                % 1-NN
                mdl = fitcknn(X_treino, Y_treino, 'NumNeighbors', 1);
                Y_pred = predict(mdl, X_teste);
                
            case 'svm_linear'
                % SVM Padrão
                mdl = fitcecoc(X_treino, Y_treino);
                Y_pred = predict(mdl, X_teste);
                
            case 'svm_rbf'
                % SVM Gaussiano
                t = templateSVM('KernelFunction', 'gaussian');
                mdl = fitcecoc(X_treino, Y_treino, 'Learners', t);
                Y_pred = predict(mdl, X_teste);
                
            otherwise
                error('Classificador inválido. Escolha: ''knn'', ''svm_linear'' ou ''svm_rbf''.');
        end
    end
end