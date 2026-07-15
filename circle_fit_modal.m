function fn = circle_fit_modal(f, H_complex, peak_f)
    % Restrição geométrica de vizinhança
    idx = f > (peak_f - 1.5) & f < (peak_f + 1.5);
    f_band = f(idx);
    H_band = H_complex(idx);
    
    x = real(H_band);
    y = imag(H_band);
    
    % Formulação matricial dos Mínimos Quadrados
    Mat = [x, y, ones(size(x))];
    rhs = x.^2 + y.^2;
    params = Mat \ rhs;
    
    % Coordenadas geométricas do centro
    xc = params(1) / 2;
    yc = params(2) / 2;
    
    % Translação paramétrica
    x_shift = x - xc;
    y_shift = y - yc;
    angulos = unwrap(atan2(y_shift, x_shift));
    
    % Estabilização polinomial para isolar a ressonância
    angulos_smooth = smoothdata(angulos, 'sgolay', 5);
    
    % Derivada analítica estabilizada
    dTheta = abs(diff(angulos_smooth));
    [~, max_idx] = max(dTheta);
    
    % Correção do deslocamento espacial do vetor diferencial
    fn = f_band(max_idx + 1); 
end