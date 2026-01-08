function complex_coeffs = complex_coeff (V_jk_t,N_harmonics)
% this function uses the fft to determine the complex coefficients of the
% fourier series for V(t)
% VALIDATED
% the fft is computed, then centered
    all_coeff = fftshift(fft(V_jk_t)/length(V_jk_t));
    mid_index = 1 + floor(length(V_jk_t)/2);
    % the required N_harmonics coefficients are assigned
    complex_coeffs = all_coeff(mid_index-N_harmonics:mid_index+N_harmonics);

end


