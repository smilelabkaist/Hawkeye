function resultFunc = sincMake_phase(sample_data_len, pad_num, total_padded_freq, maxfreq1, maxfreq2)
    fc = 24.125e9;
    c_light = physconst('LightSpeed');
    lambdac = c_light/fc;
    
    f_dist = abs(maxfreq1-maxfreq2)/2*1000;
    dist= 2*(f_dist*c_light*(sample_data_len-1)/(2*1e6*250*1e6));
    phase = 4*pi*dist/lambdac;
    phase_gap = 2*phase;
    
    x=ones(sample_data_len,1);
    x = [x ; zeros(pad_num*sample_data_len,1)];

    Y=fft(x);
    L = length(Y);

    Result = fftshift(Y/L);

    [m,zeroindex]=max(abs(Result));

    r1 = Result;
    r2 = Result*exp(1j*phase_gap);

    freq1index = interp1(total_padded_freq,1:length(total_padded_freq),maxfreq1,'nearest');
    freq2index = interp1(total_padded_freq,1:length(total_padded_freq),maxfreq2,'nearest');
    
    shift1 = freq1index - zeroindex;
    shift2 = freq2index - zeroindex;

    shift_fft1 = [zeros(shift1,1);r1(1:end-shift1)];
    shift_fft2 = [zeros(shift2,1);r2(1:end-shift2)];


    shift_fft = shift_fft1+shift_fft2;
    
    resultFunc = shift_fft;

end