%
%=====================================================================================
%       Filename:  sincMakePhase.m
% 
%    Description:  make the sinc function, based on given dataset
%        Version:  1.0
%
%         Author:  Kang Min Bae, Hankyeol Moon
%         Email :  smilelabkaist@gmail.com
%   Organization:  Smart and Mobile Systems (Smile) Lab @ KAIST 
%                  https://smile.kaist.ac.kr/
%
%   Copyright (c)  Smart and Mobile Systems (Smile) Lab @ KAIST 
%   This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
% =====================================================================================
%

function resultFunc = sincMakePhase(dataLength, numZeroPads, paddedFreqData, shiftedPoint1, shiftedPoint2)
    fc = 24.125e9;
    lightSpeed = physconst('LightSpeed');
    lambda = lightSpeed/fc;
    
    freqDistance = abs(shiftedPoint1-shiftedPoint2)/2*1000;
    distanceResult= 2*(freqDistance*lightSpeed*(dataLength-1)/(2*1e6*250*1e6));
    phase = 4*pi*distanceResult/lambda;
    phaseDifference = 2*phase;
    
    x = zeros(numZeroPads*dataLength,1);
    x(1:dataLength) = 1;

    Y=fft(x);
    L = length(Y);

    freqResult = fftshift(Y/L);

    [value,zeroIndex]=max(abs(freqResult));

    freq1Index = interp1(paddedFreqData,1:length(paddedFreqData),shiftedPoint1,'nearest');
    freq2Index = interp1(paddedFreqData,1:length(paddedFreqData),shiftedPoint2,'nearest');
    
    shift1 = freq1Index - zeroIndex;
    shift2 = freq2Index - zeroIndex;

    shiftFreqResult1 = [zeros(shift1,1);freqResult(1:end-shift1)];
    shiftFreqResult2 = [zeros(shift2,1);freqResult(1:end-shift2)*exp(1j*phaseDifference)];


    shiftFreqResult = shiftFreqResult1+shiftFreqResult2;
    
    resultFunc = shiftFreqResult;

end