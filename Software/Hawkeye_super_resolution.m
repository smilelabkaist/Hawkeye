
%
%=====================================================================================
%       Filename:  Hawkeye_super_resolution.m
% 
%    Description:  load the EVAL-TINYRAD raw data and apply Hawkeye super resolution.
%        Version:  1.0
%
%         Author:  Kang Min Bae, Hankyeol Moon
%         Email :  {bkm2259, moonkyul1}@kaist.ac.kr
%   Organization:  Smart and Mobile Systems (Smile) Lab @ KAIST 
%                  https://smile.kaist.ac.kr/
%
%   Copyright (c)  Smart and Mobile Systems (Smile) Lab @ KAIST 
% =====================================================================================
%   - The code does not apply any form of beamforming. It is designed to use a single Rx and Tx antenna.
%   - The code works best when chirp interval was integer multiple of 1/f_mod at radar.
%   - The code assumes that the f_mod is always larger than f_range. We will not provide methods to do otherwise.
% =====================================================================================
%

%% Specify super-resolution options.

numTags = 1; % Number of tags to detect. If numTags = n, the strongest n signals in tag bin will be analyzed.
plotProcess = 0; % Plot the super-resolution process. 0 for no plot, 1 for plot in every step, 2 for just the result.
numZeroPads = 1024; % Number of zeros to pad for Hawkeye super-resolution. Refer to the paper for more information.
maxRange = 10; % Default is inf. For better performance, apply the max range limit.
receiverAntenna = 1; % Select the receiver antenna to use. Options are 1 to 4.
f_modMax = 149; % Input maximum tag modulation frequency, in kHz.
f_modMin = 151; % Input minimum tag modulation frequency, in kHz.
mobileTagLocalizationInterval = 40; % Mobile tag localization interval in #chirps.
applyAdvancedSR = 0; % 0:OFF, 1:ON. Apply a slower, but more accurate super resolution algorithm. Typically recommended when range < 3 meter.
mobileTagLocalization = 0; % 0:OFF, 1:ON. Turn on if tag is mobile.

processingType = (mobileTagLocalization-1)*4+applyAdvancedSR+1;

%% Load the raw data

filename = ''; % Input the raw data file name.

disp(filename);
loadedData = load(filename);
Data=loadedData.Data;

[rows, columns] = size(Data);
N = rows - 2;
NrFrms = columns/4;
Np = 1;
fs = 1e6;

maxFrequencyRange  = maxRange/100*30;
maxFrequencyRange = maxFrequencyRange*(8192/(N+1));

totalTime = Data(1,end);

chn = zeros(NrFrms*fix((totalTime*fs+N)/NrFrms)+NrFrms,1);
chnCount = 1;

for Idx = 1:NrFrms
        chn(chnCount:chnCount+N-1) = chn(chnCount:chnCount+N-1) + Data(3:end,(receiverAntenna)+4*(Idx-1));
    chnCount = chnCount + N;
end

chnTemp = zeros(chnCount-1,1);
chnTemp(1:chnCount-1)=chn(1:chnCount-1);
chn=chnTemp;

fftResult = fft(chn);
length = double(N);
fftResultAbs = abs(fftResult/length);
fftResultFinal = fftshift(fftResultAbs);

frequencyIndex=((1:length(fftResultAbs))-round(length(fftResultAbs)/2)).';
frequencyHz=frequencyIndex/length(frequencyIndex)*fs;
frequencykHz=frequencyHz/1000;

if plotProcess==1
    figure
    plot(frequencykHz, fftResultFinal)
end

fftResultFinalComplex = fftshift(fftResult/length);

%% Run clutter rejection.

centerClutterHigh = interp1(frequencyHz,1:length(frequencyHz),100,'nearest');
centerClutterLow = interp1(frequencyHz,1:length(frequencyHz),-100,'nearest');

peakBWDown = 8*8192/(N+1);

[value, impulse] = max(fftResultFinal(centerClutterLow:centerClutterHigh));

impulse=impulse+centerClutterLow-1;
for temp_index = impulse:-NrFrms:1
    fftResultFinal(max(round(temp_index)-peakBWDown,1) : round(temp_index)+peakBWDown) = 0;
    fftResultFinalComplex(max(round(temp_index)-peakBWDown,1) : round(temp_index)+peakBWDown) = 0;
end
for temp_index = impulse+NrFrms:NrFrms:length(fftResultFinal)
    fftResultFinal(round(temp_index)-peakBWDown : min(round(temp_index)+peakBWDown,length(fftResultFinal))) = 0;
    fftResultFinalComplex(round(temp_index)-peakBWDown : min(round(temp_index)+peakBWDown,length(fftResultFinalComplex))) = 0;
end

if plotProcess == 1
    figure()
    plot(frequencykHz, fftResultFinal) 
end

%% Check the noise floor

noiseFreqLow = interp1(frequencyHz,1:length(frequencyHz),499.95e3,'nearest');
noiseFreqHigh = interp1(frequencyHz,1:length(frequencyHz),499.99e3,'nearest');
noiseLevel = max(fftResultFinal(noiseFreqLow:noiseFreqHigh));

%% Get max peak (i.e., tag), erase the peaks, find next tag.

tagData = zeros(length(fftResultFinal),numTags);

temp_fftResultFinal = fftResultFinal;

tagDataComplex = zeros(length(fftResultFinal),numTags);
tagDataComplexMobile = zeros(length(fftResultFinal),numTags);

freqLowerBound = interp1(frequencyHz,1:length(frequencyHz),(f_modMin-maxFrequencyRange-10)*1000,'nearest');
freqHigherBound = interp1(frequencyHz,1:length(frequencyHz),(f_modMax+maxFrequencyRange+10)*1000,'nearest');

tagIndexes = zeros(3,numTags);
for temp_tagIndex = 1:numTags
    peakBWUp = 8*8192/(N+1); 
    peakBWDown = 8*8192/(N+1); 
    temp_tagPeakCount = 0;
    [value, impulse] = max(temp_fftResultFinal(freqLowerBound:freqHigherBound));
    impulse = impulse + freqLowerBound - 1;

    while (temp_fftResultFinal(impulse+peakBWUp+1) * 2 > noiseLevel) || (temp_fftResultFinal(impulse) < temp_fftResultFinal(impulse+peakBWUp+1) * 2)
        peakBWUp = peakBWUp + 1;
    end
    
    while (temp_fftResultFinal(impulse-down_guard_time-1) * 2 > noiseLevel) || (temp_fftResultFinal(impulse) < temp_fftResultFinal(impulse-down_guard_time-1) * 2)
        down_guard_time = down_guard_time + 1;
    end

    if (processingType == 4)||(processingType == 5)
        peakBWUp = interp1(frequencyHz,1:length(frequencyHz),70,'nearest')-interp1(frequencyHz,1:length(frequencyHz),0,'nearest');
        down_guard_time = peakBWUp;
    end
    
    for temp_index = impulse:-NrFrms:1
        tagData(round(temp_index),temp_tagIndex) = temp_fftResultFinal(temp_index);
        tagDataComplex(round(temp_index),temp_tagIndex) = fftResultFinalComplex(temp_index);
        tagDataComplexMobile(max(round(temp_index)-down_guard_time,1) : round(temp_index)+peakBWUp,temp_tagIndex) = fftResultFinalComplex(max(round(temp_index)-down_guard_time,1) : round(temp_index)+peakBWUp);
        temp_fftResultFinal(max(round(temp_index)-down_guard_time,1) : round(temp_index)+peakBWUp) = 0;
        temp_tagPeakCount = temp_tagPeakCount+1;
    end

    tagPeakCountDown = temp_index;
    for temp_index = impulse+NrFrms:NrFrms:length(fftResultFinal)
        tagData(round(temp_index),temp_tagIndex) = temp_fftResultFinal(temp_index);
        tagDataComplex(round(temp_index),temp_tagIndex) = fftResultFinalComplex(temp_index);
        tagDataComplexMobile(round(temp_index)-down_guard_time : min(round(temp_index)+peakBWUp,length(fftResultFinal)),temp_tagIndex) = fftResultFinalComplex(round(temp_index)-down_guard_time : min(round(temp_index)+peakBWUp,length(fftResultFinal)));
        temp_fftResultFinal(round(temp_index)-down_guard_time : min(round(temp_index)+peakBWUp,length(fftResultFinal))) = 0;
        temp_tagPeakCount = temp_tagPeakCount+1;
    end
    tagPeakCountUp = temp_index;
    tagIndexes(:,temp_tagIndex) = [temp_tagPeakCount;tagPeakCountDown;tagPeakCountUp];

end

if (processingType == 4)||(processingType == 5)
    NrFrmsInkHz = frequencykHz(1+NrFrms)-frequencykHz(1);

    for mobileProcessingIndex=1:floor(NrFrms/mobileTagLocalizationInterval)+1
        chirpIndex = 1+(mobileProcessingIndex-1)*mobileTagLocalizationInterval;
        
        total_padded_fft_data = zeros(N*(numZeroPads),numTags);
        total_padded_freq = zeros(N*(numZeroPads),numTags);
        for temp_tagIndex = 1:numTags
            ifft_data = ifft(tagDataComplexMobile(:,temp_tagIndex));
            
            zero_p_current = zeros(1,N*(numZeroPads));
            zero_p_current(1:N) = ifft_data((chirpIndex-1)*N+1:chirpIndex*N);

            padded_fft = fft(zero_p_current);
    
            pad_frequencyIndex = ((1:length(padded_fft))-round(length(padded_fft)/2)).';
            pad_frequencyHz = pad_frequencyIndex/length(pad_frequencyIndex)*fs;
            pad_frequencykHz = pad_frequencyHz/1000;
            
            total_padded_fft_data(:,temp_tagIndex) = abs(padded_fft);
            total_padded_freq(:,temp_tagIndex) = pad_frequencykHz;
    
            %% find f_dist
    
            padded_fft_cut = zeros(length(padded_fft),1);
    
            freqLowerBound = interp1(pad_frequencykHz,1:length(pad_frequencykHz),(150-maxFrequencyRange-10),'nearest');
            freqHigherBound = interp1(pad_frequencykHz,1:length(pad_frequencykHz),(150+maxFrequencyRange+10),'nearest');
            padded_fft_cut(freqLowerBound:freqHigherBound) = padded_fft(freqLowerBound:freqHigherBound);
    
            [a,b,c] = findpeaks(abs(padded_fft_cut)./1e4,pad_frequencykHz);
    
            maxvalues = zeros(1,2);
            maxindex = zeros(1,2);
            for temp_index = 1:2
                [value, impulse] = max(a);
                a(impulse) = 0;
                maxvalues(1,temp_index) = value;
                maxindex(1,temp_index) = impulse;
            end
    
            maxindex2 = sort(maxindex);
    
    
            if (b(maxindex2(2))-b(maxindex2(1))) < 1.5*NrFrmsInkHz
                [nvalue1, impulse] = max(maxvalues);
                maxfreq1 = b(maxindex(impulse));
                maxfreq2 = maxfreq1;
            else
                maxfreq1 = b(maxindex2(1));
                maxfreq2 = b(maxindex2(2));
            end
    
    
            f_dist = abs(maxfreq1-maxfreq2)/2*1000;
    
            %% calculate distance
    
            c_light = physconst('LightSpeed');
    
            distance_result = f_dist*c_light*(N-1)/(2*1e6*250*1e6)

            if plotProcess > 0
                figure()
                plot(total_padded_freq(:,temp_tagIndex),total_padded_fft_data(:,temp_tagIndex)); 
                title(filename)
                %xlim([198,202]);
            end
            
    
            %% Save Peak Amplitude
    
            ref_peak_sample = min(interp1(pad_frequencykHz,1:length(pad_frequencykHz), maxfreq1,'nearest'), interp1(pad_frequencykHz,1:length(pad_frequencykHz), maxfreq2,'nearest'));
            far_peak_sample = max(interp1(pad_frequencykHz,1:length(pad_frequencykHz), maxfreq1,'nearest'), interp1(pad_frequencykHz,1:length(pad_frequencykHz), maxfreq2,'nearest'));
            f_peak_mean = (total_padded_fft_data(ref_peak_sample,temp_tagIndex)+total_padded_fft_data(far_peak_sample,temp_tagIndex))/2;

            %% Do RMS compare.
    
            if processingType == 5
    
                if maxfreq2 < maxfreq1
                   [maxfreq2, maxfreq1] = deal(maxfreq1,maxfreq2); 
                end
    
    
                freqUnit = total_padded_freq(2,temp_tagIndex)-total_padded_freq(1,temp_tagIndex);
                shift1 = maxfreq1;
                shift2 = maxfreq2;
    
                freq1index = find(total_padded_freq(:,temp_tagIndex)==maxfreq1);
                freq2index = find(total_padded_freq(:,temp_tagIndex)==maxfreq2);
                sweepUnit = 2;
                sweepEnd= 50;
                sincSweepRange=-gap_end*freqUnit:sweepUnit*freqUnit:gap_end*freqUnit;
                
                min_cnt = 0;
                min_rms_sum = inf;
                sincs_final = 0;
    
                cnt = 0;
                
                for temp_index = sincSweepRange
                    cnt = cnt+1;
                    shift1 = maxfreq1+temp_index;
                    shift2 = maxfreq2-temp_index;
                    sincs = abs(sincMake_phase(N,numZeroPads, total_padded_freq(:,temp_tagIndex),shift1, shift2));
                    sincs = (f_peak_mean/max(sincs))*sincs;
    
                    rmswindow = [freq1index-numZeroPads*6:freq2index+numZeroPads*6];
    
                    sinc_sub = total_padded_fft_data(rmswindow,temp_tagIndex)-sincs(rmswindow);
                    
                    rms_sum = rms(sinc_sub);
    
                    if rms_sum < min_rms_sum
                        min_cnt = cnt;
                        min_rms_sum = rms_sum;
                        sincs_final = sincs;
                    end
                end

                new_max_freq1 = maxfreq1+sincSweepRange(min_cnt);
                new_max_freq2 = maxfreq2-sincSweepRange(min_cnt);
    
                f_dist = abs(new_max_freq1-new_max_freq2)/2*1000;
    
                advanced_distance_result = f_dist*c_light*(N-1)/(2*1e6*250*1e6)
    
    
                if plotProcess > 0
                    figure()
                    plot(total_padded_freq(:,temp_tagIndex),total_padded_fft_data(:,temp_tagIndex));
                    hold on
                    plot(total_padded_freq(:,temp_tagIndex),sincs_final)
                    hold off
                    drawnow;
     
                end
            end
        end
    end
end

if (processingType == 1) || (processingType == 2)
    for temp_tagIndex = 1:numTags
        maxvalues = zeros(1,3);
        maxindex = zeros(1,3);
        temp_fftResultFinal = tagData(:,temp_tagIndex);
        for temp_index = 1:3
            [value, impulse] = max(temp_fftResultFinal(freqLowerBound:freqHigherBound));
            impulse = impulse + freqLowerBound - 1;
            temp_fftResultFinal(impulse-peakBWDown : impulse+peakBWDown) = 0;
            maxvalues(1,temp_index) = value;
            maxindex(1,temp_index) = impulse;
        end


    end


    %% New data for padding
    sampled_data = zeros(tagIndexes(1,numTags),numTags);
    sampled_freq = zeros(tagIndexes(1,numTags),numTags);
    for temp_tagIndex = 1:numTags
        if tagIndexes(1,temp_tagIndex) ~= round((tagIndexes(3,temp_tagIndex)-tagIndexes(2,temp_tagIndex))/NrFrms)+1
            disp("Somethings wrong - padding check")
        end
        cnt = 0;
        for j = tagIndexes(2,temp_tagIndex):NrFrms:tagIndexes(3,temp_tagIndex)
            cnt = cnt+1;
            sampled_data(cnt,temp_tagIndex)=tagDataComplex(round(j),temp_tagIndex);
            sampled_freq(cnt,temp_tagIndex) = frequencykHz(round(j));
        end
    end

    total_padded_fft_data = zeros((numZeroPads)*length(sampled_data),numTags);
    total_padded_freq = zeros((numZeroPads)*length(sampled_data),numTags);
    for temp_tagIndex = 1:numTags
        ifft_data = ifft(sampled_data(:,temp_tagIndex));

        zero_p = zeros(1,N*(numZeroPads));
        zero_p(1:N) = ifft_data;

        padded_fft = fft(zero_p);


        pad_frequencyIndex=((1:length(padded_fft))-round(length(padded_fft)/2)).';
        pad_frequencyHz=pad_frequencyIndex/length(pad_frequencyIndex)*fs;
        pad_frequencykHz=pad_frequencyHz/1000;
        

        offsetToSampleFreq = pad_frequencykHz(1)-sampled_freq(1,temp_tagIndex);
        pad_frequencykHz = pad_frequencykHz-offsetToSampleFreq;
        
        [max_impulse_value, max_impulse_index] = max(tagData(:,temp_tagIndex));
        max_impulse_freq = frequencykHz(max_impulse_index);
        candidates = find(round(abs(padded_fft),5)==round(max_impulse_value,5));
        %max 2 candidtates
        if length(candidates) == 2
           if abs(pad_frequencykHz(candidates(1)) - max_impulse_freq) > abs(pad_frequencykHz(candidates(2)) - max_impulse_freq)
              candidates(1) = candidates(2);
           end
        end

        offsetToKHzFreq = pad_frequencykHz(candidates(1)) - max_impulse_freq;
        pad_frequencykHz = pad_frequencykHz-offsetToKHzFreq;

        total_padded_fft_data(:,temp_tagIndex) = abs(padded_fft);
        total_padded_freq(:,temp_tagIndex) = pad_frequencykHz;

        %% find f_dist

        padded_fft_cut = zeros(length(padded_fft),1);

        freqLowerBound = interp1(pad_frequencykHz,1:length(pad_frequencykHz),(150-maxFrequencyRange-10),'nearest');
        freqHigherBound = interp1(pad_frequencykHz,1:length(pad_frequencykHz),(150+maxFrequencyRange+10),'nearest');
        padded_fft_cut(freqLowerBound:freqHigherBound) = padded_fft(freqLowerBound:freqHigherBound);

        [a,b,c] = findpeaks(abs(padded_fft_cut)./1e4,pad_frequencykHz);

        maxvalues = zeros(1,2);
        maxindex = zeros(1,2);
        for temp_index = 1:2
            [value, impulse] = max(a);
            a(impulse) = 0;
            maxvalues(1,temp_index) = value;
            maxindex(1,temp_index) = impulse;
        end

        maxindex2 = sort(maxindex);

        NrFrmsInkHz = frequencykHz(1+NrFrms)-frequencykHz(1);


        if (b(maxindex2(2))-b(maxindex2(1))) < 1.5*NrFrmsInkHz
            [nvalue1, impulse] = max(maxvalues);
            maxfreq1 = b(maxindex(impulse));
            maxfreq2 = maxfreq1;
        else
            maxfreq1 = b(maxindex2(1));
            maxfreq2 = b(maxindex2(2));
        end

        
        f_dist = abs(maxfreq1-maxfreq2)/2*1000;

        %% calculate distance

        c_light = physconst('LightSpeed');

        distance_result = f_dist*c_light*(tagIndexes(1,numTags)-1)/(2*1e6*250*1e6)

        %% Save Peak Amplitude

        ref_peak_sample = min(interp1(pad_frequencykHz,1:length(pad_frequencykHz), maxfreq1,'nearest'), interp1(pad_frequencykHz,1:length(pad_frequencykHz), maxfreq2,'nearest'));
        far_peak_sample = max(interp1(pad_frequencykHz,1:length(pad_frequencykHz), maxfreq1,'nearest'), interp1(pad_frequencykHz,1:length(pad_frequencykHz), maxfreq2,'nearest'));
        f_peak_mean = (total_padded_fft_data(ref_peak_sample,temp_tagIndex)+total_padded_fft_data(far_peak_sample,temp_tagIndex))/2;
            
        if processingType == 1

            if plotProcess > 0
                figure()
                plot(frequencykHz,tagData(:,temp_tagIndex),'k','LineWidth',1)
                hold on
                plot(total_padded_freq(:,temp_tagIndex),total_padded_fft_data(:,temp_tagIndex)); 
                hold off
                title(num2str((maxfreq1+maxfreq2)/2));
                %xlim([198,202]);
            end
            
        end
        %% Do RMS compare.

        if processingType == 2

            if maxfreq2 < maxfreq1
               [maxfreq2, maxfreq1] = deal(maxfreq1,maxfreq2); 
            end


            offsetToSampleFreq = total_padded_freq(2,temp_tagIndex)-total_padded_freq(1,temp_tagIndex);
            shift1 = maxfreq1;
            shift2 = maxfreq2;

            freq1index = find(total_padded_freq(:,temp_tagIndex)==maxfreq1);
            freq2index = find(total_padded_freq(:,temp_tagIndex)==maxfreq2);
            sweepUnit = 2;
            sweepEnd= 50;
            sincSweepRange=-gap_end*offsetToSampleFreq:sweepUnit*offsetToSampleFreq:gap_end*offsetToSampleFreq;
            
            min_cnt = 0;
            min_rms_sum = inf;
            sincs_final = 0;

            cnt = 0;
            
            for temp_index = sincSweepRange
                cnt = cnt+1;
                shift1 = maxfreq1+temp_index;
                shift2 = maxfreq2-temp_index;
                sincs = abs(sincMake_phase(length(sampled_data),numZeroPads, total_padded_freq(:,temp_tagIndex),shift1, shift2));
                sincs = (f_peak_mean/max(sincs))*sincs;

                rmswindow = [freq1index-numZeroPads*6:freq2index+numZeroPads*6];

                sinc_sub = total_padded_fft_data(rmswindow,temp_tagIndex)-sincs(rmswindow);
                
                rms_sum = rms(sinc_sub);

                if rms_sum < min_rms_sum
                    min_cnt = cnt;
                    min_rms_sum = rms_sum;
                    sincs_final = sincs;
                end
            end
            new_max_freq1 = maxfreq1+sincSweepRange(min_cnt);
            new_max_freq2 = maxfreq2-sincSweepRange(min_cnt);

            f_dist = abs(new_max_freq1-new_max_freq2)/2*1000;

            advanced_distance_result = f_dist*c_light*(tagIndexes(1,numTags)-1)/(2*1e6*250*1e6)


            if plotProcess > 0
                figure()
                plot(frequencykHz,tagData(:,temp_tagIndex),'k','LineWidth',1)
                hold on
                plot(total_padded_freq(:,temp_tagIndex),total_padded_fft_data(:,temp_tagIndex));
                plot(total_padded_freq(:,temp_tagIndex),sincs_final)
                hold off
                title(num2str((maxfreq1+maxfreq2)/2))
                xlim([120,180]);
                drawnow;
 
            end
        end
    end
end
