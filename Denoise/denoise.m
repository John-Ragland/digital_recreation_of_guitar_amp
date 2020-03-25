function [audio_filt] = denoise(audio)
%DENOISE [filt_aud] = denoise(audio) takes the input audio and filter's it
%using the designed filter information contained in Designed_Filter.mat

    load('Designed_Filter2.mat')
    audio_filt = filter(b,a,filter(b2,a2,filter(b3,a3,filter(b4,a4,filter(b5,a5,audio)))));
end

