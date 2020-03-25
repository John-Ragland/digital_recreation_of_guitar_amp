% DEMO
% created by John Ragland

clear
clc
close all

%Load sample audio
load('audio.mat');
Fs = 44.1e+3;

%Enter values of Tone Stack and volume
%   vol must be < 1;
l = 0.8;
m = 0.5;
t = 1;
vol = 0.99;

output = tube_pre_eq(l,m,t,vol,audio);

fprintf('Playing Original Audio:\n');
sound(audio,Fs);
pause

fprintf('Playing Output of Algorithm:\n');
sound(output,Fs);