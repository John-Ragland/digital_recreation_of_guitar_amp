function [output, output2] = tube_pre_eq(l, m, t, vol, audio)
%TUBE_PRE_EQ        Simulation of Fender Blues Jr.
%   MATLAB Simulation of Fender Blues Jr EQ System and Tube Pre Amp
%
%   [output] = Tube_Pre_Eq(l,m,t,vol,audio) returns the output of a
%   simulated tube pre amp with an EQ RC array termed 'tone stack' where
%   
%   where audio is the input audio waveform (MATLAB Variable)
%   l, m, t, vol are potentiometer posisitions for Bass, Mid, Treble and
%       Volume respectivly (between 0 and 1)
%
%       Volume may not be equal to 0 or 1
%
% John Ragland
% Auburn University
% SOURCES   1. Discretiation of the '59 Fender Bassman Tone Stack - D. T.
%               Yeh
%           2. Automated Physical Modeling of Nonlinear Audio Circuits For
%               Real-Time Audio Effects Parts I and II. - D, T. Yeh


dbstop if error
%Sets Audio to Voltage Level of Humbucker Pickup for Electric Guitar
%   Divide by 2 - humbucker
%   Divide by 4 - single coil
audio = audio/4;

%Add Gain of First Stage
audio = audio*51.2;


Fs = 44100;
Ts = 1/Fs;
c = 2/Ts;

R1 = 250e+3;
R2 = 250e+3;
R3 = 25e+3;
R4 = 100e+3;
C1 = 0.250e-9;
C2 = 22e-9;
C3 = 22e-9;

b1 = t*C1*R1 + m*C3*R3+l*(C1*R2+C2*R2)+(C1*R3 + C2*R3);

b2 = t*(C1*C2*R1*R4+C1*C3*R1*R4)-m^2*(C1*C3*R3^2+C2*C3*R3^2) ...
    + m*(C1*C3*R1*R3+C1*C3*R3^2+C2*C3*R3^2)...
    + l*(C1*C2*R1*R2 + C1*C2*R2*R4+C1*C3*R2*R4)...
    + l*m*(C1*C3*R2*R3+C2*C3*R2*R3)...
    + (C1*C2*R1*R3+C1*C2*R3*R4+C1*C3*R3*R4);

b3 = l*m*(C1*C2*C3*R1*R2*R3+C1*C2*C3*R2*R3*R4)...
    -m^2*(C1*C2*C3*R1*R3^2+C1*C2*C3*R3^2*R4)...
    +m*(C1*C2*C3*R1*R3^2+C1*C2*C3*R3^2*R4)...
    +t*C1*C2*C3*R1*R3*R4-t*m*C1*C2*C3*R1*R3*R4...
    +t*l*C1*C2*C3*R1*R2*R4;

a0 = 1;

a1 = (C1*R1+C1*R3+C2*R3+C2*R4+C3*R4)+m*C3*R3+l*(C1*R2+C2*R2);

a2 = m*(C1*C3*R1*R3-C2*C3*R3*R4+C1*C3*R3^2 ...
    + C2*C3*R3^2)+l*m*(C1*C3*R2*R3+C2*C3*R2*R3) ...
    - m^2*(C1*C3*R3^2+C2*C3*R3^2)+l*(C1*C2*R2*R4 ...
    + C1*C2*R1*R2+C1*C3*R2*R4+C2*C3*R2*R4)...
    +(C1*C2*R1*R4+C1*C3*R1*R4+C1*C2*R3*R4...
    +C1*C2*R1*R3+C1*C3*R3*R4+C2*C3*R3*R4);

a3 = l*m*(C1*C2*C3*R1*R2*R3+C1*C2*C3*R2*R3*R4)...
    -m^2*(C1*C2*C3*R1*R3^2+C1*C2*C3*R3^2*R4)...
    +m*(C1*C2*C3*R3^2*R4+C1*C2*C3*R1*R3^2 ...
    -C1*C2*C3*R1*R3*R4) + l*C1*C2*C3*R1*R2*R4 ...
    + C1*C2*C3*R1*R3*R4;

b = [b3 b2 b1 0];
a = [a3 a2 a1 a0];
ws = (0:(Fs/2)/9999:Fs/2)*2*pi;
f = ws/(2*pi);
[Hs,ws] = freqs(b,a,ws);


%% Creating of Z Domain Filter with Bilinear Transform

B0 = -b1*c - b2*c^2 - b3*c^3;
B1 = -b1*c + b2*c^2 + 3*b3*c^3;
B2 = b1*c + b2*c^2 - 3*b3*c^3;
B3 = b1*c - b2*c^2 + b3*c^3;

A0 = -a0 - a1*c - a2*c^2 - a3*c^3;
A1 = -3*a0 - a1*c + a2*c^2 + 3*a3*c^3;
A2 = -3*a0 + a1*c + a2*c^2 - 3*a3*c^3;
A3 = -a0 + a1*c - a2*c^2 + a3*c^3;

Bz = [B0 B1 B2 B3];
Az = [A0 A1 A2 A3];

[Hz,wz] = freqz(Bz,Az,length(ws));

% filter Audio
%audio = filter(Bz,Az,audio);
% ^ this operation moved to end.

clearvars -except Bz Az audio Fs Ts vol

%% Calculate Tube Simulation

%Define Project Constants and Matrices
% Volume Variable - Position of Pot 0 - 1;
load('g_map.mat');

%Circuit Constants
C1 = 100e-12;
C2 = 22e-6;
Rfat = 100+3;

% Create Equation (19) for Tube Preamp Including Volume
G = [
1/100,-1/100,0,0,0,0,0,0,-1,0,0,0;...
-1/100,	1/100+1/((1-vol)*1e+6),-1/((1-vol)*1e+6),0,0,0,0,0,0,0,C1,0;...
0,-1/((1-vol)*1e+6),1/(130e+3)+1/((vol)*1e+6)+1/((1-vol)*1e+6)+1/(10e+3),0,-1/(10e+3),0,0,0,0,0,-C1,0;...
0,0,1/(10e+3),0,-1/(10e+3),0,0,0,0,0,0,0;...
0,0,0,0,0,-1/(100e+3),0,1/(100e+3),0,0,0,0;...
0,0,0,0,0,0,1/(1.5e+3),0,0,0,0,C2;...
0,0,0,-1/Rfat,0,0,0,0,0,0,0,C2;...
0,0,0,0,0,1/(100e+3),0,-1/(100e+3),0,-1,0,0;...
1,0,0,0,0,0,0,0,0,0,0,0;...
0,0,0,0,0,0,0,1,0,0,0,0;...
0,1,-1,0,0,0,0,0,0,0,0,0;...
0,0,0,-1,0,0,1,0,0,0,0,0];

G_inv = inv(G);

M1 = zeros(12,2);
M1(11,1) = 1;
M1(12,2) = 1;

M2 = zeros(12,2);
M2(10,1) = 1;
M2(9,2) = 1;

M3 = zeros(12,2);

M3(4,1) = 1;
M3(5:6,2) = 1;
M3(6,1) = 1;

xmat = G_inv*M1;
umat = G_inv*M2;
imat = G_inv*M3;

A = xmat(11:12,:);
B = umat(11:12,:);
C = imat(11:12,:);

D = [xmat(5,:)-xmat(7,:); xmat(6,:)-xmat(7,:)];
E = [umat(5,:)-umat(7,:); umat(6,:)-umat(7,:)];
F = [imat(5,:)-imat(7,:); imat(6,:)-imat(7,:)];

L = xmat(6,:);
M = umat(6,:);
N = imat(6,:);

xnm1 = [0;0];
x2nm1 = [0;0];

H = inv((Fs*eye(2) - A));
K = D*H*C + F;

for k = 1:length(audio);
    tic
    un = [238;audio(k)];
    
    pn = Fs*D*H*xnm1 + (D*H*B + E)*un;
    pn_save(:,k) = pn;
    
    pn_map = ceil(pn_map_func(pn));
    pn_map_actual = pn_map_func(pn);
    
    x1 = pn_map(1);
    x2 = x1-1;
    y1 = pn_map(2);
    y2 = y1-1;
    
    %Bilinear Interpolate P_map
    a1 = inv([1 x1 y1 x1*y1;...
        1 x1 y2 x1*y2;...
        1 x2 y1 x2*y1;...
        1 x2 y2 x2*y2])*[g_map(1,x1,y1);g_map(1,x1,y2);g_map(1,x2,y1);g_map(1,x2,y2)];
        
    a2 = inv([1 x1 y1 x1*y1;...
        1 x1 y2 x1*y2;...
        1 x2 y1 x2*y1;...
        1 x2 y2 x2*y2])*[g_map(2,x1,y1);g_map(2,x1,y2);g_map(2,x2,y1);g_map(2,x2,y2)];
    
    in = [a1(1) + a1(2)*pn_map_actual(1) + a1(3)*pn_map_actual(2) + a1(4)*pn_map_actual(1)*pn_map_actual(2);...
        a2(1) + a2(2)*pn_map_actual(1) + a2(3)*pn_map_actual(2) + a2(4)*pn_map_actual(1)*pn_map_actual(2)];
    
    xn = Fs*H*xnm1 + H*(B*un + C*in);
    
    xnm1 = xn;

    
    output(k) = L*xn + M*un + N*in;
end

%output = filter(Bz,Az,output);

% Take away DC and normalize to 1.
output = (output(10:end) - mean(output(10:end)));
output = output/max(abs(output));
end

