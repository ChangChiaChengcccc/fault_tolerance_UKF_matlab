% white Gaussian noise generator
% dBw-dBm-W converter : http://www.elecfans.com/tools/dbm.htm
close all;

% generate white Gaussian noise with 1000*3 matrix, -10dBW(0.099watts)
noise = wgn(1000, 1, -60); %2*10^-3 position
noise = wgn(1000, 1, -50); %5*10^-3 velocity 
noise = wgn(1000, 1, -42); %0.01 omega
% variance of the matrix
var(noise)

% plot the data
plot(noise)
