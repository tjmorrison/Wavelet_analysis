%Compute Wavelet on TIR
clear; clc; close all;

machine = 'imac_local';
%add paths
switch machine
    case 'imac_local'
        %code
        addpath('/Users/travismorrison/Documents/Code/functions_library');
        %data path
        data_path ='/Users/travismorrison/Local_Data/MATERHORN/data/';
        addpath(data_path); 
    case 'imac_Achilles'
        %code
        addpath('/Users/travismorrison/Documents/Code/functions_library');
        %data path
        data_path ='/Volumes/ACHILLES/MATERHORN/data/';
        addpath(data_path); 
    case 'linux'
        %code
        addpath('/home/tjmorrison/Research/function_library');
        %data path
        data_path = addpath('/media/tjmorrison/ACHILLES/MATERHORN/data');
end
% Set plot stuff
ft_size = 25;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex'); 
set(0,'DefaultAxesFontSize',ft_size); 

%% Load TIR data
% spectral cut off to 1Hz
% Load TIR filtered and detrended data
load('/Users/travismorrison/Documents/Code/TIR_xwire_spectral_analysis/IOP9_60minDetrend_1HzFilter.mat')
TIR.nx = size(T_prime_cut,1);
TIR.ny = size(T_prime_cut,2);
TIR.freq = 20;
%% Wavelet analysis
T_patch = squeeze(mean(T_prime_cut,3));
figure_opt = 'on';
Wlet_type = 'morl'; %'morlet'

compute_1D_wavelet( T_prime_cut(:,1) ,.03, Wlet_type, figure_opt)