%Compute pod
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
% Plot xire and TIR Spectras
ft_size = 25;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex'); 
set(0,'DefaultAxesFontSize',ft_size); 
%% Load tower data
load(strcat(data_path,'tower_data/Playa_tower_raw/PlayaSpring_raw_GPF_LinDet_2013_05_24.mat'));
IOP9_start = 3236840;
IOP9_end = IOP9_start+(20*60*60); 
tower.w = rearrangeHeights(detrend(inpaint_nans(rawFlux.wPF(IOP9_start:IOP9_end,:))));
tower.T= rearrangeHeights(detrend(inpaint_nans(rawFlux.sonTs(IOP9_start:IOP9_end,:))));
clear rawFlux;
tower.num_z = 6; 

% spectral cut off to 1Hz
tower.Fs = 20;
for z = 1:tower.num_z
    [tower.T(:,z)] = spectral_cuttoff(tower.Fs,tower.T(:,z)', 0 ,1);
    [tower.w(:,z)] = spectral_cuttoff(tower.Fs,tower.w(:,z)', 0,1);
end
%% Wavelet analysis
figure_opt = 'on';
Wlet_type = 'haar';
%for z = 1:tower.num_z
compute_1D_wavelet( tower.T(1:120:end,1) ,6, Wlet_type, figure_opt)