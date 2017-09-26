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

% Load TIR data
% spectral cut off to 1Hz
% Load TIR filtered and detrended data
load(strcat(data_path, 'TIR_data/IOP9/surf_temp1.mat'))
TIR.nx = size(surf_temp1,1);
TIR.ny = size(surf_temp1,2);
TIR.freq = 20;
% Wavelet analysis
TIR.patch = squeeze(mean(surf_temp1,3));
figure_opt = 'on';
Wlet_type = 'amor'; %'morlet'

%compute_1D_wavelet( T_prime_cut(:,1) ,.03, Wlet_type, figure_opt)
%% 
Dx = .03; %m
var = TIR.patch(:,1)';

%figure() %In case we want to plot, a 'running' plot at each gate distance.

p = (nextpow2(length(var)))+1;
scales = 2.^(1:p);  %The number of scales is 'for us' given by the 
                    %length of the signal, since we want the points to
                    %be equidistant, in order to have an homogeneous
                    %representation of the spectra over all the
                    %frequencies.

Fs = scal2frq(scales,'morl',1/Dx); %Computes the "Pseudo-frequency" associated
                                 %with each Wavelet at each scale.

% ------------------------------------
%Padding the Signal, avoids end/border effects:
diff = (2^p)-length(var);
tmpU = padarray(var,[0 floor(diff/2)],'pre');
tmpdiff = diff - floor(diff/2);
tmpU = padarray(tmpU,[0 tmpdiff],'post');
% ------------------------------------

%Computing the Wavelet of the actual signal.
cwt(tmpU,Wlet_type,1/Dx);
[coefs,f] = cwt(tmpU,Wlet_type,1/Dx);
Energy = (abs(coefs)).^2;
tmpSpectra = (mean(Energy,2))';
%
%Grouping all the transforms along the different gates in a single
%Variable.

for j = 1:length(tmpSpectra)
    Spectra(j) = tmpSpectra(j);
    MSpectra(j) = f(j)'*tmpSpectra(j);
end
%
%IN case we want to plot, a 'running' plot at each gate distance.
switch figure_opt
    case 'on'
        figure()
        semilogx(f,f'.*tmpSpectra,'-k')
        ylabel('E$f$','interpreter','latex','fontsize',20)
        xlabel('$f (Hz)$','interpreter','latex','fontsize',20)
        hold on;
        
    case 'off' 
        %do nothing            
end



%% Spectra in a log-log plot:
%---------------------------------------
avgSpectra = squeeze(mean(Spectra,1));

figure()
loglog(Fs,avgSpectra,'-k')

f = 0.0007:0.001:0.06; %We add the k^-1 section.
%f = 0.01:0.001:0.06; %We add the k^-1 section.
hold on;loglog(f,0.000074e-30*(f.^(-1)),'-b')

ff = 0.05:0.01:0.25;  %We add the k^-(5/3) section.
hold on;%loglog(ff,0.001*(ff.^(-5/3)),'-r')
ylabel('$|Y(f)|$','Interpreter','latex','fontsize',14,'FontName','Arial');
xlabel('$f\,[Hz]$','Interpreter','latex','fontsize',14,'FontName','Arial');

%Premultiplied Spectra in a semilogx plot:
%------------------------------------------

figure()
semilogx(f,f.*avgSpectra,'-k')
ylabel('$f\cdot|Y(f)|$','Interpreter','latex','fontsize',14,'FontName','Arial');
xlabel('$f\,[Hz]$','Interpreter','latex','fontsize',14,'FontName','Arial');







