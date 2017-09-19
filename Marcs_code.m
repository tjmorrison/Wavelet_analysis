%Wavelet Analysis on S%% Wavelet analysis of the Original Signal:

Wlet='haar'; %'cmor2-1'; %Type of Wavelet that we want to use.
Dt = 2;      %[s] Correpsonds to the interpolated time-line.
Fs = 1/Dt;   %Sampling Frequency.


%figure() %In case we want to plot, a 'running' plot at each gate distance.

%This "i" index is used for the Radial (along the LiDar beam) distance.
% It basically controls how far downstream do we want to go for the radial
% averaging of the spectra. The furthest we go downstream, the more noise
% we introduce in the spectra.

for i = 30:40%size(RecField,1) 
%for i = 1:size(OSignal,1)
%for i =  30:40 %2:10 %13:23 %2:45 46:50 %30:40
    
    %u = RecField(i,:); %If we want to prduce the Wavelet Spectra of a
                        %given POD reconstructed Signal.
    u = OSignal(i,:);   %For the Wavelet Transform of the Origianl Signal. 
    %u = FiltData(i,:);   %For the Wavelet Transform of the HIgh&Low pass filtered Signal. 
    
    
    p = (nextpow2(length(u)))+1;
    scales = 2.^(1:p);  %The number of scales is 'for us' given by the 
                        %length of the signal, since we want the points to
                        %be equidistant, in order to have an homogeneous
                        %representation of the spectra over all the
                        %frequencies.
                        
    freq = scal2frq(scales,Wlet,Dt); %Computes the "Pseudo-frequency" associated
                                     %with each Wavelet at each scale.

    %------------------------------------
    %Padding the Signal, avoids end/border effects:
    diff = (2^p)-length(u);
    tmpU = padarray(u,[0 floor(diff/2)],'pre');
    tmpdiff = diff - floor(diff/2);
    tmpU = padarray(tmpU,[0 tmpdiff],'post');
    %------------------------------------
    
    %Computing the Wavelet of the actual signal.
    coefs = cwt(tmpU,scales,Wlet);
    Energy = (abs(coefs)).^2;
    tmpSpectra = (mean(Energy,2));
    
    %Grouping all the transforms along the different gates in a single
    %Variable.
    
    for j = 1:length(tmpSpectra)
        Spectra(i,j) = tmpSpectra(j);
        MSpectra(i,j) = freq(j)*tmpSpectra(j);
    end
    
    %IN case we want to plot, a 'running' plot at each gate distance.
    %semilogx(freq,freq'.*tmpSpectra,'-k')
    %hold on;
    %pause;
    
    
end

%Color Plot of the Premultiplied Spectra
%---------------------------------------
x = log10(freq)';            %Pseudo Frequnce axis.
y = (1:size(Spectra,1))'*dx; %Gate Distance axis 

figure()
pcolor(x,y,MSpectra);shading interp;
colorbar;
xlabel('log(f) ','fontsize',14,'FontName','Arial');
ylabel('r [m] ','fontsize',14,'FontName','Arial');

%Spectra in a log-log plot:
%---------------------------------------
avgSpectra = squeeze(mean(Spectra,1));

figure()
loglog(freq,avgSpectra,'-k')

f = 0.0007:0.001:0.06; %We add the k^-1 section.
%f = 0.01:0.001:0.06; %We add the k^-1 section.
hold on;loglog(f,0.0074*(f.^(-1)),'-b')

ff = 0.05:0.01:0.25;  %We add the k^-(5/3) section.
hold on;loglog(ff,0.001*(ff.^(-5/3)),'-r')
ylabel('$|Y(f)|$','Interpreter','latex','fontsize',14,'FontName','Arial');
xlabel('$f\,[Hz]$','Interpreter','latex','fontsize',14,'FontName','Arial');

%Premultiplied Spectra in a semilogx plot:
%------------------------------------------

figure()
semilogx(freq,freq.*avgSpectra,'-k')
ylabel('$f\cdot|Y(f)|$','Interpreter','latex','fontsize',14,'FontName','Arial');
xlabel('$f\,[Hz]$','Interpreter','latex','fontsize',14,'FontName','Arial');


% Plotting the same spectra in Wavenumbers:
%------------------------------------------
Z = 15;
meanU = mean(Umean(30:40));
wnum = (2*pi/meanU)*freq;
xAx = wnum*Z;

%Downstream check of the inclination of the beam:
gate = 30; %The good spectra was obtained betw. gates 30 and 40.
theta_shift = 2 %[deg]
Zshift = (gate*18)*sin(theta_shift*(pi/180))


figure()
loglog(xAx,avgSpectra,'-k')
ylabel('$E_{u}(k)$','Interpreter','latex','fontsize',14,'FontName','Arial');
xlabel('$kz$','Interpreter','latex','fontsize',14,'FontName','Arial');

wnum_1 = (2*pi/meanU)*f*Z;
wnum_2 = (2*pi/meanU)*ff*Z;

hold on;loglog(wnum_1,0.0074*(f.^(-1)),'-b')
hold on;loglog(wnum_2,0.001*(ff.^(-5/3)),'-r')

line([1 1], [0.01 100],'LineStyle','--','color','k','LineWidth',1)
line([0.015 0.015], [0.01 100],'LineStyle','--','color','k','LineWidth',1)


%% Pwelch of the Original Signal:

clear f Pxx PDens avgPxx

Dt = 2; %[s] Correpsonds to the interpolated time-line.
fs = 1/Dt; %DT = 2s.


windw = round(size(u,2)/4);
overlap = floor(windw/2);
    
for i = 1:size(OSignal,1)
    
    %Sig = OSignal(i,:)';
    Sig = RecField(i,:);
    
    [Pxx f] = pwelch(Sig,windw,overlap,windw,fs,'onesided');
    
    for j = 1:length(Pxx)                       
        PDens(i,j) = Pxx(j);
    end
end

gate_start = 1;
gate_end = 1;
avgPxx=mean(PDens(gate_start:gate_end,:),1);

figure()
semilogx(f',f'.*avgPxx,'b')
ylabel('$f\cdot|Y(f)|$','Interpreter','latex');
xlabel('$f\,[Hz]$','Interpreter','latex');

figure()
loglog(f',avgPxx,'b')
ylabel('$|Y(f)|$','Interpreter','latex');
xlabel('$f\,[Hz]$','Interpreter','latex');


%% Looking at the FFT of the Original Data.

for i = 1:size(OSignal,1)
    
    %u = OSignal(i,:);
    u = RecField(i,:);
    
    L = length(u);
    NFFT = 2^nextpow2(size(u,2)); % Next power of 2 from length of y
    Y = fft(u,NFFT);
    
    tmpPDens = (Y.*conj(Y))/NFFT;
    
    for j = 1:length(tmpPDens)
        PDens(i,j) = tmpPDens(j);
    end
end

Dt = 2; %[s] Correpsonds to the interpolated time-line.
Fs = 1/Dt; %DT = 2s.
f = (Fs/2)*linspace(0,1,NFFT/2+1);

avgPDens = mean(PDens,1);

figure()
loglog(f(2:end),avgPDens(2:NFFT/2+1),'-k')
ylabel('$|Y(f)|$','Interpreter','latex');
xlabel('$f\,[Hz]$','Interpreter','latex');


%% Hilbert Transform:

it = 1;
n = 5;
inputSig = squeeze(Proj(it,n,shift+1:end-shift));

N = size(inputSig,1);
w = window(@gausswin,N,2.5); 
tmp = w.*inputSig;

inputSig = tmp;


%%
HT_signal=hilbert(inputSig);
Phase=angle(HT_signal);
MODULUS=abs(HT_signal); % MODULUS^2 = ENERGY

x=real(HT_signal);
y=imag(HT_signal);


Nperiodo = size(inputSig);
Fsamp = 1/Dt;

for j=3:Nperiodo-2
    Delta(j)=(12*(x(j)^2+y(j)^2))^-1*(x(j)*(y(j-2)-8*y(j-1)+8*y(j+1)-y(j+2))-y(j)*(x(j-2)-8*x(j-1)+8*x(j+1)-x(j+2))); % Instantaneous Frequency from Eq. 7 in the paper
end
Fhil=(2*pi)^-1*Fsamp*Delta; % Instantaneous Frequency from Eq. 7 in the paper

%% HLCC:

HT_signal1 = HT_signal;
Phase1 = Phase;
MODULUS1 = MODULUS;
FHil1 = Fhil;

HT_signal2 = HT_signal;
Phase2 = Phase;
MODULUS2 = MODULUS;
FHil2 = Fhil;



tmp = (Phase1 - Phase2)./(MODULUS1.*MODULUS2);
figure();plot(tmp)

%% High & low pass filtering the Data:

%A) Low pass filter:

Fs = 1/Dt;
Fhighest = 0.06226; %for low pass filter. (Allows lower frequencies to pass and kills the higher ones). 
Flowest = 0.0009728; %for high pass filter. (Allows higher frequencies to pass, and not to the smaller ones). 

FILTER_h = 'high';
FILTER_l = 'low';

for i = 1:size(OSignal,1)

    u = OSignal(i,:)';
    
    Sig = [flipud(u) u flipud(u)];
    
    [B,A]=butter(4,Flowest/(Fs/2),FILTER_h); %4 is the order of the Butterworth filter.
    filt_data=filtfilt(B,A,Sig);
    
    
    Sig = filt_data;
    
    [B,A]=butter(4,Fhighest/(Fs/2),FILTER_l); %4 is the order of the Butterworth filter.
    filt_data=filtfilt(B,A,Sig);
    
    tmp = length(filt_data);
    for j = 1:length(u)
        FiltData(i,j) = filt_data(tmp+j);
    end
end
    onic data