% PhaseRandomWalk.m
% see [Roberts2009], [Barry1990]

clear variables
FontSize = 9;
% s = RandStream('mt19937ar','5Seed',6); RandStream.setGlobalStream(s);

%%
N = 2^8; % no. of PN processes
M = 2^16; % no. of symbols

LW = 0.1e6; % [Hz] 2 x X MHz laser linewidths (Tx & LO)

Rs = 20e9; % [Baud] symbol rate (i.e. sampling freq. of random trials: must be hifgh enough to prevent aliasing in FFT)
dt = 1/Rs; % [s] symbol duration = time step for discrete phase noise model (1 SpS)
T = M*dt;  % [s] observation duration
t = 0:dt:T-dt; %[s] time axis

varPN = 2*pi*LW*dt; %[rad²] [Roberts2009]/eq.(1), [Barry1990]/eq.(112), phase variance 
                    % proportional to LW and observation time
                    
rmsDev = 180/pi*sqrt(varPN); % [degree/symbol] RMS phase deviation     
disp(['RMS phase deviation per symbol = ',num2str(rmsDev,'%1.4f'),' [deg/symbol]'])
                    
%% generate an ensemble of N phase noise (Wiener) processes
rng(1);
phaseInc = sqrt(varPN)*randn(M,N); % [rad] Gaussian i.i.d. phase increments with variance varPN
phaseAcc = cumsum(phaseInc,1); % [rad] accumulated phase = random walk
phaseAcc(1,:) = 0;

% brownian bridge?
% missing value to multiple of 2*pi for last step
missStep  =  round( phaseAcc(end,:)/(2*pi))*2*pi - phaseAcc(end,:);
% generate linear trends to add to each random walk
for k=1:N
%       phaseAcc(:,k) = phaseAcc(:,k) + linspace(0,missStep(k),M)';
end    

varMeas = var(phaseAcc,0,2); % estimated phase variance at every time 
                             % instant accross all sampled random processes

%% instantanous frequency and PSD from phase noise --> [Barry1990]/eqn.(107)ff
wc = 0*2*pi*Rs/8; % laser carrier frequency (on freq. grid!)
Ps = 1e-3; % laser lightwave power [W]
x = sqrt(Ps)*exp(1i*(repmat(wc*t(:),1,N) + phaseAcc)); % N sample functions
% ACF = ifft(fft(x) .* conj(fft(x))); % (cyclic) autocorrelation fct. of process x

win = blackmanharris(size(x,1),'symmetric');
win_sum = trapz(1/length(win),win.^2);
x = x .* repmat(win,1,size(x,2)) /sqrt(win_sum);

PSD = fftshift( fft(x).*conj(fft(x))/length(x)^2 ); % = fft(ACF(x))

%freq = fftshift([-M/2:M/2-1]'/T);
freq = [-M/2:M/2-1]'/T; % valid only for even number of samples!

%% analytical Lorentzian PSD [Barry1990]/eq.116
% TODO: lw = max(1/T,LW); % lower bound LW to FFT frequency resolution
Sx = 2*Ps/(pi*LW)./(1+((2*pi*freq-wc)/(pi*LW)).^2) * 1/T; % [W/Hz]*1/T --> W (in df)

%% Plots
% t = 0:length(t)-1;

f = figure(1); clf
subplot(221) % Ensemble of Wiener PN processes
plot(t/1e-6,180/pi*phaseAcc), grid on
hc=gca; hc.FontSize=FontSize-1;
% xlabel('symbol index','FontSize',FontSize)
xlabel('Time [µs]','FontSize',FontSize)
ylabel('acc. phase [deg]','FontSize',FontSize)
title('Ensemble of Wiener PN processes (randwom walk)','FontSize',FontSize-1)
% axis([0 max(t) -max(abs(phaseAcc(:)))*180/pi*1.1 max(abs(phaseAcc(:)))*180/pi*1.1])
max_y = max(1.1 * max(abs(phaseAcc(:)))*180/pi, 0.1);
axis([0 max(t/1e-6) -max_y max_y])

subplot(222) % Accumulated RMS phase deviation
plot(t/1e-6,180/pi*sqrt(varPN*(0:length(t)-1)),'r:', t/1e-6,180/pi*sqrt(varMeas),'-','Linewidth',2), grid on
hc=gca; hc.FontSize=FontSize-1;
%xlabel('symbol index','FontSize',FontSize)
xlabel('Time [µs]','FontSize',FontSize)
ylabel('RMS phase deviation [deg]','FontSize',FontSize)
title('Accumulated RMS phase deviation','FontSize',FontSize-1)
% axis([0 max(t) 0 Inf])
axis([0 max(t/1e-6) 0 Inf])

subplot(223) % PSDs of all sample processes
plot(freq/1e9,10*log10(max(PSD,1e-300)/1e-3))
hold on, plot(freq/1e9, 10*log10(Sx/1e-3),'r' ,'linewidth',1.5); hold off
axis([-M/2/T/1e9 M/2/T/1e9 -70 10*log10(10*Ps/1e-3)])
hc=gca; hc.FontSize=FontSize-1;
xlabel('Frequency [GHz]','FontSize',FontSize), ylabel('PSD [dBm]','FontSize',FontSize)
title('PSDs of all sample processes','FontSize',FontSize-1)

subplot(224) % avg. PSD from all phase noise processes
%plot(freq/1e9,10*log10(mean(PSD,2)/1e-3),'linewidth',2)
p1=plot(freq/1e9, 10*log10(mean(PSD,2)/1e-3) ,'linewidth',1.5); hold on
p2=plot(freq/1e9, 10*log10(Sx/1e-3),'r' ,'linewidth',1.5); hold off
span_x = max(1000*LW,1/T);
%span_x=Rs/2;
grid on, axis([(wc/2/pi-span_x)/1e9 (wc/2/pi+span_x)/1e9 -90 10*log10(10*Ps/1e-3)])
hc=gca; hc.FontSize=FontSize-1;
% axis([-Inf Inf -70 10*log10(1.1*Ps/1e-3)])
xlabel('Frequency [GHz]','FontSize',FontSize), ylabel('PSD [dBm]','FontSize',FontSize)
title('Zoom to avg. PSD from all phase noise processes','FontSize',FontSize-1)
hl=legend([p2,p1],'Lorentzian PSD','mean(sample PSDs)','Location','south');hl.FontSize=FontSize-1;hl.Box='off';