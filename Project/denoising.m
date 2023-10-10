clc;
clear all;
clf;
%s signal %f freqquncy of the file  for this it is 44100 hz
[s,fp]=audioread('bird.wav');
%audio file is 1.78 secs long so, 1,78*44100 = 78.5kHz

t_length=length(s);
t=linspace(0,t_length/fp,length(s));
amp=10;
s(1:length(s))=amp*s(1:length(s));
BITS=5;
for n=1:BITS
    a =2^n - 1;
  
 
  signal_quantized = round( s(1:t_length) * a) / a;
  
  
 
  noise = s(1:t_length) - signal_quantized;
  
figure;
plot(t, s(1:t_length), '-b', t, signal_quantized, '-g', t, noise, '-r');
legend('Original signal', 'Quantized signal', 'Error of sampling - noise');
xlabel('Time');
ylabel('Amplitude');
title('Singular Step of Quantization ',int2str(n));
grid on;
  signal_deviation(n,:) = 20 * log10( std( s(1:t_length) ) / std(noise) );
  Quantized_signal  = strcat('Quantized_signal', int2str(n), '.png');
  print(Quantized_signal, '-dpng', '-r300');

   signal_rounded = round(s * a) / a;
   audiowrite(['bird', int2str(n), '.wav'], signal_rounded, fp);
end  


%second frame with plots

subplot(3,3,1);
 plot(t, s(1:length(s)), '-b');
 
 title('Original signal');
 xlabel('Time [s]');
 ylabel('Amplitude [-]');
 grid on;
 hold on;

  subplot(3,3,3);
 plot(signal_deviation, 'r*');
 hold on;
 plot(signal_deviation ,'-b');
 hold on;
 grid on;
  title({'SNR signal to noise ratio', 'dependent on bits of filter'});
 xlabel({'Bits of filter -', 'resolution of quantization'});
 ylabel({'value of SNR', 'in logarytmic scale [dB]'});
  n_poly = 1:1:BITS; 

   polyfit(n_poly', signal_deviation, 1);

    file1 = ['bird' , int2str(BITS), '.wav'];
 file2 = ['bird' , int2str(BITS + 1), '.wav'];
 copyfile(file1, file2);
  [s_final, fp] = audioread(file2);
 N = fp / 2;
%below function is to determine which e power of 2 to
%can have fourier transformerfourier transformer as best way to use fourier
%tranformer is in powert of 2

Nf_help = 1;
    help = 2;
    while(help < N)
        help = help* 2;
        Nf_help=Nf_help+1;
    end

  Nf = 2^Nf_help;
   Nf2 = Nf/2 + 1;
   f = linspace(0, fp/2, Nf2);

    subplot(3,3,4);
 plot(t, s_final(1:t_length));
 title({'Signal quantized','-without noise'});
 xlabel('Time [s]');
 ylabel('Amplitude [-]');
 hold on;
 grid on;
 s_fft = fft(s_final, Nf);
  s_fft_abs = abs(s_fft);

   subplot(3,3,5);
 plot(f, s_fft_abs(1:Nf2) );
 title({"Module of a signal after", "iterative sampling"}); 
 xlabel('Frequenzy [Hz]');
 ylabel({'Module of a signal after', 'iterative sampling [-]'});
 box off; grid on;
  s_fft_angle = angle(s_fft);

 subplot(3,3,6);
 plot(f, s_fft_angle(1:Nf2));
 title({"Phase of signal after", "iterative sampling"});
 xlabel('Frequenzy [Hz]');
 ylabel('Phase angle of a signal [rad]');
 box off; grid on; axis tight;      

  subplot(3,3,7);
 plot(t,  noise(1:t_length) );
 title('Noise of sampling'); 
 xlabel('Time [s]'); 
 ylabel('Amplitude [-]');
 hold on; grid on;    
 noise_fft = fft(noise, Nf);

noise_fft_abs = abs(noise_fft);   


subplot(3,3,8);
plot(f, noise_fft_abs(1:Nf2));
title("Module of noise spectrum");
xlabel('Frequenzy [Hz]');
ylabel('Module of noise spectrum [-]');
grid on; hold on; 

noise_ftt_angle = angle(noise_fft);

subplot(3,3,9);
plot(f, noise_ftt_angle(1:Nf2));
title("Phase angle of noise");
xlabel('Frequenzy [Hz]');
ylabel('Phase angle of signal [rad]');
grid on; hold on;

filename = 'Fig. 1 - Plots of signals original sampled and noise.png';
set(gcf, 'Position', [0, 0, 1600, 900]);
print(filename, '-dpng', '-r300');

freqCut = [2200 2900];
wc = freqCut / (fp/2);
firCoeff = fir1(511, wc, 'bandpass');
subplot(2,3,4);
stem(firCoeff);
title("Filter coefficients");
hold on;


s_filtered = filter(firCoeff, 1, s_final);

subplot(2,3,1);
figure;
plot(t, s_filtered, '-r', t, s_final, '-b');
xlabel('Time');
ylabel('Amplitude');
legend('s_{filtered}', 's_{final}');
title('Plot of s_{filtered} and s_{final}');
grid on;

s_filtered_fft = fft(s_filtered, Nf);

s_filtered_fft_abs = abs(s_filtered_fft);
subplot(2,3,2);
plot(f, s_filtered_fft_abs(1:Nf2));
title({"Module of signal spectrum after" , "using bandpass filer"});
xlabel('Frequenzy [Hz]');
ylabel('Module of signal spectrum after reduction [-]');
grid on; hold on;

s_filtered_fft_angle = angle(s_filtered_fft);

subplot(2,3,3);
plot(f, s_filtered_fft_angle(1:Nf2));
title({"Signal phase angle", "after filter FIR"});
xlabel('Frequenzy [Hz]');
ylabel('Signal phase angle [rad]');
grid on; hold on;

audiowrite(['bird_filtered.wav'], s_filtered, fp); 

[H,f] = freqz(firCoeff, 1, 2^Nf_help, fp);

subplot(2,3,5);
plot(f, abs(H));
title({"Frequenzy", "response"});
xlabel('Frequenzy [Hz]')'
ylabel('Frequenzy response [-]');
grid on; hold on;

subplot(2,3,6);
plot(f, angle(H));
title({"Frequenzy response of", "phase angle"});
xlabel('Frequenzy [Hz]');
ylabel('Phase angle of signal [rad]');
grid on; hold on;

filename = 'Fig. 2 - Plots of filter FIR.png';
set(gcf, 'Position', [0, 0, 1600, 900]);
print(filename, '-dpng', '-r300');

figure;
freqz(firCoeff, 1, 2^Nf_help, fp);

filename = 'Fig. 3 - Plots of frequency response of filter FIR.png';
set(gcf, 'Position', [0, 0, 1600, 900]);
print(filename, '-dpng', '-r300');

% [s,fp]=audioread('bird.wav')
% sound(s,fp);
% [s,fp]=audioread('bird_filtered.wav')
% sound(s,fp);


