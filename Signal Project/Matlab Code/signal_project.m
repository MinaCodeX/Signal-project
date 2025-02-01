%1. SIGNALS GENERATION (RAMP, UNIT STEP), (EACH TEAM WORK ON ONLY THE TWO SIGNALS THAT MENTIONED TO IN THE PROJECT SHEET) 

% Parameters
n = -10:10; % Discrete time range

% Ramp Signal
ramp_signal = max(0, n);

% Unit Step Signal
unit_step_signal = double(n >= 0);

% Plot signals
figure;
subplot(2,1,1);
stem(n, ramp_signal, 'LineWidth', 1.5);
title('Ramp Signal');
xlabel('n'); ylabel('Amplitude');

subplot(2,1,2);
stem(n, unit_step_signal, 'LineWidth', 1.5);
title('Unit Step Signal');
xlabel('n'); ylabel('Amplitude');


% 2. ADDITION, SUBTRACTION, MULTIPLICATION OF DIFFERENT SIGNALS. (CHOOSE TWO OF THEM AND WORK ON THE NEW GENERATED SIGNALS)

% Addition
addition = ramp_signal + unit_step_signal;

% Subtraction
subtraction = ramp_signal - unit_step_signal;

% Multiplication
multiplication = ramp_signal .* unit_step_signal;

% Plot results
figure;
subplot(3,1,1);
stem(n, addition, 'LineWidth', 1.5);
title('Addition of Signals');

subplot(3,1,2);
stem(n, subtraction, 'LineWidth', 1.5);
title('Subtraction of Signals');

subplot(3,1,3);
stem(n, multiplication, 'LineWidth', 1.5);
title('Multiplication of Signals');


% 3. SAMPLING OF THE SIGNALS.

% Continuous signal (sine wave example)
t = 0:0.01:2*pi; % Continuous time
x_cont = sin(t);

% Sampled signal
T_s = 0.2; % Sampling interval
t_sampled = 0:T_s:2*pi;
x_sampled = sin(t_sampled);

% Plot continuous and sampled signals
figure;
plot(t, x_cont, 'LineWidth', 1.5); hold on;
stem(t_sampled, x_sampled, 'r', 'LineWidth', 1.5);
title('Sampling of a Signal');
legend('Continuous Signal', 'Sampled Signal');
xlabel('Time'); ylabel('Amplitude');


% 4. CONVOLUTION OF TWO SIGNALS USING FORMULA AND COMPARE THE RESULTS USING "CONVO" COMMAND.

% Example signals
x1 = unit_step_signal;
x2 = ramp_signal;

% Manual convolution
n_conv = 0:length(x1) + length(x2) - 2; % Range for convolution
conv_manual = conv(x1, x2);

% MATLAB convolution
conv_builtin = conv(x1, x2);

% Plot comparison
figure;
subplot(2,1,1);
stem(n_conv, conv_manual, 'LineWidth', 1.5);
title('Manual Convolution');

subplot(2,1,2);
stem(n_conv, conv_builtin, 'r', 'LineWidth', 1.5);
title('MATLAB Convolution (conv command)');



% 5. MAKE A FILTER (FOR EXAMPLE AVERAGING FILTER) AND PASS NOISY SIGNAL THROUGH IT. YOU CAN USE "FILTER" COMMAND FOR FILTRATION OF SIGNALS.

% Noisy signal
noise = rand(1, length(n));
noisy_signal = ramp_signal + noise;

% Averaging filter
M = 3; % Filter size
h = ones(1, M) / M;

% Filter the signal
filtered_signal = filter(h, 1, noisy_signal);

% Plot noisy and filtered signals
figure;
subplot(2,1,1);
stem(n, noisy_signal, 'LineWidth', 1.5);
title('Noisy Signal');

subplot(2,1,2);
stem(n, filtered_signal, 'LineWidth', 1.5);
title('Filtered Signal');



% 6. FIND FREQUENCY RESPONSE OF A SIGNAL USING COMMAND AND FORMULA.

% Frequency response of filter
[H, w] = freqz(h, 1, 512);

% Plot frequency response
figure;
plot(w/pi, abs(H), 'LineWidth', 1.5);
title('Frequency Response of the Averaging Filter');
xlabel('Normalized Frequency (\times \pi rad/sample)');
ylabel('Magnitude');


% 7. CALCULATE DISCRETE TIME FOURIER TRANSFORM AND ITS DIFFERENT PROPERTIES (LIKE CONVOLUTION, TIME SHIFTING).

% DTFT computation
omega = linspace(-pi, pi, 512); % Frequency range
X = fftshift(fft(ramp_signal, 512)); % DTFT using FFT

% Plot DTFT
figure;
plot(omega, abs(X), 'LineWidth', 1.5);
title('Discrete-Time Fourier Transform (DTFT)');
xlabel('Frequency (rad/sample)');
ylabel('Magnitude');



% 8. DISCRETE FOURIER TRANSFORM AND ITS PROPERTIES (E.G. 4,8,16- POINT DFT, DFT WITH ZERO PADDING).

% DFT
N = 8; % Points
dft_signal = fft(ramp_signal, N);

% Plot DFT
figure;
stem(0:N-1, abs(dft_signal), 'LineWidth', 1.5);
title('Discrete Fourier Transform (DFT)');
xlabel('Frequency Index');
ylabel('Magnitude');


% 9. THEN MOVE TO FILTERS AND FILTER A SIGNAL BY PASSING IT THROUGH DIFFERENT...
% 10-CALCULATE Z-TRANS

% Z-transform using symbolic math
syms z;
z_transform = sum(ramp_signal .* z.^(-n));

% Display Z-transform
disp('Z-Transform:');
disp(z_transform);
