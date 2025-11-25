Frequency-Division-Multiplexing-
Aim:
To study and implement Frequency Division Multiplexing (FDM) and Demultiplexing using SCILAB by combining six different message signals into a single composite signal for transmission and then recovering each message signal at the receiver through demodulation and filtering.
Apparatus Required:
Computer system with SCILAB software installed.
Theory:
Frequency Division Multiplexing (FDM) is a technique in which multiple message signals are transmitted simultaneously over a single communication channel by assigning each signal a unique carrier frequency. Each message is modulated with its respective carrier so that their frequency bands do not overlap. These modulated signals are then combined to form a single multiplexed signal for transmission. At the receiver end, the signal is demultiplexed by using the same carrier frequencies for demodulation, followed by low-pass filtering to recover the original baseband signals. FDM is widely used in radio broadcasting, cable television, and satellite communication systems where efficient bandwidth utilization is essential.
Algorithm:
1.Start the program and initialize the sampling frequency fs and time vector t.
2.Generate six message signals of different frequencies (150 Hz to 900 Hz) using sine functions.

3.Assign six distinct carrier frequencies (3 kHz to 13 kHz) to each message signal to avoid overlap.

4.Modulate each message signal by multiplying it with its corresponding carrier (Amplitude Modulation).

5.Add all modulated signals to obtain the combined FDM signal representing multiple channels.

6.Plot all six message signals and the multiplexed signal for observation.

7.Demodulate each signal by multiplying the FDM signal with its corresponding carrier frequency.

8.Apply a low-pass filter to extract the original baseband signals (recovering the messages).

9.Plot the demodulated signals to verify successful recovery.

10.End the program after confirming proper multiplexing and demultiplexing operation.

Code:
clc;
clear;
close;

// ---------- settings ----------
PI = 3.14;                // use 3.14 as requested
fs = 80000;               // high sampling rate -> smooth plots
t = 0:1/fs:0.01;          // time vector

fm = [200, 400, 600, 800, 1000, 1200];       // message freqs (Hz)
fc = [5000, 7000, 9000, 11000, 13000, 15000]; // carrier freqs (Hz)

// ---------- generate message signals ----------
m = zeros(length(fm), length(t));
for i = 1:length(fm)
    m(i, :) = sin(2 * PI * fm(i) * t);
end

// ---------- multiplex (AM with cosine carriers) ----------
fdm = zeros(1, length(t));
for i = 1:length(fc)
    carrier = cos(2 * PI * fc(i) * t);
    fdm = fdm + m(i, :) .* carrier;
end

// ---------- design FIR low-pass filter (sinc * Hamming) ----------
fc_cut = 2000;           // cutoff 2 kHz (pass message band)
M = 200;                 // half filter order -> filter length = 2*M+1 = 401
n = -M:M;
wc = 2 * fc_cut / fs;    // normalized (0..1) where 1 => Nyquist

// ideal sinc (note: use PI for pi)
sinc = zeros(n);
for k = 1:length(n)
    x = 2 * fc_cut * n(k) / fs;               // x = 2*fc_cut*n/fs
    if x == 0 then
        sinc(k) = 1;
    else
        // sinc(x) = sin(pi * x) / (pi * x)
        sinc(k) = sin(PI * x) / (PI * x);
    end
end

// ideal impulse response scaled by cutoff (2*fc_cut/fs)
h_ideal = (2 * fc_cut / fs) * sinc;

// apply Hamming window
w = 0.54 - 0.46 * cos(2 * PI * (n + M) / (2*M)); // Hamming length 2M+1
h = h_ideal .* w;

// normalize filter gain at DC
h = h / sum(h);

// ---------- demodulate (coherent detection + LPF) ----------
demod = zeros(length(fm), length(t));
for i = 1:length(fc)
    x = fdm .* cos(2 * PI * fc(i) * t);     // mix to baseband
    y = conv(x, h, 'same');                 // filter (same length)
    demod(i, :) = y;
end

// ---------- plots ----------
scf(1); clf;
for i = 1:size(m,1)
    subplot(size(m,1),1,i);
    plot(t, m(i,:));
end

scf(2); clf;
plot(t, fdm);

scf(3); clf;
for i = 1:size(demod,1)
    subplot(size(demod,1),1,i);
    plot(t, demod(i,:));
end

disp("Done: smooth demodulated signals.");
Output Waveform
image image image
WhatsApp Image 2025-11-21 at 12 39 25_2aa2aaa9

WhatsApp Image 2025-11-21 at 12 39 24_6683fa45

Result:
The Frequency Division Multiplexing (FDM) and Demultiplexing of six different message signals were successfully implemented using SCILAB.All six message signals were modulated with distinct carrier frequencies and combined to form a single multiplexed signal. Upon demodulation and low-pass filtering, each original message signal was accurately recovered, verifying the correct working of the FDM and demultiplexing process.
