% HW 1, Questions 2-3

% I explicitly used Brian Hargreaves' blochsim code to work on this question: http://mrsrl.stanford.edu/~brian/bloch/

%% 2a

% constants
T1 = 1000; % ms
T2 = 100; % ms

% Bloch Equation Simulation, Excercise B-4a

TE = [2.5, 5, 10];	% ms
TR = 2*TE;	% ms
flip = pi/3;

df = -100:100; 	% Hz

Sig = zeros(length(df),length(TE));

for n=1:length(TE)

    for k=1:length(df)
		[Msig,Mss] = sssignal(flip,T1,T2,TE(n),TR(n),df(k));
		Sig(k,n)=Mss(1)+i*Mss(2);
    end

end

% ===== Plot the Results ======

subplot(2,1,1);
plot(df,abs(Sig));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
legend('TE=2.5', 'TE=5', 'TE=10');

subplot(2,1,2);
plot(df,angle(Sig));
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
axis([min(df) max(df) -pi pi]);
grid on;
legend('TE=2.5', 'TE=5', 'TE=10');

%% 2bi

flip = 10*pi/180;
T1 = 1000; % ms
T2 = 100; % ms
TR = 10; % ms
TE = 5; % ms
df = -100:100;
Msig_vec = zeros(length(df),1);
Mss_vec = zeros(length(df),3);
num_cycles = 8;

for i = 1:length(df)
    dfreq = df(i);
    [Msig,Mss] = gresignal(flip,T1,T2,TE,TR,dfreq,num_cycles);
    Msig_vec(i) = Msig;
    Mss_vec(i,:) = Mss;
end

plot(df,abs(Msig_vec))
ylim([0,1])
% the steady state signal is 0.0757

%% 2bii affect of gradient spoiler

dfreq = 0;
n_cyc = [1,2,4,8];
num_cycles_vec = 2.*n_cyc;
Msig_vec_4 = zeros(4,1);
for n = 1:4
    num_cycles = num_cycles_vec(n);
    [Msig,Mss] = gresignal(flip,T1,T2,TE,TR,dfreq,num_cycles);
    Msig_vec_4(n) = Msig;
end

figure
plot(n_cyc,abs(Msig_vec_4),'*');
ylim([0,1])

%% 2biii
% dfreq is the off-resonance
% rf_spoil is the spoiling (phi)

dfreq = 0;
rf_spoil_freq = ([0:1:180]).*pi/180;
Msig_vec_rf_spoil = zeros(length(rf_spoil_freq),1);
for j = 1:length(rf_spoil_freq)
    rf_spoil = rf_spoil_freq(j);
    %[Msig,Mss] = gssignal_iii(flip,T1,T2,TE,TR,dfreq,num_cycles,rf_spoil);
    [Msig,Mss]=spgrsignal(flip,T1,T2,TE,TR,dfreq,100,rf_spoil);

    Msig_vec_rf_spoil(j) = Msig;
end

%figure
rf_spoil_freq_deg = ([0:1:180]);

plot(rf_spoil_freq_deg,abs(Msig_vec_rf_spoil),'*');
%hold on
%plot(abs(Msig_vec))

% RF phase that best eliminates transverse magnetization should be around
% 117 degrees (not sure how to get this from the plot)

%% 3.1 (slice profile)

TBWP = 8; % number of zero crossings
sl_thk_mm = 5;
t_rf_ms = 2;
t_vec_ms = 0:0.01:2;
alpha_rad = pi/2; %10*pi/180; %30*pi/180;
%alpha_rad = 30*pi/180;
%alpha_rad = 10*pi/180;
gamma_bar_MHz_per_T = 42.6;
T1_ms = 1000;
T2_ms = 100; %2
df_Hz = 0; %200;
pos = -20:0.1:20;

rf = sinc(TBWP/2.*(t_vec_ms-1));
B1 = alpha_rad/(2*pi*gamma_bar_MHz_per_T.*sum(rf));
rf = rf.*B1*1000;
% plot(t_vec_ms,rf);
% grid on
grad_G_per_cm = 1.88.*ones(1,length(t_vec_ms)); % calculated on paper
t_grad_sec = 2.208*10^-3; 

%[msig,m]=sliceprofile(rf,grad_G_per_cm,t_vec_ms./1000,T1_ms,T2_ms,pos,df_Hz);

% add rephasing gradient
rf = [rf 0*rf];
grad_G_per_cm = [grad_G_per_cm -grad_G_per_cm/2];
t_vec_ms = [t_vec_ms t_vec_ms+2+0.01];

% 3.5 SMS
rf_SMS = rf;
slice_pos = [0, 20, 40, 60, 80]; % mm
for SMS_t = 1:length(rf)
    P = 0;
    for SMS_i = 1:5
        P = P + exp(1i*slice_pos(SMS_i)*rf(i));
    end
    rf_SMS(SMS_t) = rf(SMS_t).*P;
end

%[msig,m]=sliceprofile(rf_SMS,grad_G_per_cm,t_vec_ms./1000,T1_ms,T2_ms,pos,df_Hz);
[msig,m]=sliceprofile(rf,grad_G_per_cm,t_vec_ms./1000,T1_ms,T2_ms,pos,df_Hz);


% 3.2 plots (do for 0 and 200 Hz off-resonance)
%  m and msig are 3xN and 1xN arrays of the magnetization and signal at
%  each point in pos
% figure
% plot(t_vec_ms,rf);
% grid on
% figure
% hold on
plot(pos,abs(msig));
%plot(pos,angle(msig));
% figure
% plot(t_vec_ms,grad_G_per_cm);
% figure
% plot(pos,m(1,:));
% figure
% plot(pos,m(2,:));
% figure
% plot(pos,m(3,:));

%3.3 change the flip angle
%legend("90 deg", "30 deg", "10deg");
%fft_sp = fftshift(abs(fft(rf)));
%plot(pos,fft_sp(1:end-1));
%legend("10 deg sp", "10 deg fft", "90 deg sp", "90 deg fft");
%legend("T2=100ms", "T2=2ms");

%3.4 take out rephasing edits
%legend("no rephase", "rephase");


%% subfunctions explicitly from http://mrsrl.stanford.edu/~brian/bloch/
% made some edits to solve the homework questions, but base code from here

function Rz=zrot(phi)
Rz = [cos(phi) -sin(phi) 0;sin(phi) cos(phi) 0; 0 0 1];
end

function Rx=xrot(phi)
Rx = [1 0 0; 0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];
end

function Ry=yrot(phi)
Ry = [cos(phi) 0 sin(phi);0 1 0;-sin(phi) 0 cos(phi)];
end

function Rth=throt(phi,theta)
Rz = zrot(-theta);
Rx = xrot(phi);
Rth = inv(Rz)*Rx*Rz;
end

function [Afp,Bfp]=freeprecess(T,T1,T2,df)
%
%	Function simulates free precession and decay
%	over a time interval T, given relaxation times T1 and T2
%	and off-resonance df.  Times in ms, off-resonance in Hz.
phi = 2*pi*df*T/1000;	% Resonant precession, radians.
E1 = exp(-T/T1);	
E2 = exp(-T/T2);
Afp = [E2 0 0;0 E2 0;0 0 E1]*zrot(phi);
Bfp = [0 0 1-E1]';
end

% 
%	function [Msig,Mss] = sssignal(flip,T1,T2,TE,TR,dfreq)
% 
%	Calculate the steady state signal at TE for repeated
%	excitations given T1,T2,TR,TE in ms.  dfreq is the resonant
%	frequency in Hz.  flip is in radians.
%

function [Msig,Mss] = sssignal(flip,T1,T2,TE,TR,dfreq)

Rflip = yrot(flip);
[Atr,Btr] = freeprecess(TR-TE,T1,T2,dfreq);
[Ate,Bte] = freeprecess(TE,T1,T2,dfreq);

% Let 	M1 be the magnetization just before the tip.
%	M2 be just after the tip.
%	M3 be at TE.
%
% then
%	M2 = Rflip * M1
%	M3 = Ate * M2 + Bte
%	M1 = Atr * M3 + Btr
%
% Solve for M3...
%
%	M3 = Ate*Rflip*Atr*M3 + (Ate*Rflip*Btr+Bte)

Mss = inv(eye(3)-Ate*Rflip*Atr) * (Ate*Rflip*Btr+Bte);
Msig = Mss(1)+i*Mss(2);
end


% 
%	function [Msig,Mss] = gssignal(flip,T1,T2,TE,TR,dfreq,phi)
% 
%	Calculate the steady state signal at TE for repeated
%	excitations given T1,T2,TR,TE in ms.  dfreq is the resonant
%	frequency in Hz.  flip is in radians.
%	phi is the phase twist at the end of the sequence.

function [Msig,Mss] = gssignal(flip,T1,T2,TE,TR,dfreq,phi)

Rflip = yrot(flip);
[Atr,Btr] = freeprecess(TR-TE,T1,T2,dfreq);
[Ate,Bte] = freeprecess(TE,T1,T2,dfreq);

% 	To add the gradient spoiler twist, we just
%	multiply Atr by zrot(phi):

Atr = zrot(phi)*Atr;

% Let 	M1 be the magnetization just before the tip.
%	M2 be just after the tip.
%	M3 be at TE.
%
% then
%	M2 = Rflip * M1
%	M3 = Ate * M2 + Bte
%	M1 = Atr * M3 + Btr
%
% Solve for M3...
%
%	M3 = Ate*Rflip*Atr*M3 + (Ate*Rflip*Btr+Bte)

Mss = inv(eye(3)-Ate*Rflip*Atr) * (Ate*Rflip*Btr+Bte);
Msig = Mss(1)+i*Mss(2);

end


% 
%	function [Msig,Mss] = gssignal(flip,T1,T2,TE,TR,dfreq,phi)
% 
%	Calculate the steady state signal at TE for repeated
%	excitations given T1,T2,TR,TE in ms.  dfreq is the resonant
%	frequency in Hz.  flip is in radians.
%	phi is the phase twist at the end of the sequence.

% edited
function [Msig,Mss] = gssignal_iii(flip,T1,T2,TE,TR,dfreq,phi,rf_phase)

Rflip = throt(flip,rf_phase);
%disp(Rflip)

[Atr,Btr] = freeprecess(TR-TE,T1,T2,dfreq);
[Ate,Bte] = freeprecess(TE,T1,T2,dfreq);

% 	To add the gradient spoiler twist, we just
%	multiply Atr by zrot(phi):

Atr = zrot(phi)*Atr;

% Let 	M1 be the magnetization just before the tip.
%	M2 be just after the tip.
%	M3 be at TE.
%
% then
%	M2 = Rflip * M1
%	M3 = Ate * M2 + Bte
%	M1 = Atr * M3 + Btr
%
% Solve for M3...
%
%	M3 = Ate*Rflip*Atr*M3 + (Ate*Rflip*Btr+Bte)

Mss = inv(eye(3)-Ate*Rflip*Atr) * (Ate*Rflip*Btr+Bte);
Msig = Mss(1)+i*Mss(2);

end



% 
%	function [Msig,Mss] = gresignal(flip,T1,T2,TE,TR,dfreq)
% 
%	Calculate the steady state gradient-spoiled signal at TE for repeated
%	excitations given T1,T2,TR,TE in ms.  dfreq is the resonant
%	frequency in Hz.  flip is in radians.


function [Msig,Mss] = gresignal(flip,T1,T2,TE,TR,dfreq,num_cycles)


N = 100;
M = zeros(3,N);
phi = ([1:N]/N-0.5 ) * num_cycles*pi;


for k=1:N
	[M1sig,M1] = gssignal(flip,T1,T2,TE,TR,dfreq,phi(k));
	M(:,k)=M1;
end;


Mss = mean(M')';
Msig = Mss(1)+i*Mss(2);
end


%	function [Msig,Mss]=spgrsignal(flip,T1,T2,TE,TR,df,Nex,inc)
%
%	Function calculates the signal from an RF-spoiled sequence
%	after Nex excitations.
%
function [Msig,Mss]=spgrsignal(flip,T1,T2,TE,TR,df,Nex,inc)

if (nargin < 8)
	inc = 117/180*pi;
end;
if (nargin < 7)
	Nex = 100;
end;
if (nargin < 6)
	df = 0;
end;

Nf = 100;	% Simulate 100 different gradient-spoiled spins.
phi = [1:Nf]/Nf*2*pi;

M=zeros(3,Nf,Nex+1);

	
[Ate,Bte] = freeprecess(TE,T1,T2,df);
[Atr,Btr] = freeprecess(TR-TE,T1,T2,df);

M = [zeros(2,Nf);ones(1,Nf)];
on = ones(1,Nf);
	
Rfph = 0;
Rfinc = inc;

for n=1:Nex

	A = Ate * throt(flip,Rfph);
	B = Bte;
	M = A*M+B*on;

	Msig = mean( squeeze(M(1,:)+i*M(2,:)) ) * exp(-i*Rfph);
	Mss = M;

	M=Atr*M+Btr*on;

	for k=1:Nf
		M(:,k) = zrot(phi(k))*M(:,k);
	end;

	Rfph = Rfph+Rfinc;
	Rfinc = Rfinc+inc;
end;
end

% question 3
% Bloch Equation Simulation, Excercise F-1a
% -----------------------------------------
% 
%function [msig,m]=sliceprofile(rf,grad,t,T1,T2,pos,df)
function [msig,m]=sliceprofile(rf,grad,t,T1,T2,pos,df)

gamma = 4258;
dT = t(2)-t(1);         % s.
rfrot = 2*pi*gamma*rf*dT; % Rotation in radians.

pos = pos(:).';		% Make 1xN.
msig = 0*pos;
m = [msig;msig;msig];

for x=1:length(pos)
    tt = sprintf('Simulating postion %d of %d',x,length(pos));
    disp(tt);

    M = [0;0;1];
    [A,B] = freeprecess(1000*dT/2,T1,T2,df);

    for k = 1:length(rf)
	M = A*M+B;
	grot = 2*pi*gamma*(pos(x)/10)*grad(k)*dT/2;
	M = zrot(grot)*M;

    	M = throt(abs(rfrot(k)),angle(rfrot(k))) * M;	% RF Rotation.

	M = A*M+B;
	grot = 2*pi*gamma*(pos(x)/10)*grad(k)*dT/2;
	M = zrot(grot)*M;
    end;
    m(:,x) = M;
    msig(x) = M(1)+i*M(2);

end;
end





