% BME 599, HW2
% Problem 2

load('brain_maps.mat')

subplot(1,3,1);
imshow(M0map,[]);
title("M0 map")
colorbar
colormap jet

subplot(1,3,2);
imshow(T1map,[]);
title("T1 map")
colorbar
colormap jet

subplot(1,3,3);
imshow(T2map,[]);
title("T2 map")
colorbar
colormap jet

%% 2a

TR_vec_ms = [500, 4000, 6000]; % from lecture slides
TE_vec_ms = [15, 15, 100]; % from lecture slides
T1w_im = zeros(256,256); % short TR, short TE
PDW_im = zeros(256,256); % long TR, short TE
T2w_im = zeros(256,256); % long TR, long TE

warning off
for j = 1:256
    for k = 1:256
        [Msig,~] = sesignal(T1map(j,k),T2map(j,k),TE_vec_ms(1),TR_vec_ms(1),0,M0map(j,k));
        T1w_im(j,k) = abs(Msig);
        [Msig,~] = sesignal(T1map(j,k),T2map(j,k),TE_vec_ms(2),TR_vec_ms(2),0,M0map(j,k));
        PDW_im(j,k) = abs(Msig);
        [Msig,~] = sesignal(T1map(j,k),T2map(j,k),TE_vec_ms(3),TR_vec_ms(3),0,M0map(j,k));
        T2w_im(j,k) = abs(Msig);
    end

end

%% 2a plots
subplot(1,3,1);
imshow(PDW_im,[]);
title("PDW")

subplot(1,3,2);
imshow(T1w_im,[]);
title("T1w")

subplot(1,3,3);
imshow(T2w_im,[]);
title("T2w")

%% 2b

T1_ms_vec = [1000,1000,2000,2000];
T2_ms_vec = [50,100,50,100];
TR_ms = 3000;
ES_ms = 5;
ETL = 32;

for i = 1:4
    [Msig,Mss] = fsesignal(T1_ms_vec(i),T2_ms_vec(i),ES_ms,TR_ms,0,ETL); % size 32 (values for each echo)
end


%% subfunctions explicitly from http://mrsrl.stanford.edu/~brian/bloch/
% made some edits to solve the homework questions, but base code from here

% 
%	function [Msig,Mss] = sesignal(T1,T2,TE,TR,dfreq,m0)
% 
%	Calculate the steady state signal at TE for a spin-echo
%	sequence, given T1,T2,TR,TE in ms.  Force the
%	transverse magnetization to zero before each excitation.
%	dfreq is the resonant frequency in Hz.  flip is in radians.
%

function [Msig,Mss] = sesignal(T1,T2,TE,TR,dfreq,m0)

Rflip = yrot(pi/2);	% Rotation from excitation pulse (90)
Rrefoc = xrot(pi);	% Rotation from refocusing pulse (usually 180)

[Atr,Btr] = freeprecess(TR-TE,T1,T2,dfreq);	% Propagation TE to TR
[Ate2,Bte2] = freeprecess(TE/2,T1,T2,dfreq);	% Propagation 0 to TE/2
						% (same as TE/2 to TE)

% Neglect residual transverse magnetization prior to excitation.
%Atr = [0 0 0;0 0 0;0 0 1]*Atr;		% (Just keep Mz component)
Atr = [0 0 0;0 0 0;0 0 m0]*Atr;		% (Just keep Mz component) %%% changed from original version for the homework


% Let 	M1 be the magnetization just before the 90.
%	M2 be just before the 180.
%	M3 be at TE.
%	M4 = M1
%
% then
%	M2 = Ate2*Rflip*M1 + Bte2
%	M3 = Ate2*Rrefoc*M2 + Bte2
%	M4 = Atr * M3 + Btr
%
% Solve for M3... (Magnetization at TE)
%

Mss = inv(eye(3)-Ate2*Rrefoc*Ate2*Rflip*Atr) * (Bte2+Ate2*Rrefoc*(Bte2+Ate2*Rflip*Btr));

Msig = Mss(1)+i*Mss(2);

end

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
%	function [Msig,Mss] = fsesignal(T1,T2,TE,TR,dfreq,ETL)
% 
%	Calculate the steady state signal at TE for a multi-echo spin-echo
%	sequence, given T1,T2,TR,TE in ms.  Force the
%	transverse magnetization to zero before each excitation.
%	dfreq is the resonant frequency in Hz.  flip is in radians.
%

function [Msig,Mss] = fsesignal(T1,T2,TE,TR,dfreq,ETL)

Rflip = yrot(pi/2);	% Rotation from Excitation  (usually 90)
Rrefoc = xrot(pi);	% Rotation from Refocusing (usually 180)

[Atr,Btr] = freeprecess(TR-ETL*TE,T1,T2,dfreq);	% Propagation last echo to TR
[Ate2,Bte2] = freeprecess(TE/2,T1,T2,dfreq);	% Propagation over TE/2

% Neglect residual transverse magnetization prior to excitation.
Atr = [0 0 0;0 0 0;0 0 1]*Atr;	% Retain only Mz component.


% Since ETL varies, let's keep a "running" A and B.  We'll
% calculate the steady-state signal just after the tip, Rflip.

% Initial.
A=eye(3);
B=[0 0 0]';


% For each echo, we "propagate" A and B:
for k=1:ETL
	A = Ate2*Rrefoc*Ate2*A;			% TE/2 -- Refoc -- TE/2
	B = Bte2+Ate2*Rrefoc*(Bte2+Ate2*B);
end;


% Propagate A and B through to just after flip, and calculate steady-state.
A = Rflip*Atr*A;
B = Rflip*(Btr+Atr*B);

Mss = inv(eye(3)-A)*B;	% -- Steady state is right after 90 pulse!
M = Mss;


% Calculate signal on each echo.
for k=1:ETL
	M = Ate2*Rrefoc*Ate2*M + Bte2+Ate2*Rrefoc*Bte2;
	Mss(:,k)=M;
	Msig(k)=M(1)+i*M(2);
end;
end