% This function is to calculate the nucleation rate based on classical nucleation theory
% main Ref: S.E.Offerman et al Science 2002; S.E.Offerman et al Scientific
% Report 2016
% Suitable for Fe-C-Mn system
% CNTt: m^-3s^-1

function [CNTt EnergyB Freq Scailing]=CNT_cyclic_Gs(Temp,f_bcc,NpotT,deltaGV,t,dt,xC_F,Comp_m)

% public parameters
h=6.626e-34; % Planck constant [m2.Kg/s] or [J.s]
kB=1.381e-23; % Bolzmann constant [m2.Kg/(s2.K)] or [J/K]
QD=3.93e-19; % activation energy of iron self diffusion [J]

%%%%%%%%%%%%%% K.M. Lee, Philos Mag 2010
% Gs=0.74e7; % % misfit energy of nucleation [J/m3]
% Gs=-0.74e7; % % misfit energy of nucleation [J/m3]

% reference: M.Onink
xC=(Comp_m(1)-f_bcc*xC_F)/(1-f_bcc); % [at.%]
a_bcc=2.886*(1+17.5e-6*(Temp-800)); % [A]
a_fcc=(3.6308+0.0075*xC)*(1+(24.9-0.6*xC)*1e-6*(Temp-1000)); % [A]
Poisson=1/3; % Poisson ratio
mu=6e10; % shear modulus of Fe [N/m2]
% Gs=(2*(1+Poisson)/(9*(1-Poisson)))*mu*((a_bcc^3/2-a_fcc^3/4)/(a_fcc^3/4))^2; % [J/m3] J.D. Eshelby, 1957
Gs=(2*(1-Poisson)/(9*(1+Poisson)))*mu*((a_bcc^3/2-a_fcc^3/4)/(a_fcc^3/4))^2; % [J/m3] Cotes2004_Article_FccHcpMartensiticTransformatio
% Gs=2*mu*(1+Poisson)/(1-Poisson)*((a_bcc^3/2-a_fcc^3/4)/(a_fcc^3/4)).^2; % The strain energy of a coherent ellipsoidal precipitate 1974 [J/m3]
% Gs=6*mu*((a_bcc^3/2-a_fcc^3/4)/(a_fcc^3/4))^2; % The role of strain energy during precipitation of copper and gold from alpha iron 1962 [J/m3]
% Gs=3.8*mu*((a_bcc^3/2-a_fcc^3/4)/(a_fcc^3/4))^2; % precipitate number density in Ni-Al alloy 1970

Z=0.05; % Zeldovich factor, nearly constant
psi=5e-8; % a factor containing the information of shape of critical nucleus and interfacial energies [J3/m6]; 
          % psi_CF=3.3e-3,psi_LEA=2.1e-6,psi_exp1=5e-8,psi_exp2=1~3e-7

% these two parameters should be adjusted to match the experimental data
% by default, tao = 0
SC=1/500; % S/C, scaling factor of potential nucleation sites/nucleation scale, for isothermal and cycling [300-1000]; for cnt, [35 - 3000]
tao=0; % varies in 4~18s according to S.E.Offerman et al Scientific Report 2016
% %% Calculate incubation time according to Aaronson
% a=3.620e-10; % Jumping distance [m]
% sigma_FA=0.62; % interfacial energy between frrite and austenite [J/m^2], H.I. Aaronson Metall.Mater.Trans.A 1988
% V_at=(2.860e-10)^3/2; % average volume of an atom in the nucleus [m^3]
% D_Fe=2.86e-2*exp(-292000/(8.314*Temp)); % diffusivity of iron in fcc [m^2/s]
% tao=8*kB*Temp*a^4*sigma_FA/(V_at^2*deltaGV^2*D_Fe^2); % incubation time [s], the result is wrong

ASC=SC*exp(-tao./t);% pre-factor which could be adjusted by varing SC and tao
% ASC=1/20; % default value, can be adjusted to adjust the nucleation temperature range

if deltaGV<Gs
    Nt=0;
else
    Nt=ASC.*Z*NpotT*(1-f_bcc).*(kB.*Temp./h).*exp(-QD./(kB.*Temp)).*exp(-psi./(kB.*Temp.*(deltaGV-Gs).^2));
end

CNTt=Nt;
EnergyB=exp(-psi./(kB.*Temp.*(deltaGV-Gs).^2));
Freq=(kB.*Temp./h).*exp(-QD./(kB.*Temp));
Scailing=ASC.*Z*NpotT*(1-f_bcc);

