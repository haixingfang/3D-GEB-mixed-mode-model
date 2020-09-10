function ff=solute_drag_fun_soft_impingement(x,Temp,C0_mn,wC_A,wC_F,Ux, ...
    Xneq,Xpeq,X0,DC,kafang,Rbcc,SN,distance,hardflag)

% % only for testing
% % needs to transferred from the main program
% % % load('test_N11_5cycles_solute_InterX.mat');
% Xneq=xC_F_eq(i);
% Xpeq=xC_A_eq(i);
% X0=N_p(l,34);
% DC=D_C(i)*1e-12;
% kafang=Kafang(i);
% Rbcc=N_p(l,7);
% SN=length(N_PR{l}(:,1));
% distance=N_p(l,32);
% hardflag=N_p(l,26);
% Temp=T(i); % [K]
% C0_mn=Comp_m(2);
% wC_A=N_p(l,15);% remote C in austenite [wt.%]
% wC_F=N_p(l,17);% remote C in ferrite [wt.%]
% x=out;

% General parameters
M_Fe=56;
M_C=12;
M_Mn=55;
M_Si=28;
Uwx=1/((1/Ux-1)*M_Mn/M_Fe); % Uwx=w(mn)/w(Fe), Ux=x(Mn)/(x(Mn)+x(Fe))
R=8.314; % [J/(K.mol)]
L_int=0.5e-9; % interface thinckness [m]
delta_int=L_int/2; % half [m]
Nstep=500; % this step should at least > 500
Lint=-3*delta_int:6*delta_int/Nstep:3*delta_int; % define the positions [m], origin indexed by Nstep/2+1
% E0=[9.9 5 3].*1e3; % binding energies of Mn Ni Co [J/mol], same as H.Chen Acta Mater. 2017
% E0=[1 0.5 0.3].*1e3.*R*T; % binding energies dependent on temperature. 2017
E0=6.0e3;% [J/mol] with a uncertainty of 1.4 kJ/mol, experimentally determined by APT, Scprita 2016 by M. Goune group
E0=7e3; % modified on July 17, 2019

% mu(c)=mu0+RTln(xC)+RT(e11xC+e12xMn);
% mu(Mn)=mu0+RTln(xMn)+RT(e12xC+e22xMn);
% mu(Fe)=mu0+RTln(1-xC-xMn)+RT*Yita);
% line1=fcc;line2=bcc
Para_Fe=[28218 -8.44;24312 -8.01];
Para_C=[14547 9.12 -5.66;47969 4.89 -8.11];
Para_Mn=[-49791 -7.63 -1.06;-40813 -6.83 -8.31];

%%%%% convert concentrations
xC_F_int=Xneq/100;
xC_A=100*(wC_A/M_C)/(wC_A/M_C+(100-wC_A)*(Uwx/(1+Uwx))/M_Mn+(100-wC_A)*(1/(1+Ux))/M_Fe); % remote C in austenite [at.%]
xC_F=xC_F_int;
xMn_A=100*(100-xC_A)*Ux/(xC_A+(100-xC_A)*Ux+(100-xC_A)*(1-Ux)); % Mn concentent [mol%] in austenite
xMn_F=100*(100-xC_F)*Ux/(xC_F+(100-xC_F)*Ux+(100-xC_F)*(1-Ux)); % Mn concentent [mol%] in ferrite

%%% Diffusivities [m^2/s]
Dc_alpha=0.02e-4*exp(-10115/Temp)*exp(0.5898*(1+2/pi*atan(14.985-15309/Temp))); % J Agren Acta Metall. 1982 [m^2/s]
y_C=(xC_A/100)/(1-xC_A/100);
Dc_gamma=4.53e-7*((1+y_C*(1-y_C)*8339.9/Temp)*exp(-(1/Temp-2.221e-4)*(17767-26436*y_C)));% J Agren Script Metall. 1986 [m^2/s]
% Dc_alpha=2.2e-4*exp(-125000/(R*T)); %  M. Militzer, M.G. Mecozzi, J. Sietsma, S. van der Zwaag, Acta Mater. 54 (2006)
% Dc_gamma=0.15e-4*exp(-142000/(R*T)); %  M. Militzer, M.G. Mecozzi, J. Sietsma, S. van der Zwaag, Acta Mater. 54 (2006)
Dc_int=sqrt(Dc_alpha*Dc_gamma); % average [m^2/s]

Dmn_alpha=0.756e-4*exp(-224500/(R*Temp)); % H. Oikawa, The Technology Reports of the Tohoku University, 1982, 215
Dmn_gamma=0.178e-4*exp(-264000/(R*Temp)); % H. Oikawa, The Technology Reports of the Tohoku University, 1982, 215
Dmn_int=sqrt(Dmn_alpha*Dmn_gamma); % average [m^2/s]
% Dmn_int=0.5e-4*exp(-247650/(R*T)); % B. Zhu M. Millitzer Comp. Mater. Sci 2015

Dsi_alpha=0.4e-4*exp(-242800/(R*Temp)); % J. Mahieu, Doc Thesis: Ghent Univ. 2004
Dsi_gamma=0.2e-4*exp(-219800/(R*Temp)); % J. Mahieu, Doc Thesis: Ghent Univ. 2004
Dsi_int=sqrt(Dsi_alpha*Dsi_gamma); % average [m^2/s]

Dcr_alpha=2.33e-4*exp(-238800/(R*Temp)); % H. Liu et al Applied Surface Science,2009
Dcr_gamma=0.169e-4*exp(-263900/(R*Temp)); % H. Liu et al Applied Surface Science,2009
Dcr_int=sqrt(Dcr_alpha*Dcr_gamma); % average [m^2/s]

Mumn_gamma=Para_Mn(1,1)+R.*Temp.*log(xMn_A/100)+R.*Temp.*(Para_Mn(1,2).*(xC_A/100)+Para_Mn(1,3).*(xMn_A/100));
Mumn_alpha=Para_Mn(2,1)+R.*Temp.*log(xMn_F/100)+R.*Temp.*(Para_Mn(2,2).*(xC_F/100)+Para_Mn(2,3).*(xMn_F/100));
deltaE_Mn=(Mumn_gamma-Mumn_alpha)/2; % [J/mol]


%% Directly using the integrated results to calculate dissipation energy due to trans-diffusion G_diss
Dim_a=Dmn_int*(deltaE_Mn-E0)/(R*Temp*x(3)*delta_int); % dimensionless parameter a
Dim_b=Dmn_int*(deltaE_Mn+E0)/(R*Temp*x(3)*delta_int); % dimensionless parameter b
Dim_v=abs(x(3)*delta_int/Dmn_int); % dimensionless parameter v

G_diff1=Dim_a.^2*R*Temp*x(3)*C0_mn*delta_int/(Dmn_int*Dim_v*(1+2*Dim_a+Dim_a.^2));
G_diff1=G_diff1*(-exp(Dim_v+Dim_v*Dim_a)+exp(Dim_v+Dim_v*Dim_a)*Dim_v+ ...
    exp(Dim_v+Dim_v*Dim_a)*Dim_v*Dim_a+1)*exp(-Dim_v-Dim_v*Dim_a);

G_diff2=-Dim_b*R*Temp*x(3)*C0_mn*delta_int/(Dmn_int*Dim_v*(1+Dim_a+2*Dim_b+2*Dim_a*Dim_b+ ...
    Dim_b^2+Dim_b^2*Dim_a));
G_diff2=G_diff2*(Dim_a*exp(Dim_v+Dim_v*Dim_b)+Dim_a*Dim_b*exp(Dim_v+Dim_v*Dim_b)+ ...
    Dim_b*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)-Dim_a*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)- ...
    Dim_v*Dim_b*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)+Dim_a*exp(Dim_v+Dim_v*Dim_a)- ...
    Dim_v*Dim_b^2*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)- ...
    Dim_v*Dim_b^2*Dim_a*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)-Dim_a- ...
    Dim_v*Dim_a*Dim_b*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)-Dim_a*Dim_b- ...
    Dim_b*exp(Dim_v+Dim_v*Dim_a))*exp(-2*Dim_v-Dim_v*Dim_b-Dim_v*Dim_a);
G_diff=(G_diff1+G_diff2)/100; % [J/mol]
G_diff_min=(2*C0_mn*deltaE_Mn+R*Temp*C0_mn*(exp(-2*deltaE_Mn/(R*Temp))-1))/100; % [J/mol]

% G_friction: Energy dissipation due to friction, intrinsic interface mobility
M0_int=[2.7e-6 0.035 1.7e-5 4e-7]; % [m3.m/(J.s)]
QM_int=[145e3 147e3 140e3 140e3]; % [J/mol]
% Corresponding to
% Line 1: J. Zhu, H. Chen, Acta 2017;
% Line 2: M. Hillert et al, Scripta 2006;
% Line 3: J.J. Wits et al, Acta 2000;
% Line 4: G.P. Krielaart et al, MSE A 1997;
M_int=M0_int(1)/7.1e-6; % [mol.m/(J.s)]
Mobility=M_int*exp(-QM_int(1)/(R*Temp)); % [mol.m/(J.s)]
G_friction=x(3)./Mobility; % [J/mol]

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% chemical driving force [J/mol]
xMn_A_int=C0_mn./100; % interfacial Mn in ferrite [atomic fraction]
xMn_F_int=C0_mn./100; % interfacial Mn in austenite [atomic fraction]
xFe_A_int=1-x(1)-xMn_A_int; % interfacial Fe in austenite [atomic fraction]
xFe_F_int=1-xC_F_int-xMn_F_int; % interfacial Fe in ferrite [atomic fraction]

Mufe_gamma_int=Para_Fe(1,1)+R.*Temp.*log(xFe_A_int)+R.*Temp.*Para_Fe(1,2);
Muc_gamma_int=Para_C(1,1)+R.*Temp.*log(x(1))+R.*Temp.*(Para_C(1,2).*x(1)+Para_C(1,3).*xMn_A_int);
Mumn_gamma_int=Para_Mn(1,1)+R.*Temp.*log(xMn_A_int)+R.*Temp.*(Para_Mn(1,2).*x(1)+Para_Mn(1,3).*xMn_A_int);

Mufe_alpha_int=Para_Fe(2,1)+R.*Temp.*log(xFe_F_int)+R.*Temp.*Para_Fe(2,2);
Muc_alpha_int=Para_C(2,1)+R.*Temp.*log(xC_F_int)+R.*Temp.*(Para_C(2,2).*xC_F_int+Para_C(2,3).*xMn_F_int);
Mumn_alpha_int=Para_Mn(2,1)+R.*Temp.*log(xMn_F_int)+R.*Temp.*(Para_Mn(2,2).*xC_F_int+Para_Mn(2,3).*xMn_F_int);

% % G_chem=xC_F_int*(Muc_gamma_int-Muc_alpha_int)+xMn_F_int*(Mumn_gamma_int-Mumn_alpha_int)+ ...
% %     xFe_F_int*(Mufe_gamma_int-Mufe_alpha_int); % chemical driving force [J/mol]
% % G_chem=xFe_F_int*(Mufe_gamma_int-Mufe_alpha_int);% chemical driving force [J/mol]
% G_chem=xC_F_int*(Muc_gamma_int-Muc_alpha_int)+ ...
%     xFe_F_int*(Mufe_gamma_int-Mufe_alpha_int); % chemical driving force [J/mol], more accurate
G_chem=kafang*(Xpeq-x(1)); % simplified way to calculate driving force

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vbcc=4/3*pi*Rbcc^3; % [um^3]
Z=2*DC/x(3)*1e6; % [um]

%%%%%%%%%%%%%%%%% with soft impingement
ff(1)=30*Vbcc*(X0-Xneq)/(4*pi)*1/SN-(9*(distance-Rbcc)^3*x(2)-10*(distance-Rbcc)^3*X0+ ...
            (distance-Rbcc)^3*x(1)+20*Rbcc^2*x(2)*(distance-Rbcc)- ...
        30*Rbcc^2*X0*(distance-Rbcc)+10*Rbcc^2*x(1)*(distance-Rbcc)+25*(distance-Rbcc)^2*Rbcc*x(2)- ...
        30*(distance-Rbcc)^2*Rbcc*X0+5*(distance-Rbcc)^2*Rbcc*x(1))*1/SN;
ff(2)=(x(1)-Xneq)*(distance-Rbcc)-Z*(x(1)-x(2));
ff(3)=G_chem-G_friction-G_diff;



