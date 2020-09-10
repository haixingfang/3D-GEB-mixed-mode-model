% Oct 13, 2019
% update July 8, 2020
function [Drag Diffusion GibbsE]=loopV_solute_drag(V_previous,Temp,C0_mn,wC_A,wC_F,Ux, ...
    Xneq,Xpeq,X0,DC,kafang,Rbcc,SN,distance,hardflag,f_alpha,CyclicFlag,CycleStart,tt,TransDir)

% i=54
% clear G_chem V_int G_diff G_friction C_diff;
% % only for testing
% % needs to transferred from the main program
% % % load('test_N11_5cycles_solute_InterX.mat');
% Xneq=xC_F_eq(i);
% Xpeq=xC_A_eq(i);
% X0=Nucleated{i-1}(ll,34);
% DC=D_C(i)*1e-12;
% kafang=Kafang(i);
% % Rbcc=N_p(l,7);
% Rbcc=Nucleated{i-1}(ll,7);
% SN=length(N_PR{l}(:,1));
% distance=Nucleated{i-1}(ll,32);
% hardflag=Nucleated{i-1}(ll,26);
% Temp=T(i); % [K]
% C0_mn=Comp_m(2);
% wC_A=Nucleated{i-1}(ll,15);% remote C in austenite [wt.%]
% wC_F=Nucleated{i-1}(ll,17);% remote C in ferrite [wt.%]
% CycleStart=tcr(3);
% tt=Timer(i);
% V_previous=V0(DragIter);

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
E0=8e3;
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
y_C=xC_A/(1-xC_A);
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

%%%% the chemical potential needs to be reconsidered May 24 2018
Mumn_gamma=Para_Mn(1,1)+R.*Temp.*log(xMn_A/100)+R.*Temp.*(Para_Mn(1,2).*(xC_A/100)+Para_Mn(1,3).*(xMn_A/100));
Mumn_alpha=Para_Mn(2,1)+R.*Temp.*log(xMn_F/100)+R.*Temp.*(Para_Mn(2,2).*(xC_F/100)+Para_Mn(2,3).*(xMn_F/100));
deltaE_Mn=(Mumn_gamma-Mumn_alpha)/2; % [J/mol]
G_diff_min=(2*C0_mn*deltaE_Mn+R*Temp*C0_mn*(exp(-2*deltaE_Mn/(R*Temp))-1))/100; % [J/mol]

% calculate the potential profile for X
if CyclicFlag==1 && tt>=CycleStart
    V_int=[-logspace(log10(1e-10),log10(5e-6),20) logspace(log10(1e-10),log10(5e-6),20)]; % interface velocity [m/s]
else
    V_int=logspace(log10(1e-10),log10(5e-6),20); % interface velocity [m/s]
end
for i=1:length(V_int)
    G_diff(i)=0;
    for j=1:length(Lint)
        Dim_a(i,j)=Dmn_int*(deltaE_Mn-E0)/(R*Temp*V_int(i)*delta_int); % dimensionless parameter a
        Dim_b(i,j)=Dmn_int*(deltaE_Mn+E0)/(R*Temp*V_int(i)*delta_int); % dimensionless parameter b
        Dim_v(i,j)=V_int(i)*delta_int/Dmn_int; % dimensionless parameter v
        Dim_v(i,j)=abs(Dim_v(i,j)); % add on April 24,2020
        Dim_x(i,j)=Lint(j)/delta_int; % dimensionless variable x
       if Lint(j)<=-delta_int
           muMn_int(i,j)=Mumn_alpha; % [J/mol]
           C_mn(i,j)=C0_mn; % [mol%]
       else if Lint(j)>-delta_int && Lint(j)<=0
               muMn_int(i,j)=Mumn_alpha+deltaE_Mn-E0+(deltaE_Mn-E0)/delta_int.*Lint(j);
               C_mn(i,j)=C0_mn*(1+Dim_a(i,j)*exp(-Dim_v(i,j)*(1+Dim_a(i,j))*(Dim_x(i,j)+1)))/(1+Dim_a(i,j));
               G_diff(i)=G_diff(i)+(C_mn(i,j)-C0_mn)*(Mumn_alpha-(Mumn_alpha+deltaE_Mn-E0))/(-delta_int)*(Lint(2)-Lint(1));        
           else if Lint(j)>0 && Lint(j)<=delta_int
                   muMn_int(i,j)=Mumn_alpha+deltaE_Mn-E0+(deltaE_Mn+E0)/delta_int.*Lint(j);
                   C_mn(i,j)=C0_mn*(1+(Dim_a(i,j)*(1+Dim_b(i,j))*exp(-Dim_v(i,j)*(1+Dim_a(i,j)))/(1+Dim_a(i,j))+ ...
                       (Dim_b(i,j)-Dim_a(i,j))/(1+Dim_a(i,j)))*exp(-Dim_v(i,j)*(1+Dim_b(i,j))*Dim_x(i,j))) ...
                       /(1+Dim_b(i,j));
                   G_diff(i)=G_diff(i)+(C_mn(i,j)-C0_mn)*(Mumn_gamma-(Mumn_alpha+deltaE_Mn-E0))/(delta_int)*(Lint(2)-Lint(1));
               else if Lint(j)>delta_int
                       muMn_int(i,j)=Mumn_gamma;
                       C_mn(i,j)=C0_mn*(1+exp(-Dim_v(i,j)*Dim_x(i,j))*((Dim_a(i,j)*exp(-Dim_v(i,j)*(Dim_a(i,j)+Dim_b(i,j)+ ...
                           1))/(1+Dim_a(i,j)))+(Dim_b(i,j)-Dim_a(i,j))*exp(-Dim_v(i,j)*Dim_b(i,j))/((1+Dim_a(i,j))*(1+Dim_b(i,j)))- ...
                           Dim_b(i,j)*exp(Dim_v(i,j))/(1+Dim_b(i,j))));
                   end
               end
           end
       end
    end
    G_diff1=Dim_a(i,j).^2*R*Temp*V_int(i)*C0_mn*delta_int/(Dmn_int*Dim_v(i,j)*(1+2*Dim_a(i,j)+Dim_a(i,j).^2));
    G_diff1=G_diff1*(-exp(Dim_v(i,j)+Dim_v(i,j)*Dim_a(i,j))+exp(Dim_v(i,j)+Dim_v(i,j)*Dim_a(i,j))*Dim_v(i,j)+ ...
        exp(Dim_v(i,j)+Dim_v(i,j)*Dim_a(i,j))*Dim_v(i,j)*Dim_a(i,j)+1)*exp(-Dim_v(i,j)-Dim_v(i,j)*Dim_a(i,j));
    G_diff2=-Dim_b(i,j)*R*Temp*V_int(i)*C0_mn*delta_int/(Dmn_int*Dim_v(i,j)*(1+Dim_a(i,j)+2*Dim_b(i,j)+2*Dim_a(i,j)*Dim_b(i,j)+ ...
        Dim_b(i,j)^2+Dim_b(i,j)^2*Dim_a(i,j)));
    G_diff2=G_diff2*(Dim_a(i,j)*exp(Dim_v(i,j)+Dim_v(i,j)*Dim_b(i,j))+Dim_a(i,j)*Dim_b(i,j)*exp(Dim_v(i,j)+Dim_v(i,j)*Dim_b(i,j))+ ...
    Dim_b(i,j)*exp(2*Dim_v(i,j)+Dim_v(i,j)*Dim_b(i,j)+Dim_v(i,j)*Dim_a(i,j))-Dim_a(i,j)*exp(2*Dim_v(i,j)+Dim_v(i,j)*Dim_b(i,j)+Dim_v(i,j)*Dim_a(i,j))- ...
    Dim_v(i,j)*Dim_b(i,j)*exp(2*Dim_v(i,j)+Dim_v(i,j)*Dim_b(i,j)+Dim_v(i,j)*Dim_a(i,j))+Dim_a(i,j)*exp(Dim_v(i,j)+Dim_v(i,j)*Dim_a(i,j))- ...
    Dim_v(i,j)*Dim_b(i,j)^2*exp(2*Dim_v(i,j)+Dim_v(i,j)*Dim_b(i,j)+Dim_v(i,j)*Dim_a(i,j))- ...
    Dim_v(i,j)*Dim_b(i,j)^2*Dim_a(i,j)*exp(2*Dim_v(i,j)+Dim_v(i,j)*Dim_b(i,j)+Dim_v(i,j)*Dim_a(i,j))-Dim_a(i,j)- ...
    Dim_v(i,j)*Dim_a(i,j)*Dim_b(i,j)*exp(2*Dim_v(i,j)+Dim_v(i,j)*Dim_b(i,j)+Dim_v(i,j)*Dim_a(i,j))-Dim_a(i,j)*Dim_b(i,j)- ...
    Dim_b(i,j)*exp(Dim_v(i,j)+Dim_v(i,j)*Dim_a(i,j)))*exp(-2*Dim_v(i,j)-Dim_v(i,j)*Dim_b(i,j)-Dim_v(i,j)*Dim_a(i,j));
    G_diff_analytical(i)=(G_diff1+G_diff2)/100; % [J/mol]
    G_diff_min=(2*C0_mn*deltaE_Mn+R*Temp*C0_mn*(exp(-2*deltaE_Mn/(R*Temp))-1))/100; % [J/mol]
    %%%% Energy dissipation due to X diffusion inside interface
    G_diff(i)=-G_diff(i)/100; % convert to the right unit [J/mol]

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
    G_friction(i)=V_int(i)./Mobility; % [J/mol]

    %%%
    [C_diff(i,:)]=Mixmode_diffusion_profile_vpasolve_drag(Xneq, ...
    Xpeq,X0,DC,V_int(i),Rbcc,SN,distance,hardflag);% mix-mode update in July 2018
    C_diff(i,:)=abs(C_diff(i,:));
    xC_A_int(i)=C_diff(i,2)/100;

    %%%% chemical driving force [J/mol]
    % xMn_A_int(i)=C_mn(i,2*(fix(2*delta_int/(6*delta_int/Nstep))+1))./100; % interfacial Mn in austenite [atomic fraction]
    xMn_A_int(i)=C0_mn./100; % interfacial Mn in ferrite [atomic fraction]
    xMn_F_int(i)=C_mn(i,fix(2*delta_int/(6*delta_int/Nstep))+1)./100;% interfacial Mn in ferrite [atomic fraction]
    % xMn_F_int(i)=C0_mn./100;
    xFe_A_int(i)=1-xC_A_int(i)-xMn_A_int(i); % interfacial Fe in austenite [atomic fraction]
    xFe_F_int(i)=1-xC_F_int-xMn_F_int(i); % interfacial Fe in ferrite [atomic fraction]

    Mufe_gamma_int(i)=Para_Fe(1,1)+R.*Temp.*log(xFe_A_int(i))+R.*Temp.*Para_Fe(1,2);
    Muc_gamma_int(i)=Para_C(1,1)+R.*Temp.*log(xC_A_int(i))+R.*Temp.*(Para_C(1,2).*xC_A_int(i)+Para_C(1,3).*xMn_A_int(i));
    Mumn_gamma_int(i)=Para_Mn(1,1)+R.*Temp.*log(xMn_A_int(i))+R.*Temp.*(Para_Mn(1,2).*xC_A_int(i)+Para_Mn(1,3).*xMn_A_int(i));

    Mufe_alpha_int(i)=Para_Fe(2,1)+R.*Temp.*log(xFe_F_int(i))+R.*Temp.*Para_Fe(2,2);
    Muc_alpha_int(i)=Para_C(2,1)+R.*Temp.*log(xC_F_int)+R.*Temp.*(Para_C(2,2).*xC_F_int+Para_C(2,3).*xMn_F_int(i));
    Mumn_alpha_int(i)=Para_Mn(2,1)+R.*Temp.*log(xMn_F_int(i))+R.*Temp.*(Para_Mn(2,2).*xC_F_int+Para_Mn(2,3).*xMn_F_int(i));

%     G_chem(i)=xC_A_int(i)*(Muc_gamma_int(i)-Muc_alpha_int(i))+xMn_A_int(i)*(Mumn_gamma_int(i)-Mumn_alpha_int(i))+ ...
%         xFe_A_int(i)*(Mufe_gamma_int(i)-Mufe_alpha_int(i)); % chemical driving force [J/mol]
%     G_chem(i)=xC_F_int*(Muc_gamma_int(i)-Muc_alpha_int(i))+xMn_F_int(i)*(Mumn_gamma_int(i)-Mumn_alpha_int(i))+ ...
%         xFe_F_int(i)*(Mufe_gamma_int(i)-Mufe_alpha_int(i)); % chemical driving force [J/mol]
%     G_chem(i)=xC_F_int*(Muc_gamma_int(i)-Muc_alpha_int(i))+xMn_F_int(i)*(Mumn_gamma_int(i)-Mumn_alpha_int(i));
%     G_chem(i)=xFe_F_int(i)*(Mufe_gamma_int(i)-Mufe_alpha_int(i));% chemical driving force [J/mol]
    G_chem(i)=xC_F_int*(Muc_gamma_int(i)-Muc_alpha_int(i))+ ...
        xFe_F_int(i)*(Mufe_gamma_int(i)-Mufe_alpha_int(i)); % chemical driving force [J/mol], more accurate
%     if i>2 && abs(G_chem(i))>abs(G_chem(i-1)) && V_int(i)>0
%         G_chem(i)=G_chem(i-1)-(G_chem(i-2)-G_chem(i-1));
%     elseif i>2 && abs(G_chem(i))>abs(G_chem(i-1)) && V_int(i)<0
%         G_chem(i)=G_chem(i-1)+abs(G_chem(i-2)-G_chem(i-1));
%     end
%     if C_diff(i,2)==C_diff(i,1) && C_diff(i,2)==C_diff(i,6)
%         G_chem(i)=0;
%     end
    G_chem_kafang(i)=kafang*(Xpeq-xC_A_int(i)*100); % simplified way to calculate driving force
    G_chem(i)=kafang*(Xpeq-xC_A_int(i)*100); % simplified way to calculate driving force
    i;
end
if all(V_int>0)
    Vsol = InterX([V_int;G_chem],[V_int;G_diff+G_friction]); % [m/s, J/mol]
else
    Vsol1 = InterX([V_int(1:length(V_int)/2);G_chem(1:length(V_int)/2)], ...
        [V_int(1:length(V_int)/2);G_diff(1:length(V_int)/2)+G_friction(1:length(V_int)/2)]); % [m/s, J/mol]
    Vsol2 = InterX([V_int(length(V_int)/2+1:end);G_chem(length(V_int)/2+1:end)], ...
        [V_int(length(V_int)/2+1:end);G_diff(length(V_int)/2+1:end)+G_friction(length(V_int)/2+1:end)]); % [m/s, J/mol]
    Vsol=[Vsol1 Vsol2];
end
if ~isempty(Vsol)
    Velocity=Vsol(1,find(abs(Vsol(1,:))==min(abs(Vsol(1,:)))));
    Gchem=Vsol(2,find(abs(Vsol(1,:))==min(abs(Vsol(1,:)))));
    Gdiff = interp1(V_int,G_diff,Velocity);
    Gfriction = interp1(V_int,G_friction,Velocity);
    Gdiss = Gdiff+Gfriction;
    exitflag=1;
    Xpm=interp1(V_int,C_diff(:,1),Velocity);
    if Xpm<X0
        Xpm=X0;
    end
    Xip=interp1(V_int,C_diff(:,2),Velocity);
    DiffLL=interp1(V_int,C_diff(:,3),Velocity);
    Rbcc=interp1(V_int,C_diff(:,4),Velocity);
%     softflag=round(interp1(V_int,C_diff(:,7),Velocity));
    if Xpm>X0
        softflag=1;
    else
        softflag=0;
    end
else
    if TransDir>=0
        Index_sol=find(V_int==min(abs(V_int)));
        Velocity=0;
    else
        Index_sol=find(V_int==max(-abs(V_int)));
        Velocity=0;
    end
    Gchem=G_chem(Index_sol);
    Gdiff=G_diff(Index_sol);
    Gfriction=Velocity./Mobility; % [J/mol]
    Gdiss=Gdiff+Gfriction;
    if abs(Gchem)<abs(Gdiss)
        exitflag=1;
    else
        exitflag=0;
    end
    Xpm=C_diff(Index_sol,1);
    Xip=C_diff(Index_sol,2);
    DiffLL=C_diff(Index_sol,3);
    softflag=C_diff(Index_sol,7);
end
Drag=[Velocity Gchem Gfriction Gdiff Gdiss exitflag]';
Diffusion=[Xpm Xip DiffLL Rbcc Xneq Xpeq softflag]';
GibbsE=[V_int;G_chem;G_diff;G_friction;G_diff+G_friction]';

% plot
plotflag=0;
if plotflag==1
if CyclicFlag==1 && tt>=CycleStart
    figure;
    as=1;
    bs=fix(length(V_int)/2);
    subplot(1,2,1);
    hold all;
    plot(-V_int(as:bs),G_diff(as:bs),'ro-');
    plot(-V_int(as:bs),G_friction(as:bs),'bx-');
    plot(-V_int(as:bs),G_diff(as:bs)+G_friction(as:bs),'kd-');
    plot(-V_int(as:bs),G_chem(as:bs),'m^-');
    legend('G_{diff}','G_{friction}','G_{total}','G_{chem}');
    xlabel('v_{int} (m/s)');
    ylabel('Gibbs energy (J/mol)');
    set(gca,'XScale','log');
    box on;set(gca,'fontsize',12);
    set(gca,'linewidth',1.5);
    title('(a) V_{int} < 0');
    
    as=fix(length(V_int)/2)+1;
    bs=length(V_int);
    subplot(1,2,2);
    hold all;
    plot(V_int(as:bs),G_diff(as:bs),'ro-');
    plot(V_int(as:bs),G_friction(as:bs),'bx-');
    plot(V_int(as:bs),G_diff(as:bs)+G_friction(as:bs),'kd-');
    plot(V_int(as:bs),G_chem(as:bs),'m^-');
    legend('G_{diff}','G_{friction}','G_{total}','G_{chem}');
    xlabel('v_{int} (m/s)');
    ylabel('Gibbs energy (J/mol)');
    set(gca,'XScale','log');
    box on;set(gca,'fontsize',12);
    set(gca,'linewidth',1.5);
    title('b) V_{int} > 0');
else
    figure;
    hold all;
    plot(V_int,G_diff,'ro-');
    plot(V_int,G_friction,'bx-');
    plot(V_int,G_diff+G_friction,'kd-');
    plot(V_int,G_chem,'m^-');
    legend('G_{diff}','G_{friction}','G_{total}','G_{chem}');
    xlabel('v_{int} (m/s)');
    ylabel('Gibbs energy (J/mol)');
    set(gca,'XScale','log');
    box on;set(gca,'fontsize',12);
    set(gca,'linewidth',1.5);
end
end
