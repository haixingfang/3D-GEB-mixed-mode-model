function [Xip G_chem G_chem_kafang G_diff]=loop_Xip(Temp,C0_mn,wC_A,wC_F,Ux,Xneq,Xpeq,X0,kafang)

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

Para_Fe=[28218 -8.44;24312 -8.01];
Para_C=[14547 9.12 -5.66;47969 4.89 -8.11];
Para_Mn=[-49791 -7.63 -1.06;-40813 -6.83 -8.31];

Dmn_alpha=0.756e-4*exp(-224500/(R*Temp)); % H. Oikawa, The Technology Reports of the Tohoku University, 1982, 215
Dmn_gamma=0.178e-4*exp(-264000/(R*Temp)); % H. Oikawa, The Technology Reports of the Tohoku University, 1982, 215
Dmn_int=sqrt(Dmn_alpha*Dmn_gamma); % average [m^2/s]

V_int=1e-7; % interface velocity [m/s]
xC_A_int=linspace(X0,Xpeq,10)/100;
% xC_A_int=linspace(X0,2*Xpeq,20)/100;
Xip=xC_A_int*100;
for i=1:length(xC_A_int)   
    %%%%% convert concentrations
    xC_F_int=Xneq/100;
    xC_A=100*(wC_A/M_C)/(wC_A/M_C+(100-wC_A)*(Uwx/(1+Uwx))/M_Mn+(100-wC_A)*(1/(1+Ux))/M_Fe); % remote C in austenite [at.%]
    xC_F=xC_F_int;
    xMn_A=100*(100-xC_A)*Ux/(xC_A+(100-xC_A)*Ux+(100-xC_A)*(1-Ux)); % Mn concentent [mol%] in austenite
    xMn_F=100*(100-xC_F)*Ux/(xC_F+(100-xC_F)*Ux+(100-xC_F)*(1-Ux)); % Mn concentent [mol%] in ferrite

    Mumn_gamma=Para_Mn(1,1)+R.*Temp.*log(xMn_A/100)+R.*Temp.*(Para_Mn(1,2).*(xC_A/100)+Para_Mn(1,3).*(xMn_A/100));
    Mumn_alpha=Para_Mn(2,1)+R.*Temp.*log(xMn_F/100)+R.*Temp.*(Para_Mn(2,2).*(xC_F/100)+Para_Mn(2,3).*(xMn_F/100));
    deltaE_Mn=(Mumn_gamma-Mumn_alpha)/2; % [J/mol]
    G_diff(i)=0;
    for j=1:length(Lint)
        Dim_a(i,j)=Dmn_int*(deltaE_Mn-E0)/(R*Temp*V_int*delta_int); % dimensionless parameter a
        Dim_b(i,j)=Dmn_int*(deltaE_Mn+E0)/(R*Temp*V_int*delta_int); % dimensionless parameter b
        Dim_v(i,j)=V_int*delta_int/Dmn_int; % dimensionless parameter v
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
    
    %%%% Energy dissipation due to X diffusion inside interface
    G_diff(i)=-G_diff(i)/100; % convert to the right unit [J/mol]
    if V_int>0
        G_diff(i)=abs(G_diff(i));
    else
        G_diff(i)=-abs(G_diff(i));
    end

    %%%% chemical driving force [J/mol]
    xMn_A_int(i)=C_mn(i,2*(fix(2*delta_int/(6*delta_int/Nstep))+1))./100; % interfacial Mn in austenite [atomic fraction]
%     xMn_A_int(i)=C0_mn./100; % interfacial Mn in ferrite [atomic fraction]
    xMn_F_int(i)=C_mn(i,fix(2*delta_int/(6*delta_int/Nstep))+1)./100;% interfacial Mn in ferrite [atomic fraction]
%     xMn_F_int(i)=C0_mn./100;
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
%     G_chem(i)=xFe_F_int(i)*(Mufe_gamma_int(i)-Mufe_alpha_int(i));% chemical driving force [J/mol]
    G_chem(i)=xC_F_int*(Muc_gamma_int(i)-Muc_alpha_int(i))+ ...
        xFe_F_int(i)*(Mufe_gamma_int(i)-Mufe_alpha_int(i)); % chemical driving force [J/mol], more accurate
    G_chem_kafang(i)=kafang*(Xpeq-xC_A_int(i)*100); % simplified way to calculate driving force
end
