function [Drag Diffusion Gdiff0]=solute_drag_fsolver(V_previous,Temp,C0_mn,wC_A,wC_F,Ux, ...
    Xneq,Xpeq,X0,DC,kafang,Rbcc,SN,distance,hardflag,f_alpha,Xip0,CyclicFlag,CycleStart,tt,TransDir)

% % only for testing
% % needs to transferred from the main program
% % % load('test_N11_5cycles_solute_InterX.mat');
% i=130
% Xneq=xC_F_eq(i);
% Xpeq=xC_A_eq(i);
% X0=Nucleated{i-1}(ll,34);
% DC=D_C(i)*1e-12;
% kafang=Kafang(i);
% Rbcc=Nucleated{i-1}(ll,7);
% SN=length(N_PR{l}(:,1));
% distance=Nucleated{i-1}(ll,32);
% hardflag=Nucleated{i-1}(ll,26);
% Temp=T(i); % [K]
% C0_mn=Comp_m(2);
% wC_A=N_p(l,15);% remote C in austenite [wt.%]
% wC_F=N_p(l,17);% remote C in ferrite [wt.%]
% V_previous=v_t(i-1,l)/1e6; % [m/s]
% if V_previous==0 || abs(V_previous)>1e-6
%   V_previous=5e-8; % [m/s]
% end
% clear V_int;
% clear Dim_a Dim_b Dim_v G_diff1 G_diff2 G_diff G_friction;
% clear X1 X2 X3;
% clear G_chem xFe_A_int Mufe_gamma_int Muc_gamma_int;
% clear Z;
% Xip0=DiffInfo{l}(i-1,1,2);
% CycleStart=tcr(3);
% tt=Timer(i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
E0=7.0e3; % modified on July 17, 2019

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
G_diff_min=(2*C0_mn*deltaE_Mn+R*Temp*C0_mn*(exp(-2*deltaE_Mn/(R*Temp))-1))/100; % [J/mol]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ff=@(x)solute_drag_fun(x,Temp,C0_mn,wC_A,wC_F,Ux,Xneq,Xpeq,X0,DC,kafang,Rbcc,SN,distance,hardflag);
% opts=optimset('Display','iter','Algorithm',{'levenberg-marquardt',0.005},'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-5);% option settings for fsolve
% opts=optimset('Display','notify','Algorithm',{'levenberg-marquardt',0.005},'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-5);% option settings for fsolve
opts=optimset('Display','off','Algorithm',{'levenberg-marquardt',0.005},'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-10);% option settings for fsolve
SolFlag=0;
v0=[V_previous 5e-7 1e-8 1e-9];
jj=1;
while SolFlag~=1
    x0=[X0+(Xpeq-X0).*rand(1,1) distance.*rand(1,1) v0(jj)];
    [out,fval,exitflag]=fsolve(ff,x0,opts); % [Xip DiffLL Velocity]
    if jj==length(v0)
        Xip=[];
        DiffLL=[];
        Velocity=[];
        Xpm=X0;
        SolFlag=1;
    end
    if out(1)>X0 && (out(2)>0 && out(2)<distance) && exitflag>0
        Xip=out(1);
        DiffLL=out(2);
        Velocity=out(3);
        Xpm=X0;
        SolFlag=1;
    end
    if (out(1)>X0 && out(1)<Xpeq) && (out(2)>0 && out(2)<distance) && exitflag<=0 && Rbcc<0.5
        Xip=out(1);
        DiffLL=out(2);
        Velocity=out(3);
        Xpm=X0;
        SolFlag=1;
    end
    jj=jj+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(DiffLL)
    if DiffLL<=(distance-Rbcc) %&& hardflag==0
        softflag=0;
    else
        ff=@(x)solute_drag_fun_soft_impingement(x,Temp,C0_mn,wC_A,wC_F,Ux,Xneq,Xpeq,X0,DC,kafang,Rbcc,SN,distance,hardflag);
%         opts=optimset('Display','iter','Algorithm',{'levenberg-marquardt',0.005},'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-6);% option settings for fsolve
        opts=optimset('Display','off','Algorithm',{'levenberg-marquardt',0.005},'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-6);% option settings for fsolve
        SolFlag=0;
        v0=[V_previous 1e-9 1e-8 1e-7 5e-7];
        jj=1;
        temp_out=[];
        while SolFlag~=1
            x0=[X0+(Xpeq-X0).*rand(1,1) X0 v0(jj)];
            [out,fval,exitflag]=fsolve(ff,x0,opts); % [Xip Xpm Velocity]
            if jj==length(v0)
                Xip=[];
                Xpm=[];
                Velocity=[];
                DiffLL=distance-Rbcc;
                SolFlag=1;
            end
            if out(1)>X0 && (out(2)>0 && out(2)<0.8*4.6) && exitflag>0 && ~(out(1)>Xpeq && out(2)<out(1)) && ...
                ~(out(3)<0 && (CyclicFlag<1 || (CyclicFlag==1 && tt<CycleStart)))
                Xip=out(1);
                Xpm=out(2);
                Velocity=out(3);
                DiffLL=distance-Rbcc;
                SolFlag=1;
            end
            temp_out=[temp_out;out];
            jj=jj+1;
        end
        softflag=1;
        if isempty(Xip)
            if Rbcc<distance
                Xip=(X0-Rbcc^3/distance^3.*Xneq)./(1-Rbcc^3/distance^3);
                if Xip>Xpeq && TransDir==1 && tt<CycleStart
                    Xip=Xpeq;
                end
            else
                Xip=Xpeq;
            end
            DiffLL=0;
            Xpm=Xip;
            softflag=1;
            Gchem=kafang*(Xpeq-Xip); % [J/mol]
            if Gchem<0
                G_diff_min=minGdiff(Dmn_int,deltaE_Mn,E0,R,Temp,delta_int,C0_mn,-1e-10);
            else
                G_diff_min=minGdiff(Dmn_int,deltaE_Mn,E0,R,Temp,delta_int,C0_mn,1e-10);
            end
            if abs(Gchem)>abs(G_diff_min)
                syms V_int; % [m/s]
                syms Dim_a Dim_b Dim_v G_diff1 G_diff2 G_diff G_friction;
                %% Directly using the integrated results to calculate dissipation energy due to trans-diffusion G_diss
                Dim_a=Dmn_int*(deltaE_Mn-E0)/(R*Temp*V_int*delta_int); % dimensionless parameter a
                Dim_b=Dmn_int*(deltaE_Mn+E0)/(R*Temp*V_int*delta_int); % dimensionless parameter b
                Dim_v=abs(V_int*delta_int/Dmn_int); % dimensionless parameter v
                Dim_v=abs(Dim_v); % add on April 24,2020
                G_diff1=Dim_a.^2*R*Temp*V_int*C0_mn*delta_int/(Dmn_int*Dim_v*(1+2*Dim_a+Dim_a.^2));
                G_diff1=G_diff1*(-exp(Dim_v+Dim_v*Dim_a)+exp(Dim_v+Dim_v*Dim_a)*Dim_v+ ...
                    exp(Dim_v+Dim_v*Dim_a)*Dim_v*Dim_a+1)*exp(-Dim_v-Dim_v*Dim_a);
                G_diff2=-Dim_b*R*Temp*V_int*C0_mn*delta_int/(Dmn_int*Dim_v*(1+Dim_a+2*Dim_b+2*Dim_a*Dim_b+ ...
                    Dim_b^2+Dim_b^2*Dim_a));
                G_diff2=G_diff2*(Dim_a*exp(Dim_v+Dim_v*Dim_b)+Dim_a*Dim_b*exp(Dim_v+Dim_v*Dim_b)+ ...
                Dim_b*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)-Dim_a*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)- ...
                Dim_v*Dim_b*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)+Dim_a*exp(Dim_v+Dim_v*Dim_a)- ...
                Dim_v*Dim_b^2*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)- ...
                Dim_v*Dim_b^2*Dim_a*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)-Dim_a- ...
                Dim_v*Dim_a*Dim_b*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)-Dim_a*Dim_b- ...
                Dim_b*exp(Dim_v+Dim_v*Dim_a))*exp(-2*Dim_v-Dim_v*Dim_b-Dim_v*Dim_a);
                G_diff=(G_diff1+G_diff2)/100; % [J/mol]
                if Gchem<0
                    G_diff_min=minGdiff(Dmn_int,deltaE_Mn,E0,R,Temp,delta_int,C0_mn,-1e-10);
                else
                    G_diff_min=minGdiff(Dmn_int,deltaE_Mn,E0,R,Temp,delta_int,C0_mn,1e-10);
                end
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
                G_friction=V_int./Mobility; % [J/mol]
                if Gchem>=0
                    [sol1]=vpasolve(Gchem-G_friction-G_diff,V_int,[1e-10 4e-5],'random',true);
                else
                    [sol1]=vpasolve(Gchem-G_friction-G_diff,V_int,[-4e-5 -1e-10],'random',true);
                end
                Velocity=double(sol1);
                exitflag=1;
                if length(Velocity)>1
                    SolIndex=find(abs(Velocity)==min(abs(Velocity)));
                    Velocity=Velocity(SolIndex);
                end
                if abs(Velocity)>abs((Gchem-G_diff_min)*Mobility)
                    Velocity=(Gchem-G_diff_min)*Mobility;
                end
                if isempty(sol1)
                    Velocity=0;
                    exitflag=0;
                end
            else
                Velocity=0;
                exitflag=0;
            end
        end
    end
else
        ff=@(x)solute_drag_fun_soft_impingement(x,Temp,C0_mn,wC_A,wC_F,Ux,Xneq,Xpeq,X0,DC,kafang,Rbcc,SN,distance,hardflag);
        opts=optimset('Display','off','Algorithm',{'levenberg-marquardt',0.005},'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-6);% option settings for fsolve
        SolFlag=0;
        v0=[V_previous 1e-6 1e-7 1e-8 1e-9];
        jj=1;
        temp_out=[];
        while SolFlag~=1
            x0=[X0+(Xpeq-X0).*rand(1,1) X0 v0(jj)];
            [out,fval,exitflag]=fsolve(ff,x0,opts); % [Xip Xpm Velocity]
            if jj==length(v0)
                Xip=[];
                Xpm=[];
                Velocity=[];
                DiffLL=distance-Rbcc;
                SolFlag=1;
            end
            if out(1)>X0 && (out(2)>0 && out(2)<0.8*4.6) && exitflag>0 && ~(out(1)>Xpeq && out(2)<out(1)) && ...
                ~(out(3)<0 && (CyclicFlag<1 || (CyclicFlag==1 && tt<CycleStart)))
                Xip=out(1);
                Xpm=out(2);
                Velocity=out(3);
                DiffLL=distance-Rbcc;
                SolFlag=1;
            end
            temp_out=[temp_out;out];
            jj=jj+1;
        end
        softflag=1;
        if isempty(Xip)
            if Rbcc<distance
                Xip=(X0-Rbcc^3/distance^3.*Xneq)./(1-Rbcc^3/distance^3);
                if Xip>Xpeq && TransDir==1 && tt<CycleStart
                    Xip=Xpeq;
                end
            else
                Xip=Xpeq;
            end
            DiffLL=0;
            Xpm=Xip;
            softflag=1;
            Gchem=kafang*(Xpeq-Xip); % [J/mol]
            if Gchem<0
                G_diff_min=minGdiff(Dmn_int,deltaE_Mn,E0,R,Temp,delta_int,C0_mn,-1e-10);
            else
                G_diff_min=minGdiff(Dmn_int,deltaE_Mn,E0,R,Temp,delta_int,C0_mn,1e-10);
            end
            if abs(Gchem)>abs(G_diff_min)
                syms V_int; % [m/s]
                syms Dim_a Dim_b Dim_v G_diff1 G_diff2 G_diff G_friction;
                %% Directly using the integrated results to calculate dissipation energy due to trans-diffusion G_diss
                Dim_a=Dmn_int*(deltaE_Mn-E0)/(R*Temp*V_int*delta_int); % dimensionless parameter a
                Dim_b=Dmn_int*(deltaE_Mn+E0)/(R*Temp*V_int*delta_int); % dimensionless parameter b
                Dim_v=abs(V_int*delta_int/Dmn_int); % dimensionless parameter v
                Dim_v=abs(Dim_v); % add on April 24,2020
                G_diff1=Dim_a.^2*R*Temp*V_int*C0_mn*delta_int/(Dmn_int*Dim_v*(1+2*Dim_a+Dim_a.^2));
                G_diff1=G_diff1*(-exp(Dim_v+Dim_v*Dim_a)+exp(Dim_v+Dim_v*Dim_a)*Dim_v+ ...
                    exp(Dim_v+Dim_v*Dim_a)*Dim_v*Dim_a+1)*exp(-Dim_v-Dim_v*Dim_a);
                G_diff2=-Dim_b*R*Temp*V_int*C0_mn*delta_int/(Dmn_int*Dim_v*(1+Dim_a+2*Dim_b+2*Dim_a*Dim_b+ ...
                    Dim_b^2+Dim_b^2*Dim_a));
                G_diff2=G_diff2*(Dim_a*exp(Dim_v+Dim_v*Dim_b)+Dim_a*Dim_b*exp(Dim_v+Dim_v*Dim_b)+ ...
                Dim_b*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)-Dim_a*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)- ...
                Dim_v*Dim_b*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)+Dim_a*exp(Dim_v+Dim_v*Dim_a)- ...
                Dim_v*Dim_b^2*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)- ...
                Dim_v*Dim_b^2*Dim_a*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)-Dim_a- ...
                Dim_v*Dim_a*Dim_b*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)-Dim_a*Dim_b- ...
                Dim_b*exp(Dim_v+Dim_v*Dim_a))*exp(-2*Dim_v-Dim_v*Dim_b-Dim_v*Dim_a);
                G_diff=(G_diff1+G_diff2)/100; % [J/mol]
                if Gchem<0
                    G_diff_min=minGdiff(Dmn_int,deltaE_Mn,E0,R,Temp,delta_int,C0_mn,-1e-10);
                else
                    G_diff_min=minGdiff(Dmn_int,deltaE_Mn,E0,R,Temp,delta_int,C0_mn,1e-10);
                end
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
                G_friction=V_int./Mobility; % [J/mol]
                if Gchem>=0
                    [sol1]=vpasolve(Gchem-G_friction-G_diff,V_int,[1e-10 4e-5],'random',true);
                else
                    [sol1]=vpasolve(Gchem-G_friction-G_diff,V_int,[-4e-5 -1e-10],'random',true);
                end
                Velocity=double(sol1);
                if length(Velocity)>1
                    SolIndex=find(abs(Velocity)==min(abs(Velocity)));
                    Velocity=Velocity(SolIndex);
                    exitflag=1;
                end
                if abs(Velocity)>abs((Gchem-G_diff_min)*Mobility)
                    Velocity=(Gchem-G_diff_min)*Mobility;
                end
                if isempty(sol1)
                    Velocity=0;
                    exitflag=0;
                end
            else
                Velocity=0;
                exitflag=0;
            end
            
        end
end

%% Directly using the integrated results to calculate dissipation energy due to trans-diffusion G_diss
clear V_int Dim_a Dim_b Dim_v G_diff1 G_diff2 G_diff G_friction;
Dim_a=Dmn_int*(deltaE_Mn-E0)/(R*Temp*Velocity*delta_int); % dimensionless parameter a
Dim_b=Dmn_int*(deltaE_Mn+E0)/(R*Temp*Velocity*delta_int); % dimensionless parameter b
Dim_v=abs(Velocity*delta_int/Dmn_int); % dimensionless parameter v
Gdiff1=Dim_a.^2*R*Temp*Velocity*C0_mn*delta_int/(Dmn_int*Dim_v*(1+2*Dim_a+Dim_a.^2));
Gdiff1=Gdiff1*(-exp(Dim_v+Dim_v*Dim_a)+exp(Dim_v+Dim_v*Dim_a)*Dim_v+ ...
    exp(Dim_v+Dim_v*Dim_a)*Dim_v*Dim_a+1)*exp(-Dim_v-Dim_v*Dim_a);
Gdiff2=-Dim_b*R*Temp*Velocity*C0_mn*delta_int/(Dmn_int*Dim_v*(1+Dim_a+2*Dim_b+2*Dim_a*Dim_b+ ...
    Dim_b^2+Dim_b^2*Dim_a));
Gdiff2=Gdiff2*(Dim_a*exp(Dim_v+Dim_v*Dim_b)+Dim_a*Dim_b*exp(Dim_v+Dim_v*Dim_b)+ ...
    Dim_b*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)-Dim_a*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)- ...
    Dim_v*Dim_b*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)+Dim_a*exp(Dim_v+Dim_v*Dim_a)- ...
    Dim_v*Dim_b^2*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)- ...
    Dim_v*Dim_b^2*Dim_a*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)-Dim_a- ...
    Dim_v*Dim_a*Dim_b*exp(2*Dim_v+Dim_v*Dim_b+Dim_v*Dim_a)-Dim_a*Dim_b- ...
    Dim_b*exp(Dim_v+Dim_v*Dim_a))*exp(-2*Dim_v-Dim_v*Dim_b-Dim_v*Dim_a);
Gdiff=(Gdiff1+Gdiff2)/100; % [J/mol]

Gdiff0=(2*C0_mn*deltaE_Mn+R*Temp*C0_mn*(exp(-2*deltaE_Mn/(R*Temp))-1))/100; % dissipation when v_int = 0[J/mol]
if Xip>Xpeq
    Gdiff0=minGdiff(Dmn_int,deltaE_Mn,E0,R,Temp,delta_int,C0_mn,-1e-10);
end
if isnan(Gdiff) % Gdiff is NaN
    Gdiff=(2*C0_mn*deltaE_Mn+R*Temp*C0_mn*(exp(-2*deltaE_Mn/(R*Temp))-1))/100; % [J/mol]
    if Xip>Xpeq
        Gdiff=minGdiff(Dmn_int,deltaE_Mn,E0,R,Temp,delta_int,C0_mn,-1e-10);
    end
end
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
Gfriction=Velocity./Mobility; % [J/mol]
Gdiss=Gfriction+Gdiff; % [J/mol]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gchem=kafang*(Xpeq-Xip); % [J/mol]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Drag=[Velocity Gchem Gfriction Gdiff Gdiss exitflag]';
Diffusion=[Xpm Xip DiffLL Rbcc Xneq Xpeq softflag]';

