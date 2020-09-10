%This function return the diffusion profile of C based on quadratic
%diffusion profile

function [Diffusion]=Mixmode_diffusion_profile_vpasolve(Xneq,Xpeq,X0,DC,Mob,kafang,Rbcc,SN,distance,hardflag)
% Diffusion=[Xpm Xip DiffLL Rbcc Xneq Xpeq];
% mol%, mol%, um, um, mol%, mol%

% clear all;
% % Fe-1.0 at.%C as an example
% T=1050; % [K]
% X0=1; % [at.%]
% distance=10; % [um]
% V0=4/3*pi*distance^3; % [um^3]
% Xneq=0.0658; % [at.%]
% Xpeq=2.0619; % [at.%]
% DC=1.14e-12; % [m2/s]
% kafang=110; % [J/at.%]
% Mob=5.4e-8; % [mol.m/(J.s)]
% feq=0.1383;
% feq=0.3;
% Vbcc=V0*feq; % [um^3]
% % Abcc=4*pi*(3*Vbcc/(4*pi))^(2/3); % [um^2]
% Rbcc=(3*Vbcc/(4*pi))^(1/3); % [um]
% SN=1;

% %%%%%%%%%%%%%%%%%%%%%%%%%% test the solution by inputting known, updated
% % on May 24, 2018
% Xneq=xC_F_eq(i);
% Xpeq=xC_A_eq(i);
% X0=N_p(l,34);
% DC=D_C(i)*1e-12;
% Mob=Mobility(i);
% kafang=Kafang(i);
% Rbcc=N_p(l,7);
% SN=length(N_PR{l}(:,1));
% distance=N_p(l,32);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vbcc=4/3*pi*Rbcc^3; % [um^3]
Z=2*DC/(Mob*kafang)*1e6; % [um*at.%]

syms X1; % unknown interfacial C concentration in austenite
syms X2; % unknown diffusion length

[sol1 sol2]=vpasolve(30*Vbcc*(X0-Xneq)/(4*pi*(X1-X0))*1/SN-(X2^3+5*Rbcc*X2^2+10*Rbcc^2*X2)*1/SN, ...
    (Xpeq-X1)*(X1-Xneq)*X2-Z*(X1-X0), ...
    [X1 X2],[0.05,(X0-0.9*0.05)/(1-0.9);0,distance],'random',true);
Xip=double(sol1);
DiffLL=double(sol2);
Xpm=X0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FF=@(x)mixmode_fun(x,X0,Xneq,Xpeq,Rbcc,Z,SN);
% opts=optimset('Display','iter','Algorithm',{'levenberg-marquardt',0.005},'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-6);% option settings for fsolve
% x0(1)=Xpeq; % initial Xip [at.%]
% % x0(2)=distance-Rbcc; % initial diffusion length [um]
% x0(2)=0.05;
% [sol1,fval,exitflag]=fsolve(FF,x0,opts);
% Xip=sol1(1);
% DiffLL=sol1(2);
% Xpm=X0;

if DiffLL<=(distance-Rbcc) && hardflag==0
    softflag=0;
else
    syms X3;
%     Vbcc=4/3*pi*Rbcc^3;
    [sol1 sol2 sol3]=vpasolve(30*Vbcc*(X0-Xneq)/(4*pi)*1/SN-(9*X2^3*X3-10*X2^3*X0+X2^3*X1+20*Rbcc^2*X3*X2- ...
        30*Rbcc^2*X0*X2+10*Rbcc^2*X1*X2+25*X2^2*Rbcc*X3-30*X2^2*Rbcc*X0+5*X2^2*Rbcc*X1)*1/SN, ...
        (Xpeq-X1)*(X1-Xneq)*X2-Z*(X1-X3), ...
        X2-(distance-Rbcc),[X1,X2,X3],[0.05,(X0-0.9*0.05)/(1-0.9);0,distance;0.05,(X0-0.9*0.05)/(1-0.9)],'random',true);
    Xip=double(sol1);
    DiffLL=double(sol2);
    Xpm=double(sol3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     SS=@(xx)mixmode_fun_softimpingement(xx,X0,Xneq,Xpeq,Rbcc,Z,distance,SN)
%     opts=optimset('Display','iter','Algorithm',{'levenberg-marquardt',0.005},'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-6);% option settings for fsolve
%     xx0(1)=(X0-(Rbcc/distance)^3*Xneq)/(1-(Rbcc/distance)^3); % initial Xip [at.%]
%     xx0(2)=distance-Rbcc; % initial diffusion length [um]
%     xx0(3)=(X0-(Rbcc/distance)^3*Xneq)/(1-(Rbcc/distance)^3); % initial central C content[at.%]
%     [sol2,fval,exitflag]=fsolve(SS,xx0,opts);
%     Xip=sol2(1);
%     DiffLL=sol2(2);
%     Xpm=sol2(3);
    softflag=1;
end
Diffusion=[Xpm Xip DiffLL Rbcc Xneq Xpeq softflag]';
Xip=(X0-Rbcc^3/distance^3.*Xneq)./(1-Rbcc^3/distance^3);

end



