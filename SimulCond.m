function [Comp T_a3 T_a1 Tc a_eq b_eq kafang_p T_wC Temp_eq Feq_max RT hold_t cycN cycT CyclicFlag]=SimulCond()

CyclicFlag=1; % 1-activate the cycling; 0-isothermal; -1-continuous cooling

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Fe-0.1C-0.17Mn
% Comp=[0.1 0.17 0];% C Mn Si in wt.%
% T_a3=1132;   % A3-temperature under PE [K]
% T_a1=[995 994]; % A1-temperature othorequilibrium, A1+, A1-, [K]
% Tc=1043;     % Curie temperature [K]
% a_eq=[1.50936E-05 -3.69719E-02 2.26099E+01]; % PE
% b_eq=[-8.94382E-05 1.06755E-01]; % PE
% kafang_p=[0 -0.1572 269.44]; % kafang(T)=a*T^2+b*T+c [J/(mol*mol%)]
% T_wC=[6.05521E+01 -2.59977E+02 1.15750E+03]; % TA3=a*wC^2+b*wC+c under PE [K] (wC in wt.%)
% Temp_eq=[1136 995]; % A3-ortho A1+ [K]
% Feq_max=0.985; % maximum ferrite fraction
% %%%%% cycling settings
% RT=[20*60 10*60 10*60 20*60]./60;% First cooling rate; cycling heating rate; cycling cooling rate; final cooling rate [oC/s]
% hold_t=[3 0 0].*60; % holding time for first duration, duration after heating, duration after cooling [s]
% cycN=3; % number of cycles
% cycT=[800 854]+273; % cycling temperature [K]

% Fe-0.1C-0.5Mn
Comp=[0.1 0.49 0];% C Mn Si in wt.% the same as in M.Militzer Acta Mater. 2006 Phase field modelling of austenite to ferrite transformation
T_a3=1118;   % A3-temperature othorequilibrium=1125K, Para-equilibrium=1117K, linerization of PE from M.G.Meccozi=1106K [K]
T_a1=[986 983]; % A1-temperature othorequilibrium, A1+, A1-, [K]
Tc=1043-5*0.49;     % Curie temperature [K]
a_eq=[1.55572e-5 -3.78004e-2 2.29105e1]; % PE
b_eq=[-9.09127e-5 1.06277e-1]; % PE
kafang_p=[0 -0.1661 277.66]; % kafang(T)=a*T^2+b*T+c [J/(mol*mol%)]
T_wC=[5.76665e1 -2.48065e2 1.14063e3]; % TA3=a*wC^2+b*wC+c under PE [K] (wC in wt.%)
Temp_eq=[1124 986]; % A3-ortho A1+ [K]
Feq_max=0.985; % maximum ferrite fraction
%%%%% cycling settings
RT=[20*60 10*60 10*60 20*60]./60;% First cooling rate; cycling heating rate; cycling cooling rate; final cooling rate [oC/s]
hold_t=[3 0 0].*60; % holding time for first duration, duration after heating, duration after cooling [s]
cycN=1; % number of cycles
cycT=[785 842]+273; % cycling temperature [K]
% 
% % Fe-0.1C-1.0Mn
% Comp=[0.1 1.0 0];% C Mn Si in wt.%
% T_a3=1096;   % A3-temperature under PE [K]
% T_a1=[973 966]; % A1-temperature othorequilibrium, A1+, A1-, [K]
% Tc=1043-5;     % Curie temperature [K]
% a_eq=[1.62248E-05 -3.89762E-02 2.33212E+01]; % PE
% b_eq=[-9.00827E-05 1.02737E-01]; % PE
% kafang_p=[0 -0.1892 297.6]; % kafang(T)=a*T^2+b*T+c [J/(mol*mol%)]
% T_wC=[5.31265E+01 -2.30403E+02 1.11586E+03]; % TA3=a*wC^2+b*wC+c under PE [K] (wC in wt.%)
% Temp_eq=[1107 973]; % A3-ortho A1+ [K]
% Feq_max=0.985; % maximum ferrite fraction
% %%%%% cycling settings
% RT=[20*60 10*60 10*60 20*60]./60;% First cooling rate; cycling heating rate; cycling cooling rate; final cooling rate [oC/s]
% hold_t=[3 0 0].*60; % holding time for first duration, duration after heating, duration after cooling [s]
% cycN=3; % number of cycles
% cycT=[767 822]+273; % cycling temperature [K]
% 
% % Fe-0.1C-1.5Mn
% Comp=[0.1 1.5 0];% C Mn Si in wt.%
% T_a3=1073;   % A3-temperature under PE [K]
% T_a1=[959 947]; % A1-temperature othorequilibrium, A1+, A1-, [K]
% Tc=1043-10;     % Curie temperature [K]
% a_eq=[1.67951E-05 -3.99534E-02 2.36303E+01]; % PE
% b_eq=[-8.88324E-05 9.89725E-02]; % PE
% kafang_p=[0 -0.2463 350.74]; % kafang(T)=a*T^2+b*T+c [J/(mol*mol%)]
% T_wC=[4.81792E+01 -2.12909E+02 1.09194E+03]; % TA3=a*wC^2+b*wC+c under PE [K] (wC in wt.%)
% Temp_eq=[1092 959]; % A3-ortho A1+ [K]
% Feq_max=0.985; % maximum ferrite fraction
% %%%%% cycling settings
% RT=[20*60 10*60 10*60 20*60]./60;% First cooling rate; cycling heating rate; cycling cooling rate; final cooling rate [oC/s]
% hold_t=[3 0 0].*60; % holding time for first duration, duration after heating, duration after cooling [s]
% cycN=3; % number of cycles
% cycT=[747 797]+273; % cycling temperature [K]
% 
% % Fe-0.1C-2.0Mn
% Comp=[0.1 2.0 0];% C Mn Si in wt.%
% T_a3=1051;   % A3-temperature under PE [K]
% T_a1=[947 926]; % A1-temperature othorequilibrium, A1+, A1-, [K]
% Tc=1043-15;     % Curie temperature [K]
% a_eq=[1.73870E-05 -4.09646E-02 2.39511E+01]; % PE
% b_eq=[-8.82830E-05 9.60369E-02]; % PE
% kafang_p=[0 -0.2643 362.77]; % kafang(T)=a*T^2+b*T+c [J/(mol*mol%)]
% T_wC=[4.63613E+01 -2.01776E+02 1.07100E+03]; % TA3=a*wC^2+b*wC+c under PE [K] (wC in wt.%)
% Temp_eq=[1077 926]; % A3-ortho A1+ [K]
% Feq_max=0.985; % maximum ferrite fraction
% %%%%% cycling settings
% RT=[20*60 10*60 10*60 20*60]./60;% First cooling rate; cycling heating rate; cycling cooling rate; final cooling rate [oC/s]
% hold_t=[30 0 0].*60; % holding time for first duration, duration after heating, duration after cooling [s]
% cycN=3; % number of cycles
% %cycT=[728 777]+273; % cycling temperature [K]
% cycT=[680 730]+273; % cycling temperature [K]
% %cycT=[666 710]+273;
% 
% % Fe-0.1C-2.5Mn
% Comp=[0.1 2.5 0];% C Mn Si in wt.%
% T_a3=1033;   % A3-temperature under PE [K]
% T_a1=[934 903]; % A1-temperature othorequilibrium, A1+, A1-, [K]
% Tc=1043-20;     % Curie temperature [K]
% a_eq=[1.79303E-05 -4.18770E-02 2.42206E+01]; % PE
% b_eq=[-8.74750E-05 9.30111E-02]; % PE
% kafang_p=[0 -0.2624 356.69]; % kafang(T)=a*T^2+b*T+c [J/(mol*mol%)]
% T_wC=[4.40999E+01 -1.90293E+02 1.05051E+03]; % TA3=a*wC^2+b*wC+c under PE [K] (wC in wt.%)
% Temp_eq=[1063 934]; % A3-ortho A1+ [K]
% Feq_max=0.985; % maximum ferrite fraction
% %%%%% cycling settings
% RT=[20*60 10*60 10*60 20*60]./60;% First cooling rate; cycling heating rate; cycling cooling rate; final cooling rate [oC/s]
% hold_t=[3 0 0].*60; % holding time for first duration, duration after heating, duration after cooling [s]
% cycN=3; % number of cycles
% cycT=[710 757]+273; % cycling temperature [K]
% 
% % Fe-0.25C-0.17Mn
% Comp=[0.25 0.17 0];% C Mn Si in wt.%
% T_a3=1089;   % A3-temperature under PE [K]
% T_a1=[997 995]; % A1-temperature othorequilibrium, A1+, A1-, [K]
% Tc=1043;     % Curie temperature [K]
% a_eq=[1.55725E-05 -3.79852E-02 2.31433E+01]; % PE
% b_eq=[-9.42502E-05 1.11739E-01]; % PE
% kafang_p=[0 -0.2435 359.46]; % kafang(T)=a*T^2+b*T+c [J/(mol*mol%)]
% T_wC=[6.05521E+01 -2.59970E+02 1.15749E+03]; % TA3=a*wC^2+b*wC+c under PE [K] (wC in wt.%)
% Temp_eq=[1091 997]; % A3-ortho A1+ [K]
% Feq_max=0.985; % maximum ferrite fraction
% %%%%% cycling settings
% RT=[20*60 10*60 10*60 20*60]./60;% First cooling rate; cycling heating rate; cycling cooling rate; final cooling rate [oC/s]
% hold_t=[3 0 0].*60; % holding time for first duration, duration after heating, duration after cooling [s]
% cycN=3; % number of cycles
% cycT=[725 802]+273; % cycling temperature [K]
% 
% % Fe-0.247C-2.06Mn
% Comp=[0.247 2.06 0];% C Mn Si in wt.% of the sample collected from Ancelor Mittar
% T_a3=1019;   % A3-temperature othorequilibrium=1045K, Para-equilibrium=1019K,experimental show A3=1034K
% T_a1=[963 938]; % A1+=963 K, A1-=938 K
% Tc=1034; % ND experiment show Curie temperature [K]
% a_eq=[1.74985e-5 -4.11663e-2 2.40289e1]; % PE 
% b_eq=[-8.85649e-5 9.60186e-2]; % PE   
% % a_eq=[0.0000215981 -0.05069811 29.64121]; % LE 
% % b_eq=[-0.0001131724 0.1247617]; % LE
% kafang_p=[0 -0.2249 323.53]; % kafang(T)=a*T^2+b*T+c [J/(mol*mol%)]
% T_wC=[4.69338e1 -2.01709e2 1.06888e3]; % TA3=a*wC^2+b*wC+c under PE [K] (wC in wt.%)
% Temp_eq=[1045 963];% A1-=938K; A3_para=1019K
% Feq_max=0.962; % maximum ferrite fraction
% %%%%% cycling settings
% RT=[20*60 10 10 20*60]./60;% First cooling rate; cycling heating rate; cycling cooling rate; final cooling rate [oC/s]
% hold_t=[2 1 1].*60; % holding time for first duration, duration after heating, duration after cooling [s]
% cycN=1; % number of cycles
% cycT=[710 730]+273; % cycling temperature [K]
% 
% % Fe-0.023C-0.17Mn
% Comp=[0.023 0.17 0];% C Mn Si in wt.%
% T_a3=1163;   % A3-temperature under PE [K]
% T_a1=[994 994]; % A1-temperature othorequilibrium, A1+, A1-, [K]
% Tc=1043;     % Curie temperature [K]
% a_eq=[1.55700E-05 -3.79799E-02 2.31405E+01]; % PE
% b_eq=[-9.42383E-05 1.11729E-01]; % PE
% kafang_p=[0 -0.1201 230.15]; % kafang(T)=a*T^2+b*T+c [J/(mol*mol%)]
% T_wC=[6.05520E+01 -2.59981E+02 1.15751E+03]; % TA3=a*wC^2+b*wC+c under PE [K] (wC in wt.%)
% Temp_eq=[1166 994]; % A3-ortho A1+ [K]
% Feq_max=0.985; % maximum ferrite fraction
% %%%%% cycling settings
% RT=[20*60 10*60 10*60 20*60]./60;% First cooling rate; cycling heating rate; cycling cooling rate; final cooling rate [oC/s]
% hold_t=[3 0 0].*60; % holding time for first duration, duration after heating, duration after cooling [s]
% cycN=3; % number of cycles
% cycT=[860 885]+273; % cycling temperature [K]

% % Fe-0.05C-2Mn
% Comp=[0.05 2.0 0];% C Mn Si in wt.%
% T_a3=1064;   % A3-temperature othorequilibrium=1090K, Para-equilibrium=1064K
% T_a1=[935 920]; % A1+=963 K, A1-=938 K
% Tc=1034; % ND experiment show Curie temperature [K]
% a_eq=[1.73880E-05 -4.09653E-02 2.39509E+01]; % PE 
% b_eq=[-8.82479E-05 9.60050E-02]; % PE   
% kafang_p=[0 -0.2318 331.11]; % kafang(T)=a*T^2+b*T+c [J/(mol*mol%)]
% T_wC=[4.63905E+01 -2.01848E+02 1.07104E+03]; % TA3=a*wC^2+b*wC+c under PE [K] (wC in wt.%)
% Temp_eq=[1090 920];
% Feq_max=0.992; % maximum ferrite fraction
% %%%%% cycling settings
% RT=[20*60 10 10 20*60]./60;% First cooling rate; cycling heating rate; cycling cooling rate; final cooling rate [oC/s]
% hold_t=[2 1 1].*60; % holding time for first duration, duration after heating, duration after cooling [s]
% cycN=1; % number of cycles
% cycT=[725 789]+273; % cycling temperature [K]

if CyclicFlag==0 % isothermal
    RT(2)=0;
    RT(3)=0;
    cycN=0;
    hold_t(1)=10*60; % [s]
end
if CyclicFlag==-1 % continuous cooling
    RT(2)=0;
    RT(3)=0;
    cycN=0;
    hold_t(1)=0;
    RT(1)=10; % [K/s]
    RT(4)=10; % [K/s]
end

end
