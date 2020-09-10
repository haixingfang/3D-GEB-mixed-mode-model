function FF=mixmode_fun(x,C0,Cneq,Cpeq,Rbcc,Z,SN)
% This function is to calculate the interfatial C
% concentration and diffusion length x(2)under nonoverlapping conditionVbcc=/
Vbcc=4/3*pi*Rbcc^3;
FF(1)=30*Vbcc*(C0-Cneq)/(4*pi*(x(1)-C0))*1/SN-(x(2)^3+5*Rbcc*x(2)^2+10*Rbcc^2*x(2))*1/SN; % mass conservation
FF(2)=(Cpeq-x(1))*(x(1)-Cneq)*x(2)-Z*(x(1)-C0); % no accumulation on the interface
end

