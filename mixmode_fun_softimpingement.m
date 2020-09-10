function SS=mixmode_fun_softimpingement(xx,C0,Cneq,Cpeq,Rbcc,Z,distance,SN)
% This function is to calculate the interfatial C
% concentration, diffusion length x(2)and middle C content under overlapping condition
Vbcc=4/3*pi*Rbcc^3;
SS(1)=30*Vbcc*(C0-Cneq)/(4*pi)*1/SN-(9*xx(2)^3*xx(3)-10*xx(2)^3*C0+xx(2)^3*xx(1)+20*Rbcc^2*xx(3)*xx(2)- ...
    30*Rbcc^2*C0*xx(2)+10*Rbcc^2*xx(1)*xx(2)+25*xx(2)^2*Rbcc*xx(3)-30*xx(2)^2*Rbcc*C0+5*xx(2)^2*Rbcc*xx(1))*1/SN; % mass conservation
SS(2)=(Cpeq-xx(1))*(xx(1)-Cneq)*xx(2)-Z*(xx(1)-xx(3)); % no accumulation on the interface
SS(3)=xx(2)-(distance-Rbcc);
end