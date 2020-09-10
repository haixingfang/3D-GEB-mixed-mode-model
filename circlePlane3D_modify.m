function H = circlePlane3D( center, normal, radious, theintv, normalon, style,Lb )
%CIRCLEPLANE3D Summary of this function goes here
%--------------------------------------------------------------------------
%Generate a circle plane in 3D with the given center and radious
%The plane is defined by the normal vector
%theintv is the interval theta which allow you to control your polygon
%shape
% Example:,
%
%   circlePlane3D([0 0 0], [1 -1 2], 5, 0.2, 1, [0 0 1], '-'); 
%   circlePlane3D([3 3 -3],[0 1 1], 3, 0.1, 1, 'y', '-');
%   
%   Cheng-Yuan Wu <ieda_wind@hotmail.com>
%   Version 1.00
%   Aug, 2012
%   Modified on June 27, 2015
%--------------------------------------------------------------------------
%generate circle polygon
t = 0:theintv:2*pi;
flag=find(abs(center)==Lb/2);

if flag==1
   fx=center(1)+zeros(1,length(t));
   fy=radious*cos(t)+center(2);
   fz=radious*sin(t)+center(3);
else if flag==2
   fx=center(1)+radious*cos(t);
   fy=center(2)+zeros(1,length(t));
   fz=radious*sin(t)+center(3);
    else if flag==3
   fx=center(1)+radious*cos(t);
   fy=radious*sin(t)+center(2);
   fz=center(3)+zeros(1,length(t));
        end
    end
end

% color=zeros(size(fx,2),3);
% for i=1:size(fx,2)
%     color(i,3)=1;
% end
H = fill3(fx, fy, fz,'b');
set(H,'edgecolor','none','FaceLighting','phong');

%set(H,'FaceColor','flat','FaceVertexCData',color);


% if normalon == 1
%     hold on;
%     H = plot3([center(1) center(1)+normal(1)],[center(2) center(2)+normal(2)],[center(3) center(3)+normal(3)],style);
% end
