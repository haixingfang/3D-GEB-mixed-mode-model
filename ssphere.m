function ssphere(x,y,z,r,Lb)%function ssphere
[x0,y0,z0]=sphere(20);
xx=x+x0*r;
yy=y+y0*r;
zz=z+z0*r;
figure(1);

h=surf(xx,yy,zz);
%colormap(hot);
colormap([1 0 0]);   % [R G B]
axis equal;
set(gca,'DataAspectRatio',[1 1 1])    
set(gca,'XDir','rev','YDir','rev'); 
shading Interp;
%set(gca,'BoxStyle','full','Box','on');
set(gca,'Box','on');
lightangle(45,45);
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.3;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';
grid off;
%axis([-Lb/2 Lb/2 -Lb/2 Lb/2 -Lb/2 Lb/2]);
xlim([-Lb/2 Lb/2]);
ylim([-Lb/2 Lb/2]);
zlim([-Lb/2 Lb/2]);
xlabel('x (\mum)','FontSize',16);
ylabel('y (\mum)','FontSize',16);
zlabel('z (\mum)','FontSize',16);
% x=3;
% y=10;
% z=-5;
% r=30;
% xx=x+x0*r;
% yy=y+y0*r;
% zz=z+z0*r;
% 
% h=surf(xx,yy,zz);
% patch(surf2patch(h));

% V = surface2volume2(xx,yy,zz,0.05); % Convert face to volume data set
% isoval=0.5;
% p = patch(isosurface(V,isoval), ...
%     'FaceColor','red','EdgeColor','red', ...
%     'AmbientStrength',0.3, ...
%     'DiffuseStrength',0.8, ...
%     'SpecularStrength',0.7, ...
%     'SpecularExponent',25, ...
%     'BackFaceLighting','unlit');
% %isonormals(V,p);
% 
% % Create the isocaps
% patch(isocaps(V,isoval),...
%  'FaceColor','interp',...
%  'EdgeColor','none')
% colormap hsv
% 
% % Define the view
% daspect([1,1,1])
% axis tight
% view(3)



%plot3(xx,yy,zeros(size(zz)));
%plot3(zeros(size(xx)),yy,zz);
%fill3(xx,zeros(size(yy)),zz);