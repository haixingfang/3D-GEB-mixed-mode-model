function [Plane0]=CutOffPlot(sphere0,Lb)
% sphere0=[x,y,z,r] the coordinates of the sphere center and the radius
% Lb is the box size

%   PLANE  : [x0 y0 z0  dx1 dy1 dz1  dx2 dy2 dz2]

Plane(1,1,:)=[0 0 Lb/2];
Plane(1,2,:)=[Lb/2 Lb/2 Lb/2]-[0 0 Lb/2]; % Direction vector1
Plane(1,3,:)=[-Lb/2 Lb/2 Lb/2]-[0 0 Lb/2]; % Direction vector2

Plane(2,1,:)=[0 0 -Lb/2];
Plane(2,2,:)=[Lb/2 Lb/2 -Lb/2]-[0 0 -Lb/2]; % Direction vector1
Plane(2,3,:)=[-Lb/2 Lb/2 -Lb/2]-[0 0 -Lb/2]; % Direction vector2

Plane(3,1,:)=[-Lb/2 0 0];
Plane(3,2,:)=[-Lb/2 Lb/2 Lb/2]-[-Lb/2 0 0]; % Direction vector1
Plane(3,3,:)=[-Lb/2 -Lb/2 Lb/2]-[-Lb/2 0 0]; % Direction vector2

Plane(4,1,:)=[Lb/2 0 0];
Plane(4,2,:)=[Lb/2 Lb/2 Lb/2]-[Lb/2 0 0]; % Direction vector1
Plane(4,3,:)=[Lb/2 -Lb/2 Lb/2]-[Lb/2 0 0]; % Direction vector2

Plane(5,1,:)=[0 Lb/2 0];
Plane(5,2,:)=[Lb/2 Lb/2 Lb/2]-[0 Lb/2 0]; % Direction vector1
Plane(5,3,:)=[-Lb/2 Lb/2 Lb/2]-[0 Lb/2 0]; % Direction vector2

Plane(6,1,:)=[0 -Lb/2 0];
Plane(6,2,:)=[Lb/2 -Lb/2 Lb/2]-[0 -Lb/2 0]; % Direction vector1
Plane(6,3,:)=[-Lb/2 -Lb/2 Lb/2]-[0 -Lb/2 0]; % Direction vector2


theta=0:pi/50:2*pi;
for i=1:6
    Plane00=permute(Plane(i,:,:),[3 2 1]);
    Plane0=reshape(Plane00,[1 9]);
    [center0 Rc0 nor0]=intersectPlaneSphere(Plane0, sphere0);
    %center(i,:)=center0;
    %plotCircle3D(center(i,:),nor0,Rc0);
    if Rc0>0
       H=circlePlane3D_modify(center0,nor0,Rc0,pi/50,1,'',Lb);
    end
    hold on;
end
xlabel('x (\mum)','FontSize',16);
ylabel('y (\mum)','FontSize',16);
zlabel('z (\mum)','FontSize',16);
axis([-Lb/2 Lb/2 -Lb/2 Lb/2 -Lb/2 Lb/2]);
box on;

