function [position_MXY position_MXZ position_MYZ position_MO]=MirrorPoints(position,Lb)
% This function is to find four mirror points of interest for original
% point position: N0*3 array
% Lb is the box size
% Returns the coordinates of the mirror points:  
% model

% Find mirror points
for i=1:length(position)
    if position(i,1)>=0&&position(i,1)<=Lb/2
       if position(i,2)>=0&&position(i,2)<=Lb/2  
           if position(i,3)>=0&&position(i,3)<=Lb/2  % The 1st Octant
             position_MXY(i,1)=position(i,1);
             position_MXY(i,2)=position(i,2);
             position_MXY(i,3)=position(i,3)-Lb;
             
             position_MXZ(i,1)=position(i,1);
             position_MXZ(i,2)=position(i,2)-Lb;
             position_MXZ(i,3)=position(i,3);
             
             position_MYZ(i,1)=position(i,1)-Lb;
             position_MYZ(i,2)=position(i,2);
             position_MYZ(i,3)=position(i,3);
             
             position_MO(i,1)=position(i,1)-Lb;
             position_MO(i,2)=position(i,2)-Lb;
             position_MO(i,3)=position(i,3)-Lb;
           else                                    % The 5th Octant
             position_MXY(i,1)=position(i,1);
             position_MXY(i,2)=position(i,2);
             position_MXY(i,3)=position(i,3)+Lb;
             
             position_MXZ(i,1)=position(i,1);
             position_MXZ(i,2)=position(i,2)-Lb;
             position_MXZ(i,3)=position(i,3);
             
             position_MYZ(i,1)=position(i,1)-Lb;
             position_MYZ(i,2)=position(i,2);
             position_MYZ(i,3)=position(i,3);
             
             position_MO(i,1)=position(i,1)-Lb;
             position_MO(i,2)=position(i,2)-Lb;
             position_MO(i,3)=position(i,3)+Lb;
           end
       else
           if position(i,3)>=0&&position(i,3)<=Lb/2  % The 4th Octant
             position_MXY(i,1)=position(i,1);
             position_MXY(i,2)=position(i,2);
             position_MXY(i,3)=position(i,3)-Lb;
             
             position_MXZ(i,1)=position(i,1);
             position_MXZ(i,2)=position(i,2)+Lb;
             position_MXZ(i,3)=position(i,3);
             
             position_MYZ(i,1)=position(i,1)-Lb;
             position_MYZ(i,2)=position(i,2);
             position_MYZ(i,3)=position(i,3);
             
             position_MO(i,1)=position(i,1)-Lb;
             position_MO(i,2)=position(i,2)+Lb;
             position_MO(i,3)=position(i,3)-Lb;
           else                                    % The 8th Octant
             position_MXY(i,1)=position(i,1);
             position_MXY(i,2)=position(i,2);
             position_MXY(i,3)=position(i,3)+Lb;
             
             position_MXZ(i,1)=position(i,1);
             position_MXZ(i,2)=position(i,2)+Lb;
             position_MXZ(i,3)=position(i,3);
             
             position_MYZ(i,1)=position(i,1)-Lb;
             position_MYZ(i,2)=position(i,2);
             position_MYZ(i,3)=position(i,3);
             
             position_MO(i,1)=position(i,1)-Lb;
             position_MO(i,2)=position(i,2)+Lb;
             position_MO(i,3)=position(i,3)+Lb;
           end
       end
    else
        if position(i,2)>=0&&position(i,2)<=Lb/2  
           if position(i,3)>=0&&position(i,3)<=Lb/2  % The 2nd Octant
             position_MXY(i,1)=position(i,1);
             position_MXY(i,2)=position(i,2);
             position_MXY(i,3)=position(i,3)-Lb;
             
             position_MXZ(i,1)=position(i,1);
             position_MXZ(i,2)=position(i,2)-Lb;
             position_MXZ(i,3)=position(i,3);
             
             position_MYZ(i,1)=position(i,1)+Lb;
             position_MYZ(i,2)=position(i,2);
             position_MYZ(i,3)=position(i,3);
             
             position_MO(i,1)=position(i,1)+Lb;
             position_MO(i,2)=position(i,2)-Lb;
             position_MO(i,3)=position(i,3)-Lb;
           else                                    % The 6th Octant
             position_MXY(i,1)=position(i,1);
             position_MXY(i,2)=position(i,2);
             position_MXY(i,3)=position(i,3)+Lb;
             
             position_MXZ(i,1)=position(i,1);
             position_MXZ(i,2)=position(i,2)-Lb;
             position_MXZ(i,3)=position(i,3);
             
             position_MYZ(i,1)=position(i,1)+Lb;
             position_MYZ(i,2)=position(i,2);
             position_MYZ(i,3)=position(i,3);
             
             position_MO(i,1)=position(i,1)+Lb;
             position_MO(i,2)=position(i,2)-Lb;
             position_MO(i,3)=position(i,3)+Lb;
           end
       else
           if position(i,3)>=0&&position(i,3)<=Lb/2  % The 3rd Octant
             position_MXY(i,1)=position(i,1);
             position_MXY(i,2)=position(i,2);
             position_MXY(i,3)=position(i,3)-Lb;
             
             position_MXZ(i,1)=position(i,1);
             position_MXZ(i,2)=position(i,2)+Lb;
             position_MXZ(i,3)=position(i,3);
             
             position_MYZ(i,1)=position(i,1)+Lb;
             position_MYZ(i,2)=position(i,2);
             position_MYZ(i,3)=position(i,3);
             
             position_MO(i,1)=position(i,1)+Lb;
             position_MO(i,2)=position(i,2)+Lb;
             position_MO(i,3)=position(i,3)-Lb;
           else                                    % The 7th Octant
             position_MXY(i,1)=position(i,1);
             position_MXY(i,2)=position(i,2);
             position_MXY(i,3)=position(i,3)+Lb;
             
             position_MXZ(i,1)=position(i,1);
             position_MXZ(i,2)=position(i,2)+Lb;
             position_MXZ(i,3)=position(i,3);
             
             position_MYZ(i,1)=position(i,1)+Lb;
             position_MYZ(i,2)=position(i,2);
             position_MYZ(i,3)=position(i,3);
             
             position_MO(i,1)=position(i,1)+Lb;
             position_MO(i,2)=position(i,2)+Lb;
             position_MO(i,3)=position(i,3)+Lb;
           end
             end
    end
end
