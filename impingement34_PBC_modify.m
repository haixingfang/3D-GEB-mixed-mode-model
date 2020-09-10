%This function return the volume of the sphere i in which the
%impingement of three grains considered;
%The overlay of three spheres has been subtracted when calculating the
%volume i;
%This function is specially designed for ferrite_PBC model
function [V Tri impinge4_flag]=impingement34_PBC_modify(x,numb,flag4)
%x is the N0*12 array contains all the information
%numb is the numbers of the lines
for i=1:numb
    if flag4(i)==0
    impinge4_flag(i)=0; % The flag of 3-spheres impingement with initial value 0
    V(i)=4/3*pi*x(i,4)^3;Tri(i,:)=0;
  for j=1:numb
      
        if i~=j
           d12(i,j)=((x(i,1)-x(j,1))^2+(x(i,2)-x(j,2))^2+(x(i,3)-x(j,3))^2)^0.5;%distance between grain 1 and 2
           r12(i,j)=x(i,4)+x(j,4);%sum of the radius of grain 1 and 2
             if x(i,4)~=0
               if d12(i,j)<r12(i,j) %impingement
                  cos_beta=(x(i,4)^2+d12(i,j)^2-x(j,4)^2)/(2*x(i,4)*d12(i,j));
                  segment12(i,j)=1/3*pi*x(i,4)^3*((1-cos_beta)^2)*(2+cos_beta); %overlay volume of i due to j             
               else %no impingement
                  segment12(i,j)=0;
               end
            else segment12(i,j)=0;
            end
     for k=j:numb
         
         if j~=k&&i~=k
           d13(i,k)=((x(i,1)-x(k,1))^2+(x(i,2)-x(k,2))^2+(x(i,3)-x(k,3))^2)^0.5;%distance between grain 1 and 3
           d23(j,k)=((x(j,1)-x(k,1))^2+(x(j,2)-x(k,2))^2+(x(j,3)-x(k,3))^2)^0.5;%distance between grain 2 and 3
           r13(i,k)=x(i,4)+x(k,4);%sum of the radius of grain 1 and 3
           r23(j,k)=x(j,4)+x(k,4);%sum of the radius of grain 2 and 3
           if x(i,4)~=0
               if d12(i,j)<r12(i,j)&&d13(i,k)<r13(i,k)&&d23(j,k)<r23(j,k) %impingement
                  if d12(i,j)<max(x(i,4),x(j,4))||d13(i,k)<max(x(i,4),x(k,4))||d23(j,k)<max(x(j,4),x(k,4))
                      Tri(i,k)=0;
                  else
                  [triple overlay flag]=trioverlay(x(i,1),x(i,2),x(i,3),x(i,4),x(j,1),x(j,2),x(j,3),x(j,4),x(k,1),x(k,2),x(k,3),x(k,4));
                  Tri(i,k)=overlay;
                  Mark(i,k)=flag;
                  end
               %impinge3_flag(i)=1; % 1 indicates the grain i had impinged with two other grains
       %%%%Following is to judge if 4-spheres impingement occurs or not
       for l=k:numb
         if l~=k&&l~=j&&l~=i
           d14(i,l)=((x(i,1)-x(l,1))^2+(x(i,2)-x(l,2))^2+(x(i,3)-x(l,3))^2)^0.5;%distance between grain 1 and 4
           d24(j,l)=((x(j,1)-x(l,1))^2+(x(j,2)-x(l,2))^2+(x(j,3)-x(l,3))^2)^0.5;%distance between grain 2 and 4
           d34(k,l)=((x(k,1)-x(l,1))^2+(x(k,2)-x(l,2))^2+(x(k,3)-x(l,3))^2)^0.5;%distance between grain 3 and 4
           r14(i,l)=x(i,4)+x(l,4);%sum of the radius of grain 1 and 4
           r24(j,l)=x(j,4)+x(l,4);%sum of the radius of grain 2 and 4
           r34(k,l)=x(k,4)+x(l,4);%sum of the radius of grain 3 and 4
           if x(i,4)~=0
               if d12(i,j)<r12(i,j)&&d13(i,k)<r13(i,k)&&d23(j,k)<r23(j,k)&&d14(i,l)<r14(i,l)&&d24(j,l)<r24(j,l)&&d34(k,l)<r34(k,l) %impingement of 4 grains
                  impinge4_flag(i)=1; % 1 indicates the grain i had impinged with two other grains
               end
           end
         end
         end
         %%%%%%%%
               else Tri(i,k)=0;
               end
           else Tri(i,k)=0;
           end
            V(i)=V(i)+Tri(i,k); 
         %else
          %  V(i)=V(i);
         end   
     end
             
     else segment12(i,j)=0;
     end    
        V(i)=V(i)-abs(segment12(i,j));    
   end
      if V(i)<0
         V(i)=0;
      end
    else
        impinge4_flag(i)=1; % The flag of 3-spheres impingement with initial value 0
        V(i)=4/3*pi*x(i,4)^3;
        Tri(i,:)=0;
    end
               
end
end