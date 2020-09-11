% This function can be only used when each pair of the spheres among the three spheres are intersected
% with each other
function [V_ABC overlay flag]=trioverlay(x1,y1,z1,alpha,x2,y2,z2,belta,x3,y3,z3,gamma)
%V_ABC is the triple intersection volume;
%overlay is the volume of partial triple intersection contributed by sphere-1
%This function is a combination of function_trilateration and
%function_contriA
V_ABC=0;overlay=0;flag=0;
if alpha==0||belta==0||gamma==0
    V_ABC=0;
else if alpha>0&&belta>0&&gamma>0
a=sqrt((x2-x3)^2+(y2-y3)^2+(z2-z3)^2);%distance of sphere B and C
b=sqrt((x1-x3)^2+(y1-y3)^2+(z1-z3)^2);%distance of sphere A and C
c=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);%distance of sphere A and B
e1=(belta^2-gamma^2)/a^2;
e2=(gamma^2-alpha^2)/b^2;
e3=(alpha^2-belta^2)/c^2;
w2=(alpha^2*a^2+belta^2*b^2+gamma^2*c^2)*(a^2+b^2+c^2)-2*(alpha^2*a^4+belta^2*b^4+gamma^2*c^4)+a^2*b^2*c^2*(e1*e2+e2*e3+e3*e1-1);%w2 is w^2
q1=a*(b^2+c^2-a^2+belta^2+gamma^2-2*alpha^2+e1*(b^2-c^2));
q2=b*(c^2+a^2-b^2+gamma^2+alpha^2-2*belta^2+e2*(c^2-a^2));
q3=c*(a^2+b^2-c^2+alpha^2+belta^2-2*gamma^2+e3*(a^2-b^2));
cos_AB=(alpha^2+c^2-belta^2)/(2*alpha*c);%cosine of intersection angle between A and B
cos_BA=(belta^2+c^2-alpha^2)/(2*belta*c);%cosine of intersection angle between B and A
AB=1/3*pi*alpha^3*(1-cos_AB)^2*(2+cos_AB);
BA=1/3*pi*belta^3*(1-cos_BA)^2*(2+cos_BA);
V_AB=AB+BA;%Intersection volume of A and B

cos_AC=(alpha^2+b^2-gamma^2)/(2*alpha*b);%cosine of intersection angle between A and C
cos_CA=(gamma^2+b^2-alpha^2)/(2*gamma*b);%cosine of intersection angle between C and A
AC=1/3*pi*alpha^3*(1-cos_AC)^2*(2+cos_AC);
CA=1/3*pi*gamma^3*(1-cos_CA)^2*(2+cos_CA);
V_AC=AC+CA;%Intersection volume of A and C

cos_BC=(belta^2+a^2-gamma^2)/(2*belta*a);%cosine of intersection angle between B and C
cos_CB=(gamma^2+a^2-belta^2)/(2*gamma*a);%cosine of intersection angle between C and B
BC=1/3*pi*belta^3*(1-cos_BC)^2*(2+cos_BC);
CB=1/3*pi*gamma^3*(1-cos_CB)^2*(2+cos_CB);
V_BC=BC+CB;%Intersection volume of B and C

if w2>0  %three spheres intersect in two common points
w=sqrt(w2);
V1=w/12;%volume of the tetrahedron PABC
%{
V2=1/3*alpha^3*(atan(2*w/q2)+atan(2*w/q3)-atan(b*w/alpha/q2*(1-e2))-atan(c*w/alpha/q3*(1+e3)));%volume of the overlay part of A
V3=1/3*belta^3*(atan(2*w/q1)+atan(2*w/q3)-atan(a*w/belta/q1*(1-e1))-atan(c*w/belta/q3*(1+e3)));%volume of the overlay part of B
V4=1/3*gamma^3*(atan(2*w/q1)+atan(2*w/q2)-atan(a*w/gamma/q1*(1-e1))-atan(b*w/gamma/q2*(1+e2)));%volume of the overlay part of C
V5=(1/2/pi)*V_AB*atan(2*w/q3);%volume of the overlay part of A and B encolosed
V6=(1/2/pi)*V_AC*atan(2*w/q2);%volume of the overlay part of A and C encolosed
V7=(1/2/pi)*V_BC*atan(2*w/q1);%volume of the overlay part of B and C encolosed
%}
at(1)=atan(2*w/q1);
at(2)=atan(2*w/q2);
at(3)=atan(2*w/q3);
at(4)=atan(b*w/alpha/q2*(1-e2));
at(5)=atan(c*w/alpha/q3*(1+e3));
at(6)=atan(c*w/belta/q3*(1-e3));
at(7)=atan(a*w/belta/q1*(1+e1));
at(8)=atan(a*w/gamma/q1*(1-e1));
at(9)=atan(b*w/gamma/q2*(1+e2));
for i=1:9
    if at(i)<0
        at(i)=at(i)+pi;
    else at(i)=at(i);
    end
end
V2=a/4*(belta^2+gamma^2-a^2*(1/6-e1^2/2))*at(1);%volume of the overlay part of A
V3=b/4*(gamma^2+alpha^2-b^2*(1/6-e2^2/2))*at(2);%volume of the overlay part of B
V4=c/4*(alpha^2+belta^2-c^2*(1/6-e3^2/2))*at(3);%volume of the overlay part of C
V5=1/3*alpha^3*(at(4)+at(5));%volume of the overlay part of A and B encolosed
V6=1/3*belta^3*(at(6)+at(7));%volume of the overlay part of A and C encolosed
V7=1/3*gamma^3*(at(8)+at(9));%volume of the overlay part of B and C encolosed
V_ABC=2*(V1-V2-V3-V4+V5+V6+V7);
V=V_AB+V_AC+V_BC-2*V_ABC;

%Following is to calculate the half-value of the distance between two
%common points shared by the three spheres. Only under this contidition, it
%is necessary to calculate this value

%radius of the circle for the pairwise intersection
r12=sqrt(abs((2*(alpha^2*belta^2+alpha^2*c^2+belta^2*c^2)-(alpha^4+belta^4+c^4))/(4*c^2)));
r21=r12;
r13=sqrt(abs((2*(alpha^2*gamma^2+alpha^2*b^2+gamma^2*b^2)-(alpha^4+gamma^4+b^4))/(4*b^2)));
r31=r13;
r23=sqrt(abs((2*(belta^2*gamma^2+belta^2*a^2+gamma^2*a^2)-(belta^4+gamma^4+a^4))/(4*a^2)));
r32=r23;

%eg.b>=max(alpha,belta),d=d12,otherwise d=-d12
%distance between the plane of the rij-disc and the center of the i-sphere
   d12=abs((alpha^2-belta^2+c^2)/(2*c));
   d21=c-d12;
   d13=abs((alpha^2-gamma^2+b^2)/(2*b));
   d31=b-d13;
   d23=abs((belta^2-gamma^2+a^2)/(2*a));
   d32=a-d23;

%half-length of the distance between the two common points
H1=((alpha^2*belta^2+gamma^2*c^2)*(a^2+b^2-c^2)+(alpha^2*gamma^2+belta^2*b^2)*(a^2+c^2-b^2))/(2*(c^2*b^2+c^2*a^2+b^2*a^2)-(a^4+b^4+c^4));
H2=((belta^2*gamma^2+alpha^2*a^2)*(c^2+b^2-a^2)-(alpha^2*a^4+belta^2*b^2+gamma^2*c^4)-a^2*b^2*c^2)/(2*(c^2*b^2+c^2*a^2+b^2*a^2)-(a^4+b^4+c^4));
H=sqrt(abs(H1+H2));

overlay=(3*alpha-H)/(3*(alpha+belta+gamma)-3*H)*V_ABC;%overlay volume contributed by sphere-1
flag=1;
 else if w2==0   %three spheres intersect in a single point
        V_ABC=0;
        V=V_AB+V_AC+V_BC;
        overlay=0;
        flag=2;
    else if w2<0   %three spheres have no common point
        t=sqrt((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c));%quantity t related to the distance 
        t1=sqrt((a+belta+gamma)*(-a+belta+gamma)*(a-belta+gamma)*(a+belta-gamma));%quantity t-alpha,belta,gamma
        t2=sqrt((alpha+b+gamma)*(-b+alpha+gamma)*(b-alpha+gamma)*(b+alpha-gamma));%quantity t-belta,gamma,alpha
        t3=sqrt((alpha+belta+c)*(-c+alpha+belta)*(c-alpha+belta)*(c+alpha-belta));%quantity t-gamma,alpha,belta
        
        p1=((b^2-c^2+belta^2-gamma^2)^2+(t-t1)^2)/(4*a^2)-alpha^2;
        p2=((b^2-c^2+belta^2-gamma^2)^2+(t+t1)^2)/(4*a^2)-alpha^2;
        
        p3=((a^2-c^2+alpha^2-gamma^2)^2+(t-t2)^2)/(4*b^2)-belta^2;
        p4=((a^2-c^2+alpha^2-gamma^2)^2+(t+t2)^2)/(4*b^2)-belta^2;
        
        p5=((a^2-b^2+alpha^2-belta^2)^2+(t-t3)^2)/(4*c^2)-gamma^2;
        p6=((a^2-b^2+alpha^2-belta^2)^2+(t+t3)^2)/(4*c^2)-gamma^2;
         if p1>0&&p2>0&&p3>0&&p4>0&&p5>0&&p6>0   %each circle of pairwise intersection lies outside the third sphere
            V_ABC=0;
            V=V_AB+V_AC+V_BC;
            overlay=0;
            flag=3;
         %Following discussing the third sphere contains the intersection of
         %the other two spheres      
         else if p1<0&&p2<0&&p3>0&&p4>0&&p5>0&&p6>0  %the sphere centered at A contains the pairwise intersection of the other two spheres
            V_ABC=V_BC;
            V=V_AB+V_AC+V_BC-2*V_ABC;
            overlay=0;
            flag=4;
         else if p1>0&&p2>0&&p3<0&&p4<0&&p5>0&&p6>0  %the sphere centered at B contains the pairwise intersection of the other two spheres
            V_ABC=V_AC;
            V=V_AB+V_AC+V_BC-2*V_ABC;
            overlay=AC;
            flag=5;
         else if p1>0&&p2>0&&p3>0&&p4>0&&p5<0&&p6<0  %the sphere centered at C contains the pairwise intersection of the other two spheres
            V_ABC=V_AB;
            V=V_AB+V_AC+V_BC-2*V_ABC;
            overlay=AB;
            flag=6;
         %Following discussing the center of the third sphere locates in the intersection of
         %the other two spheres  
         else if p1>0&&p2>0&&p3<0&&p4<0&&p5<0&&p6<0  %the sphere centered at A located in the pairwise intersection of B and C
            V_ABC=V_AB+V_AC-4/3*pi*alpha^3;
            V=V_AB+V_AC+V_BC-2*V_ABC;
            overlay=V_ABC-BA-CA;
            flag=7;
         else if p1<0&&p2<0&&p3>0&&p4>0&&p5<0&&p6<0  %the sphere centered at B located in the pairwise intersection of A and C
            V_ABC=V_AB+V_BC-4/3*pi*belta^3;
            V=V_AB+V_AC+V_BC-2*V_ABC;
            overlay=AB;
            flag=8;
         else if p1<0&&p2<0&&p3<0&&p4<0&&p5>0&&p6>0  %the sphere centered at C located in the pairwise intersection of B and A
            V_ABC=V_AC+V_BC-4/3*pi*gamma^3;
            V=V_AB+V_AC+V_BC-2*V_ABC;
            overlay=AC;
            flag=9;
             end
             end
             end
            end
            end
            end
        end
        end
     end   
  end
end
    end
end
   