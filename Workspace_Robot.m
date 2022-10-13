clc
clear all
% a1 = 0;     alpha1 = pi/2;   d1 = 0.0345;
% a2 = 0;     alpha2 = -pi/2;  d2 = 0;
% a3 =-0.045; alpha3 = pi/2;   d3 = 0.55;
% a4 = 0.045; alpha4 = -pi/2;  d4= 0;
% a5 = 0;   alpha5 = -pi/2;  d5 = 0.30;
% a6 = 0;   alpha6 = pi/2;   d6 = 0;
% a7 = 0;     alpha7 = 0;      d7 =0.060;
% n=7;
%% Enter the robot
n=input("Enter number of links: ");
dh=zeros(n,5); %1=a 2=alp 3=d 4=theta

for i=1:n
    fprintf("Entering values for link %d\n",i);
    dh(i,2) = input("Enter value for a: ");
    dh(i,3) = input("Enter value for alpha: ");
    dh(i,1) = input("Enter 0 for prismatic or enter 1 for revolute: ");
    if dh(i,1) == 0
        dh(i,5) = input("enter value for theta: ");
    elseif dh(i,1) == 1
        dh(i,4) = input("enter value for d: ");
    end
end
%% Building the robot
clc
for i=1:n
    if dh(i,1)==0
        L(i)=Link('a',dh(i,2),'alpha',dh(i,3),'theta',dh(i,5))
    elseif dh(i,1)==1
        L(i)=Link('a',dh(i,2),'alpha',dh(i,3),'d',dh(i,4))
    end
end
myrobot=SerialLink(L,'name','Robot')

%%
q=sym('q',[1 n], 'real')
T0E=vpa(myrobot.fkine(q'),3)

%%
t=zeros(n,2);
for i=1:n
   t(i,1)=input("min");
   t(i,2)=input("max");
end   
%%
T=zeros(500,1);
N = 500;
for i=1:n 
T(1:500,i) = t(i,1) + (t(i,2)-t(i,1))*rand(N,1);
end
%%
for i = 1:N
% A1 = TransMat(a1,alpha1,d1,t1(i));
% A2 = TransMat(a2,alpha2,d2,t2(i));
% A3 = TransMat(a3,alpha3,d3,t3(i));
% A4 = TransMat(a4,alpha4,d4,t4(i));
% A5 = TransMat(a5,alpha5,d5,t5(i));
% A6 = TransMat(a6,alpha6,d6,t6(i));
% A7 = TransMat(a7,alpha7,d7,t7(i));
T_0E=vpa(subs(T0E,q,T(i,:)),3);
 X=T_0E(1,4);Y=T_0E(2,4);Z=T_0E(3,4);
 plot3(X,Y,Z,'b.')
 hold on;
i
end
