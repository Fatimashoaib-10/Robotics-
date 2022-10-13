clc
clear all
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

%% Calculating transformation matrix
clc
Q=sym('q',[1 n], 'real')
Qd=sym('qd',[1 n], 'real')
syms theta alpha a d
%T0E=myrobot.fkine(q)
Trans=[cos(theta) -sin(theta)*cos(alpha) sin(theta)*sin(alpha) a*cos(theta);...
    sin(theta) cos(theta)*cos(alpha) -cos(theta)*sin(alpha) sin(theta)*a;...
    0 sin(alpha) cos(alpha) d;0 0 0 1]
T=cell(1,n);
for i=1:n
    T{1,i}=sym2cell(simplify(subs(Trans,{theta, alpha, a, d},[Q(i),dh(i,3),dh(i,2),dh(i,4)])))
end
%% Storing z
z=sym('z',[3,n+1], 'real');
z(:,1)=[0;0;1];
for i=1:n
    t=cell2sym(T{1,i});
    z(:,i+1)=t(1:3,3);
end
%% Storing P
p= sym ('p', [3,n+1],'real')
p(:,1)=[0;0;0];
for i=1:n
    t=cell2sym(T{1,i});
    p(:,i+1)=t(1:3,4);
end
%% Calculating the Jacobian (position)
Jp=sym ('Jp',[3 n], 'real');
for i=1:n
    
    if dh(i,1)==0 %check for prismatic
        Jp(:,i)=z(:,i);
    elseif dh(i,1)==1 %check for revolute
        Jp(:,i)=cross(z(:,i),(p(:,n)-p(:,i)));
    end
end 
%% Calculating the Jacobian (orientation)
Jo=sym('Jo', [3, n], 'real');
for i=1:n
    
    if dh(i,1)==0 %check for prismatic
        Jo(:,i)=0;
    elseif dh(i,1)==1 %check for revolute
        Jo(:,i)=z(:,i);
    end
end 
%% Check for singularity
J=[Jp;Jo]; % merging Jp and Jo to obtain a 6xn matrix 
Ja=licols(J');
Ja=Ja'; % the square matrix
%% 
T0E=myrobot.fkine(Q'); %Calculating the e.e transformation and otaining the x,y,z positions as well as the orientation of the end effector 
R0E=T0E(1:3,1:3);
%the below 3 phi's are the orientation w.r.t all the 3 axes
phi_x=atan2(R0E(3,1),R0E(3,3));
phi_y=atan2(-R0E(3,1),sqrt(R0E(3,2)^2+R0E(3,3)^2));
phi_z=atan2(R0E(2,1),R0E(1,1));
xe=[T0E(1:3,4);phi_x;phi_y;phi_z];

%%

if n<6

pd= sym('pd',[n,1]);
pd_dot=sym('pd_dot',[n,1]);
    for i=1:n
     % the first three elements are the x,y,z positions, the second three are the x,y,z orientaion
     str=input("Enter the End effector position/oriention :",'s');
     f=inline(str,'t');
     pd(i,1)=inline2sym(f);
    end
    for i=1:n
         % the first three elements are the x,y,z positions, the second three are the x,y,z orientaion
    str=input("Enter the End effector position/oriention :",'s');
     f=inline(str,'t');
     pd_dot(i,1)=inline2sym(f); %this is the E.E velocity 
    end
    % for robots with more than 6D0F, the maximum number of E.E element can
    % always be 6.(First three elements are the position, the second three
    % elements are the position)
 elseif n>=6
pd= sym('pd',[6,1]);
pd_dot=sym('pd_dot',[6,1]);
     for i=1:6
         % the first three elements are the x,y,z positions, the second three are the x,y,z orientaion
    str=input("Enter the End effector position/oriention :",'s');
     f=inline(str,'t');
     pd(i,1)=inline2sym(f);
     end
     for i=1:6
         % the first three elements are the x,y,z positions, the second three are the x,y,z orientaion
     str=input("Enter the End effector position/oriention :",'s');
     f=inline(str,'t');
     pd_dot(i,1)=inline2sym(f);
     end
 end 
 %%
 t=input("Enter the simulation time: "); %% user needs to input this time 
 t=0:1:t;
 pd=subs(pd,'t',t);
 pd_dot=subs(pd_dot,'t',t);
 
 %% calculating the q and q dot 
q= zeros(n, length(t));
q(:,1) = [pi -pi/2 4]; %ask the user for initial conditions 
if n<=6
K=1*eye(n,n);
x_e=zeros(n,length(t));
xe=xe(1:n)
 
for i = 1:length(t)
    x_e(:,i)=subs(xe,Q,q(:,i)')
    J_a = subs(Ja,Q,q(:,i)')
    e(:,i)=pd(:,i) - x_e(:,i)
    xd_dot=pd_dot(:,i)
    q_dot = inv(J_a)*(xd_dot+K*e(:,i));
    q(:,i+1)=q(:,i)+q_dot*0.001;
end

elseif n>6
K=1*eye(6,6);
x_e=zeros(n,length(t));
for i = 1:length(t)
    x_e(:,i)=subs(xe,Q',q(:,i));
    J_a = subs(J,Q',q(:,i));
    e(:,i)=pd(:,i) - x_e(:,i);
    xd_dot=pd_dot(:,i);
    q_dot = pinv(J_a)*(xd_dot+K*e(:,i));
    q(:,i+1)=q(:,i)+q_dot*0.001;
end
end 

