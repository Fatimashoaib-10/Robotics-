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
end;
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
J=[Jp;Jo]
Ja=licols(J');
Ja=Ja' % the square matrix
%%
v=det(Ja);
eqn= v==0;
for i=1:n
    fprintf('q(%d)=',i);
    sol(i) = solve(v,Q(1,i))
end 
