function Equations = Manipulator_Dynamics(G,Gravity_Value)

%% Copying the values from the robot parameters to individual matricies
n = size(G,1);
dh = zeros(size(G,1),5);
for i = 1:size(G,1)
    if isequal(G{i,2},'Prismatic')
        dh(i,1) = 0;
    elseif isequal(G{i,2},'Revolute')
        dh(i,1) = 1;
    end
end

dh(:,2) = cell2mat(G(:,3));     % a
dh(:,3) = cell2mat(G(:,4));     % alpha
dh(:,4) = cell2mat(G(:,5));     % d
dh(:,5) = cell2mat(G(:,6));     % theta

ml = cell2mat(G(:,7));          % mass of link
mm = cell2mat(G(:,8));          % mass of motor
Il = {G{:,9}};                  % Inertia of link
Im = cell2mat(G(:,10));         % Inertia of motor
for i = 1:size(G,1)
    pl(1:3,i) = G{i,11};        % center of mass [3 x n] matrix
end
gr = cell2mat(G(:,12));         % Gear Ratio
g0 = Gravity_Value;

%% Building the Robot
for i=1:n
    if dh(i,1)==0
        L(i)=Link('a',dh(i,2),'alpha',dh(i,3),'theta',dh(i,5));
    elseif dh(i,1)==1
        L(i)=Link('a',dh(i,2),'alpha',dh(i,3),'d',dh(i,4));
    end
end
myrobot = SerialLink(L,'name','Robot');

%% Calculating the Transformation Matrix
Q=sym('q',[1 n], 'real');
Qd=sym('qd',[1 n], 'real');
syms theta a d alph
%T0E=myrobot.fkine(q)
Trans=[cos(theta) -sin(theta)*cos(alph) sin(theta)*sin(alph) a*cos(theta);...
    sin(theta) cos(theta)*cos(alph) -cos(theta)*sin(alph) sin(theta)*a;...
    0 sin(alph) cos(alph) d;0 0 0 1];
T=cell(1,n);
for i=1:n
    if dh(i,1)==1
    T{1,i}=sym2cell(simplify(subs(Trans,{theta, alph, a, d},[Q(i),dh(i,3),dh(i,2),dh(i,4)])));
    elseif dh(i,1)==0
        T{1,i}=sym2cell(simplify(subs(Trans,{theta, alph, a, d},[dh(i,5),dh(i,3),dh(i,2),Q(i)])));
    end
end

%% Storing z
z=sym('z',[3,n+1], 'real');
%z=zeros(3,n);
z(:,1)=[0;0;1];
for i=1:n
    t=cell2sym(T{1,i});
    z(:,i+1)=t(1:3,3);
end

%% Storing P
p= sym ('p', [3,n+1],'real');
p(:,1)=[0;0;0];
for i=1:n
    t=cell2sym(T{1,i});
    p(:,i+1)=t(1:3,4);
end
%% Storing Rotation
R=cell(1,n);
for i=1:n
    t=cell2sym(T{1,i});
    R{1,i}=sym2cell(t(1:3,1:3));
end

%% Calculating Jp of links
Jpl=cell(1,n);
Jpl_m= sym ('Jpl_m', [3,n]);
for i=1:n
    
    if dh(i,1)==0 %check for prismatic
        for j=1:n
            if j<=i
                Jpl_m(:,j)=z(:,j);
            elseif j>i
                Jpl_m(:,j)=0;
            end
        end
    elseif dh(i,1)==1 %check for revolute
        for j=1:n
            if j<=i
                Jpl_m(:,j)=cross(z(:,j),(pl(:,i)-p(:,j)));
            elseif j>i
                Jpl_m(:,j)=0;
            end
        end
       
    end
     Jpl{1,i}=sym2cell(Jpl_m);
end

%% Calculting Jp of motors
Jpm=cell(1,n);
Jpm_m= sym ('Jpm_m', [3,n], 'real');
for i=1:n
    
    if dh(i,1)==0 %check for prismatic
        for j=1:n
            if j<=i
                Jpm_m(:,j)=z(:,j);
            elseif j>i
                Jpm_m(:,j)=0;
            end
        end
    elseif dh(i,1)==1 %check for revolute
        for j=1:n
            if j<=i
                Jpm_m(:,j)=cross(z(:,j),(p(:,i)-p(:,j)));
            elseif j>i
                Jpm_m(:,j)=0;
            end
        end
        
    end
    Jpm{1,i}=sym2cell(Jpm_m);
end

%% Calculating Jo for links
Jol=cell(1,n);
Jol_m= sym ('Jpm_m', [3,n], 'real');
for i=1:n
    
    if dh(i,1)==0 %check for prismatic
        for j=1:n
            Jol_m(:,j)=0;
        end
    elseif dh(i,1)==1 %check for revolute
        for j=1:n
            if j<=i
                Jol_m(:,j)=z(:,j);
            elseif j>i
                Jol_m(:,j)=0;
            end
        end
        
    end
    Jol{1,i}=sym2cell(Jol_m);
end

%% Calculating Jo for motors
Jom=cell(1,n);
Jom_m= sym ('Jom_m', [3,n], 'real');
for i=1:n
    for j=1:n
        if j<i
            Jom_m(:,j)=Jol_m(:,j);
        elseif j==i
            Jom_m(:,j)=gr(i)*z(:,j);
        elseif j>i
            Jom_m(:,j)=0;
        end
    end
    Jom{1,i}=sym2cell(Jom_m);
end

%% Calculation of Mass Matrix
for i=1:n
    B=(ml(i)*cell2sym(Jpl{1,i})'*cell2sym(Jpl{1,i})+...
        cell2sym(Jol{1,i})'*cell2sym(R{1,i})*Il{1,i}*cell2sym(R{1,i})'*cell2sym(Jol{1,i})+...
        mm(i)*cell2sym(Jpm{1,i})'*cell2sym(Jpm{1,i})+...
        cell2sym(Jom{1,i})'*cell2sym(R{1,i})*Im(i)*cell2sym(R{1,i})'*cell2sym(Jom{1,i}));
end

%% Calculation of C
C = sym('C',[n], 'real');

for i=1:n
    for j=1:n
        C(i,j)=0;
        for k=1:n
            cijk=0.5*(diff(B(i,j),Q(k))+diff(B(i,k),Q(j))-diff(B(j,k),Q(i)));
            C(i,j)=C(i,j)+cijk*Qd(k);
        end
    end
end

%% Gravity matrix
g=sym('g',[n,1], 'real');
for i=1:n
    g(i)=0;
    for j=1:n
        g(i)=g(i)+(ml(j)*g0'*cell2sym(Jpl{j}(:,i))+mm(j)*g0'*cell2sym(Jpm{j}(:,i)));
    end
end

%% Final equation of motion
Qdd=sym('qdd', [n,1],'real');
%equation 
Equations = vpa(B*Qdd+C*Qd'+g, 3);
% this equation is equal to some value of torque which the user inputs so 
% when displaying make sure to display it as vpa(B*Qdd+C*Qd'+g, 3)=[ Torque vector]

end



