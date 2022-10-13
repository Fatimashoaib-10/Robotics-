function singout = Singularity(G)

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

%% Building the robot
clc
for i=1:n
    if dh(i,1)==0
        L(i)=Link('a',dh(i,2),'alpha',dh(i,3),'theta',dh(i,5));
    elseif dh(i,1)==1
        L(i)=Link('a',dh(i,2),'alpha',dh(i,3),'d',dh(i,4));
    end
end
myrobot = SerialLink(L,'name','Robot');

%% Calculating transformation matrix
clc
Q=sym('q',[1 n], 'real');
Qd=sym('qd',[1 n], 'real');
syms a alph d theta
%T0E=myrobot.fkine(q)
Trans=[cos(theta) -sin(theta)*cos(alph) sin(theta)*sin(alph) a*cos(theta);...
    sin(theta) cos(theta)*cos(alph) -cos(theta)*sin(alph) sin(theta)*a;...
    0 sin(alph) cos(alph) d;0 0 0 1];
T=cell(1,n);
for i=1:n
    T{1,i}=sym2cell(simplify(subs(Trans,{theta, alph, a, d},[Q(i),dh(i,3),dh(i,2),dh(i,4)])));
end

%% Storing z
z=sym('z',[3,n+1], 'real');
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

v=det(Ja);
eqn= v==0;
for i=1:n
    fprintf('q(%d)=',i);
    singout(i) = double(solve(v,Q(1,i)))
end 
end
