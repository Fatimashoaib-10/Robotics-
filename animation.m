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
q=jtraj([0,0],[0,pi/2],10);
myrobot.plot(q)
%%
