function myrobot = CreateRobot(G)

for i = 1:size(G,1)

    if isequal(G{i,2},'Prismatic')
        L(i) = Link('a',G{i,3},'alpha',G{i,4},'theta',G{i,6});
    elseif isequal(G{i,2},'Revolute')
        L(i) = Link('a',G{i,3},'alpha',G{i,4},'d',G{i,5});
    end   
end
myrobot = SerialLink(L,'name','Robot')
end