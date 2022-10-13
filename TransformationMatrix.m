function T = TransformationMatrix(a, alpha, d, theta)

c1 = cosd(theta);
c2 = cosd(alpha);
s1 = sind(theta);
s2 = sind(alpha);
T = [c1 (-s1)*c2 s1*s2 a*c1; s1 c1*c2 (-c1)*s2 a*s1; 0 s2 c2 d; 0 0 0 1];

end

