function thetaT = transformParams(theta)

thetaT = NaN(size(theta));
thetaT(1,:) = log(theta(1,:)./(1-theta(1,:)));
%thetaT(1,:) = log((theta(1,:)+1)./(1-theta(1,:)));
thetaT(2,:) = log(theta(2,:)./(1-theta(2,:)));
thetaT(3,:) = log(theta(3,:));
if size(theta,1) > 3
   thetaT(4,:) = log(theta(4,:));
end

end