function theta = invTransformParams(thetaT)

theta = NaN(size(thetaT));
theta(1,:) = exp(thetaT(1,:))./(1+exp(thetaT(1,:)));
%theta(1,:) = -1 + 2*exp(thetaT(1,:))./(1+exp(thetaT(1,:)));
if any(thetaT(1,:)>500)
    theta(1,thetaT(1,:)>500) = 1;
end
theta(2,:) = exp(thetaT(2,:))./(1+exp(thetaT(2,:)));
if any(thetaT(2,:)>500)
    theta(2,thetaT(2,:)>500) = 1;
end

theta(3,:) = exp(thetaT(3,:));
if size(thetaT,1) > 3
    theta(4,:) = exp(thetaT(4,:));
end

end