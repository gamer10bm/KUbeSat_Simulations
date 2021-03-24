function [X,Y] = OneStepMethodXVec(MethodStepHandle,FunctionHandle,...
    InitialY,XVec)
%This finds y at the end of the Xvec using XVec as the x_n points
    %The function handle should have at least one argument being either x
    %or y
        %If both, use form @(x,y)

dy = FunctionHandle;
%CurrentY is the Y_n at this point or the initial y
yn = InitialY;
%Average the difference between the elements of the Xvec for the stepsize
h = mean(diff(XVec));
for x_id = 1:length(XVec)
    xn = XVec(x_id);
    [y_n1] = MethodStepHandle(dy,h,xn,yn);
    yn = y_n1;
    xn = xn+h;
    Y(:,x_id) = yn;
    X(x_id) = xn;
end
end
