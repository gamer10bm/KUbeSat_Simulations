function [roll, pitch, yaw, Rotmat] = quaternion2euler(quatvec)
% [roll, pitch, yaw, Rotmat] = quaternion2euler(quatvec)
%
% This function takes the input of a normalized quaternion vector and
% calculates the euler angles and rotation matrix
%
% Inputs:   quatvec - 4x1 Vector of normalized quaternion
%
% Outputs: roll -      Roll angle in degrees
%               pitch -    Pitch angle in degrees
%               yaw -      Yaw angle in degrees
%               Rotmat - 3x3 rotation matrix for rotating a 3x1 vector
%                                  based on quatvec
%
% Created by Bailey Miller 3/14/2021
% See also: make3drect.m, extract3dshape.m
%

%% Compile quaternion rotation matrix
q0 = quatvec(1); q1 = quatvec(2); q2 = quatvec(3); q3 = quatvec(4);
%Using https://thepoorengineer.com/en/quaternion/ 
Rotmat = [q0^2+q1^2-q2^2-q3^2, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2);
    2*(q1*q2+q0*q3), q0^2-q1^2+q2^2-q3^2, 2*(q2*q3-q0*q1);
    2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), q0^2-q1^2-q2^2+q3^2];

%Find roll (beta)
c20 = Rotmat(3,1);
%Check for different cases
if -c20 > 0.99999
    yaw = 0;
    pitch = pi/2;
    roll = atan2(Rotmat(1,2),Rotmat(1,3));
elseif c20 > 0.9999
    yaw = 0;
    pitch = -pi/2;
    roll = atan2(-Rotmat(1,2),-Rotmat(1,3));
else
    yaw = atan2(Rotmat(2,1),Rotmat(1,1));
    pitch = asin(-Rotmat(3,1));
    roll = atan2(Rotmat(3,2),Rotmat(3,3));
end
%Convert to degrees
roll = rad2deg(roll);
pitch = rad2deg(pitch);
yaw = rad2deg(yaw);
end