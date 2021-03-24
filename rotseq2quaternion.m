function [quatvec_wxyz] = rotseq2quaternion(rotvec_rad,seqvec_xyz,roll_pitch_yaw_vec_rad)
% [quatvec_xyzw] = rotseq2quaternion(rotvec_rad,seqvec_xyz,roll_pitch_yaw_vec_rad)
%
% Given a vector of rotations with a sequence of corresponding axes of
% rotation, calculate the resulting quaternion
%
%
%

if ~exist('roll_pitch_yaw_vec_rad','var')
    roll_pitch_yaw_vec_rad = [0,0,0];
end

if length(rotvec_rad) ~= length(seqvec_xyz)
    warning('First two inputs need to be the same length')
end
%Initialize matrix functions
rotx = @(roll)[1, 0, 0; 0, cos(-roll), sin(-roll); 0, -sin(-roll), cos(-roll)];
roty = @(pitch)[cos(-pitch), 0, -sin(-pitch); 0, 1, 0; sin(-pitch), 0, cos(-pitch)];
rotz = @(yaw) [cos(-yaw), sin(-yaw), 0; -sin(-yaw), cos(-yaw), 0; 0, 0, 1];
%Do initial rotation based roll_pitch_yaw_vec
roll = roll_pitch_yaw_vec_rad(1);
pitch = roll_pitch_yaw_vec_rad(2);
yaw = roll_pitch_yaw_vec_rad(3);
% Find transformation components in ZYX sequence
C1 = rotx(roll);
C2 = roty(pitch);
C3 = rotz(yaw);

% complete transformation from inertial to body frame
rotmat1 = C3*C2*C1;

%Do the entire sequence of rotations
% rotmat1 = eye(3);
for r_id = length(rotvec_rad):-1:1
    switch seqvec_xyz(r_id)
        case 'x'
            rotmat1 = rotx(rotvec_rad(r_id))*rotmat1;
        case 'y'
            rotmat1 = roty(rotvec_rad(r_id))*rotmat1;
        case 'z'
            rotmat1 = rotz(rotvec_rad(r_id))*rotmat1;
    end
end
%Calculate resulting quaternion
% Using method in https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
m11 = rotmat1(1,1); m12 = rotmat1(1,2); m13 = rotmat1(1,3);
m21 = rotmat1(2,1); m22 = rotmat1(2,2);m23=rotmat1(2,3);
m31 = rotmat1(3,1); m32 = rotmat1(3,2); m33 =rotmat1(3,3);
trot = trace(rotmat1);
if trot>0
    s = sqrt(trot+1)*2;
    qw = 0.25*s; qx = (m32-m23)/s; qy = (m13-m31)/s;
    qz = (m21-m12)/s;
elseif m11>m22 && m11>m33
    s = sqrt(1+m11-m22-m33)*2;
    qw = (m32-m23)/s; qx = 0.25*s; qy = (m12+m21)/s;
    qz = (m13+m31)/s;
elseif m22 > m33
    s = sqrt(1+m22-m11-m33)*2;
    qw=(m13-m31)/s; qx=(m13+m31)/s; qy=0.25*s;
    qz =(m23+m32)/s;
else
    s = sqrt(1+m33-m11-m22)*2;
    qw=(m21-m12)/s;
    qx=(m13-m31)/s;
    qy=(m23+m32)/s;
    qz=0.25*s;
end

quatvec_wxyz = [qw;qx;qy;qz];
end