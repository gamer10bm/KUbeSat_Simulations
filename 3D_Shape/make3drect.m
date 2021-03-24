function [globalmat, spacevec_m] = make3drect(xyz_sizevec_m, xyz_strtvec_m, globalmax_m, global_N)
% [globalmat, spacevec_m] = make3drect(xyz_sizevec_m, xyz_strtvec_m, globalmax_m, global_N)
% 
% This function builds a 3d logic matrix that can be used to represent a
% rectangular prism.
%
% Inputs:   xyz_sizevec_m:  Vector of prism x y and z sizes in meters
%               xyz_strtvec_m:  Vector of prism start dimensions in meters
%               globalmax_m:      Maximum global space dimension in meters
%               global_N:            Length of one side of cubic globalmat
%
% Outputs:      globalmat:     Global matrix with rectangular prism within
%                    spacevec_m:  Vector of meters along each globalmat side
%
% Created by Bailey Miller 3/13
%
% See also plot3dshape.m
%% Input Check
if ~exist('globalmax_m','var')
    maxsize = max(xyz_sizevec_m);
    maxstrt = max(xyz_strtvec_m);
    globalmax_m = maxstrt*2+maxsize;
end
if ~exist('global_N','var')
    dim_scale = globalmax_m;
    dN_to_m = 0.05; %meters per step
    while dim_scale < 1
        dN_to_m = dN_to_m/10;
        dim_scale = dim_scale*10;
    end
    global_N = floor(globalmax_m/dN_to_m);
else
    dN_to_m = globalmax_m/global_N;
end
%% Convert physical units to global index vectors
xstrt_m = xyz_strtvec_m(1); ystrt_m = xyz_strtvec_m(2); zstrt_m = xyz_strtvec_m(3);
xsize_m = xyz_sizevec_m(1); ysize_m = xyz_sizevec_m(2); zsize_m = xyz_sizevec_m(3);
% Get starting position in terms of global index
xstrt_global = floor(xstrt_m/dN_to_m)+1;
ystrt_global = floor(ystrt_m/dN_to_m)+1;
zstrt_global = floor(zstrt_m/dN_to_m)+1;
if xstrt_m == 0
    xstrt_global = 1;
end
if ystrt_m == 0
    ystrt_global = 1;
end
if zstrt_m ==0
    zstrt_global = 1;
end

%Get size in terms of global index
xsize_global = ceil(xsize_m/dN_to_m);
ysize_global = ceil(ysize_m/dN_to_m);
zsize_global = ceil(zsize_m/dN_to_m);

%Make vectors of indices
xinds = xstrt_global:xstrt_global+xsize_global;
yinds = ystrt_global:ystrt_global+ysize_global;
zinds = zstrt_global:zstrt_global+zsize_global;
% if xinds(end) > global_N+1
%     xinds = xstrt_global:global_N+1;
% end
% if yinds(end) > global_N+1
%     yinds = ystrt_global:global_N+1;
% end
% if zinds(end) > global_N+1
%     zinds = zstrt_global:global_N+1;
% end

%Load global shape matrix
globalmat = zeros(global_N+2,global_N+2,global_N+2);
globalmat(xinds, yinds, zinds) = 1;

%Generate x y and z vectors of physical space
spacevec_m = (0:global_N+1)*dN_to_m;
end