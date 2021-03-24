function plot3dshape(globalmat,spacevec_m,plotnum,fignum)
if ~exist('plotnum','var')
    plotnum = 1;
end
if ~exist('fignum','var')
    fignum = 1;
    figure(fignum); clf;
end
%% Find the YZ faces along X
edge_dings_xsweep = struct('x_id',[],'y_id',[],'z_id',[]);
%Search through each x slice
for z_id = 1:size(globalmat,3)
    slicenow = reshape(globalmat(:,:,z_id),[length(globalmat) length(globalmat)]);
    %Search each y slice of the xy mat
    for y_id = 1:size(slicenow,2)
        yvec_m = slicenow(:,y_id);
        %Determine if shape is in this slice
        if any(yvec_m==plotnum)
            %Check start and end for edges
            x_inds = [];
            x_inds_strt = [];  x_inds_end = [];
            if yvec_m(1) == plotnum
                x_inds = [1];
            end
            if yvec_m(end) == plotnum
                x_inds = [x_inds length(yvec_m)];
            end
            %Check for all starting edges
            x_indvec = 2:length(yvec_m);
            x_inds_strt = x_indvec(and(yvec_m(1:end-1)==0,yvec_m(2:end)==plotnum));
            x_indvec = 1:length(yvec_m)-1;
            x_inds_end = x_indvec(and(yvec_m(1:end-1)==plotnum,yvec_m(2:end)==0));
            %Add edge_dings
            x_inds = [x_inds x_inds_strt x_inds_end];
            for x_id = x_inds
                if isempty(edge_dings_xsweep(end).x_id)
                    edge_dings_xsweep(end) = struct('x_id',x_id,'y_id',y_id,'z_id',z_id);
                else
                    edge_dings_xsweep(end+1) = struct('x_id',x_id,'y_id',y_id,'z_id',z_id);
                end
            end
        end
    end
end

%Grab each face along x
% yzfaces = struct('xvec',[],'yvec',[],'zmat',[]);
% xvec_m = unique([edge_dings_xsweep.x_id]);
% for x = xvec_m
%     %Extract yvec and zvec for this x
%     yvec_m = unique([edge_dings_xsweep([edge_dings_xsweep.x_id]==x).y_id]);
%     zvec_m = unique([edge_dings_xsweep([edge_dings_xsweep.x_id]==x).z_id]);
%     %Determine if one  face or multiple faces in y direction
%     if any(diff(yvec_m)>1)
%         big_diff_inds_y = find(diff(yvec_m)>1);
%     end
%     if any(diff(zvec_m)>1)
%         big_diff_inds_z = find(diff(zvec_m)>1);
%     end
% end
        


%Try plotting edges
if 1
    xyztable = reshape(cell2mat(struct2cell(edge_dings_xsweep)),[3 length(edge_dings_xsweep)]);
    xvec_m = spacevec_m(xyztable(1,:)); 
    yvec_m = spacevec_m(xyztable(2,:)); 
    zvec_m = spacevec_m(xyztable(3,:));
    figure(fignum)
    scatter3(xvec_m,yvec_m,zvec_m,'filled')
    hold on
    xlim([0 spacevec_m(end)]); ylim([0 spacevec_m(end)]); zlim([0 spacevec_m(end)]);
%     keyboard
end

%% Find XY faces along z
edge_dings_zsweep = struct('x_id',[],'y_id',[],'z_id',[]);
%Search through each z slice
for x_id = 1:size(globalmat,1)
    slicenow = reshape(globalmat(x_id,:,:),[length(globalmat) length(globalmat)]);
    %Search each y slice of the xy mat
    for y_id = 1:size(slicenow,2)
        yvec_m = slicenow(y_id,:);
        %Determine if shape is in this slice
        if any(yvec_m==plotnum)
            %Check start and end for edges
            z_inds = [];
            z_inds_strt = [];  z_inds_end = [];
            if yvec_m(1) == plotnum
                z_inds = [1];
            end
            if yvec_m(end) == plotnum
                z_inds = [z_inds length(yvec_m)];
            end
            %Check for all starting edges
            x_indvec = 2:length(yvec_m);
            z_inds_strt = x_indvec(and(yvec_m(1:end-1)==0,yvec_m(2:end)==plotnum));
            x_indvec = 1:length(yvec_m)-1;
            z_inds_end = x_indvec(and(yvec_m(1:end-1)==plotnum,yvec_m(2:end)==0));
            %Add edge_dings
            z_inds = [z_inds z_inds_strt z_inds_end];
            for z_id = z_inds
                if isempty(edge_dings_zsweep(end).x_id)
                    edge_dings_zsweep(end) = struct('x_id',x_id,'y_id',y_id,'z_id',z_id);
                else
                    edge_dings_zsweep(end+1) = struct('x_id',x_id,'y_id',y_id,'z_id',z_id);
                end
            end
        end
    end
end

%Try plotting edges
if 1
    xyztable = reshape(cell2mat(struct2cell(edge_dings_zsweep)),[3 length(edge_dings_zsweep)]);
    xvec_m = spacevec_m(xyztable(1,:)); 
    yvec_m = spacevec_m(xyztable(2,:)); 
    zvec_m = spacevec_m(xyztable(3,:));
    figure(fignum)
    scatter3(xvec_m,yvec_m,zvec_m,'filled')
    hold on
    xlim([0 spacevec_m(end)]); ylim([0 spacevec_m(end)]); zlim([0 spacevec_m(end)]);
%     keyboard
end

%% Find the XZ faces along y
edge_dings_ysweep = struct('x_id',[],'y_id',[],'z_id',[]);
%Search through each y slice
for z_id = 1:size(globalmat,2)
    slicenow = reshape(globalmat(:,:,z_id),[length(globalmat) length(globalmat)]);
    %Search each y slice of the xy mat
    for x_id = 1:size(slicenow,2)
        xvec_m = slicenow(x_id,:);
        %Determine if shape is in this slice
        if any(xvec_m==plotnum)
            %Check start and end for edges
            y_inds = [];
            y_inds_strt = [];  y_inds_end = [];
            if xvec_m(1) == plotnum
                y_inds = [1];
            end
            if xvec_m(end) == plotnum
                y_inds = [y_inds length(xvec_m)];
            end
            %Check for all starting edges
            y_indvec = 2:length(xvec_m);
            y_inds_strt = y_indvec(and(xvec_m(1:end-1)==0,xvec_m(2:end)==plotnum));
            y_indvec = 1:length(xvec_m)-1;
            y_inds_end = y_indvec(and(xvec_m(1:end-1)==plotnum,xvec_m(2:end)==0));
            %Add edge_dings
            y_inds = [y_inds y_inds_strt y_inds_end];
            for y_id = y_inds
                if isempty(edge_dings_ysweep(end).x_id)
                    edge_dings_ysweep(end) = struct('x_id',x_id,'y_id',y_id,'z_id',z_id);
                else
                    edge_dings_ysweep(end+1) = struct('x_id',x_id,'y_id',y_id,'z_id',z_id);
                end
            end
        end
    end
end

%Try plotting edges
if 1
    xyztable = reshape(cell2mat(struct2cell(edge_dings_ysweep)),[3 length(edge_dings_ysweep)]);
    xvec_m = spacevec_m(xyztable(1,:)); 
    yvec_m = spacevec_m(xyztable(2,:)); 
    zvec_m = spacevec_m(xyztable(3,:));
    figure(fignum)
    scatter3(xvec_m,yvec_m,zvec_m,'filled')
    xlim([0 spacevec_m(end)]); ylim([0 spacevec_m(end)]); zlim([0 spacevec_m(end)]);
%     keyboard
end

%Use corners to find edges

%Use edges to make surfaces

%Plot surfaces
end