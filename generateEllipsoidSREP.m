function generateEllipsoidSREP( nRows, nCols, dataDir, flag )
close all;
fid = fopen(fullfile(dataDir,'radii.txt'),'rt');
radii = fscanf(fid, '%f,%f,%f\n');
radii_mat = reshape(radii, 3, []);
idx = size(radii_mat,2);
current_radii = radii_mat(:,idx);
current_radii_sorted = sort(current_radii);
%%
mesh_prefix = sprintf('%04d',idx);
mesh_deformed = read_vtkMesh(fullfile(dataDir, strcat(mesh_prefix,'.vtk')));
mesh_deformed_cog = mean(mesh_deformed.points);
mesh_deformed_points_centered = mesh_deformed.points - repmat(mesh_deformed_cog, size(mesh_deformed.points) ./ size(mesh_deformed_cog));
mesh_secondMoment = mesh_deformed_points_centered'*mesh_deformed_points_centered;
[V,~] = eigs(mesh_secondMoment);
%%
nCrestPoints = nRows*2 + (nCols-2)*2;
rx = current_radii_sorted(3);
ry = current_radii_sorted(2);
rz = current_radii_sorted(1);

ELLIPSE_SCALE = 0.90;

mrx_o = (rx*rx - rz*rz) / rx;
mrx = mrx_o * ELLIPSE_SCALE;
% mrx = rx - 1/ELLIPSE_SCALE * rz * rz / rx;
mry_o = (ry*ry - rz*rz) / ry;
mry = mry_o * ELLIPSE_SCALE;
% mry = ry - 1/ELLIPSE_SCALE * rz * rz / ry;
%% visualization
tempTheta = linspace(0,2*pi,1000);
tempX = mrx*cos(tempTheta);
tempY = mry*sin(tempTheta);
figure(99), hold on
plot(tempX,tempY,'b');

tempX = mrx_o*cos(tempTheta);
tempY = mry_o*sin(tempTheta);
plot(tempX,tempY,'k');
% tempX = linspace(-(mrx^2 - mry^2) / mrx, (mrx^2 - mry^2) / mrx, 1000);
% tempY = zeros(size(tempX));
% plot(tempX, tempY, 'b')
% axis([-33.0737   33.0737  -26.0856   26.0856])
% view(-30,80)
axis equal
axis vis3d
box on
%% new sampling
deltaTheta = 2*pi/nCrestPoints;
if ~mod(nRows,2)
    initialDisplacement = -deltaTheta/2;
else
    initialDisplacement = 0;
end
rowIndex = 1;
colIndex = 1;
medial_x = zeros(nRows, nCols);
medial_y = zeros(nRows, nCols);

for i = 0:nCrestPoints-1
    theta_ = pi - deltaTheta*floor(nRows/2) - initialDisplacement - deltaTheta*i;
    
    x_ = mrx*cos(theta_);
    y_ = mry*sin(theta_);
    plot(x_,y_,'b.','MarkerSize',30)
    
    xo_ = mrx_o*cos(theta_);
    yo_ = mry_o*sin(theta_);
    plot(xo_,yo_,'k.','MarkerSize',30)
    line([x_,xo_],[y_,yo_],'color','k')
    
    medial_x(rowIndex, colIndex) = x_;
    medial_y(rowIndex, colIndex) = y_;
    
    if i < nCols - 1
        colIndex = colIndex + 1;
    elseif i < nCols - 1 + nRows - 1
        rowIndex = rowIndex + 1;
    elseif i < nCols -1 + nRows -1 + nCols - 1
        colIndex = colIndex - 1;
    else
        rowIndex = rowIndex - 1;
    end
    if (i < nCols - 1 && i > 0) || (i > nCols + nRows - 2 && i < 2*nCols + nRows - 3)
        mx_ = (mrx^2 - mry^2) * cos(theta_) / mrx;
        my_= 0;
        %         line([mx_, x_], [my_, y_])
        %         plot(mx_, my_, 'r.')
        dx_ = (x_ - mx_);
        dy_ = (y_ - my_);
        line([x_, mx_],[y_,my_],'color','b')
        
        tau_ = linspace(0,1,floor(nRows/2)+1);
        if ~mod(nRows,2)
            tau_ = tau_(2:end-1);
        else
            tau_ = tau_(1:end-1);
        end
        if ~isempty(tau_)
            tempX_ = mx_ + tau_ * dx_;
            tempY_ = my_ + tau_ * dy_;
            plot(tempX_,tempY_,'r.','MarkerSize',30)
            
            if i < nCols - 1
                medial_x(2:numel(tempX_)+1,colIndex-1) = fliplr(tempX_);
                medial_y(2:numel(tempX_)+1,colIndex-1) = fliplr(tempY_);
            else
                medial_x(numel(tempX_)+1:end-1,colIndex+1) = (tempX_);
                medial_y(numel(tempX_)+1:end-1,colIndex+1) = (tempY_);
            end
            %             plot(tempX_, tempY_, 'r.')
        end
    end
    %     pause;
end
medial_z = zeros(nRows, nCols);
%%
up_x = zeros(nRows, nCols);
up_y = zeros(nRows, nCols);
up_z = zeros(nRows, nCols);

down_x = zeros(nRows, nCols);
down_y = zeros(nRows, nCols);
down_z = zeros(nRows, nCols);

crest_x = zeros(nRows, nCols);
crest_y = zeros(nRows, nCols);
crest_z = zeros(nRows, nCols);

medial_pdm_x = [];
medial_pdm_y = [];
medial_pdm_z = [];

boundary_pdm_x = [];
boundary_pdm_y = [];
boundary_pdm_z = [];
%%
figure, hold on
axis vis3d
for i = 1:nRows
    for j = 1:nCols
        mx = medial_x(i,j);
        my = medial_y(i,j);
        plot(mx,my,'r.')
        sB = my * mrx_o;
        cB = mx * mry_o;
        l = sqrt(sB*sB + cB*cB);
        if l ~= 0
            sB_n = sB / l;
            cB_n = cB / l;
        else
            sB_n = sB;
            cB_n = cB;
        end
        
        cA = l / (mrx_o * mry_o);
        if cA > 1
            warning('cA is greater than 1\n')
            keyboard
        end
        sA = sqrt(1 - cA*cA);
        
        sx = rx * cA * cB_n - mx;
        sy = ry * cA * sB_n - my;
        sz = rz * sA;
        
        l1 = sqrt(sx*sx + sy*sy + sz*sz);
        if l1 ~= 0
            sx_n = sx / l1;
            sy_n = sy / l1;
            sz_n = sz / l1;
        else
            sx_n = sx;
            sy_n = sy;
            sz_n = sz;
        end
        bx = (sx + mx);
        by = (sy + my);
        bz = (sz);
        if abs(bx*bx/rx/rx + by*by/ry/ry + bz*bz/rz/rz - 1) > 1e-15
            warning('The point is not on the surface of the ellipsoid\n');
        end
        
        medial_pdm_x = cat(1, medial_pdm_x, mx);
        medial_pdm_y = cat(1, medial_pdm_y, my);
        medial_pdm_z = cat(1, medial_pdm_z, 0);
        
        boundary_pdm_x = cat(1, boundary_pdm_x, bx);
        boundary_pdm_y = cat(1, boundary_pdm_y, by);
        boundary_pdm_z = cat(1, boundary_pdm_z, bz);
        boundary_pdm_x = cat(1, boundary_pdm_x, bx);
        boundary_pdm_y = cat(1, boundary_pdm_y, by);
        boundary_pdm_z = cat(1, boundary_pdm_z, -bz);
        plot3(mx, my, 0, 'r.')
        up_x(i,j) = bx;
        up_y(i,j) = by;
        up_z(i,j) = bz;
        line([mx, bx], [my, by], [0, bz],'Color','c','LineWidth',2)
        down_x(i,j) = bx;
        down_y(i,j) = by;
        down_z(i,j) = -bz;
        line([mx, bx], [my, by], [0, -bz],'Color','m','LineWidth',2)
        
        if i == 1 || i == nRows || j == 1 || j == nCols
            cx = rx * cB_n - mx;
            cy = ry * sB_n - my;
            cz = 0;
            
            v = [cx;cy;cz];
            v2 = [sx;sy;0];
            v3 = norm(v)*normc(v2);
            bx = (v3(1) + mx);
            by = (v3(2) + my);
            bz = v3(3);            
%             bx = (cx + mx);
%             by = (cy + my);
%             bz = 0;
            
            if abs(bx*bx/rx/rx + by*by/ry/ry + bz*bz/rz/rz - 1) > 1e-15
                warning('The point is not on the surface of the ellipsoid\n');
%                keyboard
            end
            line([mx, bx], [my, by], [0, bz],'Color','r','LineWidth',2)
            crest_x(i,j) = bx;
            crest_y(i,j) = by;
            crest_z(i,j) = bz;
            boundary_pdm_x = cat(1, boundary_pdm_x, bx);
            boundary_pdm_y = cat(1, boundary_pdm_y, by);
            boundary_pdm_z = cat(1, boundary_pdm_z, bz);
        end
    end
end
if flag
    filename = fullfile(dataDir,strcat(mesh_prefix,sprintf('_%d_%d_pre',nRows,nCols),'.m3d'));
    medialPoints = cat(3, medial_x, medial_y, medial_z);
    upPoints = cat(3, up_x, up_y, up_z);
    downPoints = cat(3, down_x, down_y, down_z);
    crestPoints = cat(3, crest_x, crest_y, crest_z);
    writeSREP( filename, nRows, nCols, medialPoints, upPoints, downPoints, crestPoints )
end

%%
pdm_x = cat(1, medial_pdm_x, boundary_pdm_x);
pdm_y = cat(1, medial_pdm_y, boundary_pdm_y);
pdm_z = cat(1, medial_pdm_z, boundary_pdm_z);
pdm_mat = cat(2, pdm_x, pdm_y, pdm_z);
pdm_secondMoment = pdm_mat'*pdm_mat;
[V2, ~] = eigs(pdm_secondMoment);
rot = V * V2';
if det(rot) < 0
    rot = -rot;
    warning('Determinant of rot is negative')
end
%%
up_pdm = cat(2, up_x(:), up_y(:), up_z(:));
down_pdm = cat(2, down_x(:), down_y(:), down_z(:));
crest_pdm = cat(2, crest_x(:), crest_y(:), crest_z(:));
medial_pdm = cat(2, medial_x(:), medial_y(:), medial_z(:));
% transform points
transformed_up_pdm = (rot * up_pdm')' + repmat(mesh_deformed_cog, size(up_pdm) ./ size(mesh_deformed_cog));
transformed_down_pdm = (rot * down_pdm')' + repmat(mesh_deformed_cog, size(down_pdm) ./ size(mesh_deformed_cog));
transformed_crest_pdm = (rot * crest_pdm')' + repmat(mesh_deformed_cog, size(crest_pdm) ./ size(mesh_deformed_cog));
transformed_medial_pdm = (rot * medial_pdm')' + repmat(mesh_deformed_cog, size(medial_pdm) ./ size(mesh_deformed_cog));
%%
transformed_medial_x = reshape( transformed_medial_pdm(:,1), nRows, nCols);
transformed_medial_y = reshape( transformed_medial_pdm(:,2), nRows, nCols);
transformed_medial_z = reshape( transformed_medial_pdm(:,3), nRows, nCols);

transformed_up_x = reshape(transformed_up_pdm(:,1), nRows, nCols);
transformed_up_y = reshape(transformed_up_pdm(:,2), nRows, nCols);
transformed_up_z = reshape(transformed_up_pdm(:,3), nRows, nCols);

transformed_down_x = reshape(transformed_down_pdm(:,1), nRows, nCols);
transformed_down_y = reshape(transformed_down_pdm(:,2), nRows, nCols);
transformed_down_z = reshape(transformed_down_pdm(:,3), nRows, nCols);

transformed_crest_x = reshape(transformed_crest_pdm(:,1), nRows, nCols);
transformed_crest_y = reshape(transformed_crest_pdm(:,2), nRows, nCols);
transformed_crest_z = reshape(transformed_crest_pdm(:,3), nRows, nCols);
%%
%filename = fullfile(dataDir,strcat(mesh_prefix,sprintf('_%d_%d',nRows,nCols),'.m3d'));
filename = fullfile(dataDir,strcat(mesh_prefix,'.m3d'));
medialPoints = cat(3, transformed_medial_x, transformed_medial_y, transformed_medial_z);
% upPoints = cat(3, transformed_up_x, transformed_up_y, transformed_up_z);
% downPoints = cat(3, transformed_down_x, transformed_down_y, transformed_down_z);
downPoints = cat(3, transformed_up_x, transformed_up_y, transformed_up_z);
upPoints = cat(3, transformed_down_x, transformed_down_y, transformed_down_z);
crestPoints = cat(3, transformed_crest_x, transformed_crest_y, transformed_crest_z);
writeSREP( filename, nRows, nCols, medialPoints, upPoints, downPoints, crestPoints )
end