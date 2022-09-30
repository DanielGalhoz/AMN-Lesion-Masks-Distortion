step = 0.0005;

all_files_st = dir('Masks');
all_files = {all_files_st.name};
all_files([all_files_st.isdir]) = [];

for file = string(all_files)
    
    [~,name,ext] = fileparts(file);
    eye = extractBetween(name,'_','_MASK');
    if length(eye)~=1
        eye = extractBetween(name,'_','_mask');
    end
    while any(strfind(eye,'_'))
        eye = extractAfter(eye,'_');
    end
    if length(eye)~=1
        error('Error parsing filename: %s.',name)
    elseif ~strcmp(eye,'OD') && ~strcmp(eye,'OS')
        error('Error parsing filename: %s.',name)
    end

    mask = imread(file);
    if size(mask,1) ~= size(mask,2)
        if size(mask,1) == 1536 || size(mask,2) == 1536
            %fprintf('N rows != N columns\n')
            mask = imresize(mask,[1536 1536]);
        else
            error('N rows != N columns')
        end
    end
    im_size = size(mask,1);
    if im_size == 1536
        res = 0.0053086597472429276;
    elseif im_size == 768
        res = 0.01110921148210764;
    else
        error('Unregistered image size')
    end
    if size(mask,3) == 3
        mask = rgb2gray(mask);
    end
    mask_bw = imbinarize(mask);

    %figure
    %imshow(mask_bw);

    masks = {};
    mask_rest = mask_bw;
    center_x = NaN;
    center_y = NaN;
    i = 1;
    while any(any(mask_rest))
        mask_border = bwareafilt(mask_rest,1,'largest');
        
        mask_fill = imfill(mask_border,'holes');
        if any(mask_fill-mask_border,'all')
            masks{end+1} = mask_border;
        else
            stats = regionprops(mask_fill);
            centroid = stats.Centroid;
            center_x = round(centroid(1));
            center_y = round(centroid(2));
        end 
        
        mask_rest = imbinarize(mask_rest - mask_border);
        i = i+1; 
        % figure
        % imshow(mask_rest);
    end
    
    if isnan(center_x)
            center_x = round(im_size/2);
            center_y = round(im_size/2);
    end

    disp = calc_disp_map(0,eye);
    orig_masks = cell(size(masks));
    dist_masks = cell(size(masks));
    for m = 1:length(masks)
        mask_border = masks{m};
        mask_NaN_border = NaN(im_size,im_size);
        mask_area = imfill(mask_border,'holes');
        mask_NaN_area = NaN(im_size,im_size);
        for i = 1:size(mask_border,1)
            for j = 1:size(mask_border,2)
                if mask_border(i,j) == 1
                    mask_NaN_border(i,j) = 1;
                end
            end
        end
        for i = 1:size(mask_area,1)
            for j = 1:size(mask_area,2)
                if mask_area(i,j) == 1
                    mask_NaN_area(i,j) = 1;
                end
            end
        end


        border_cord_x = NaN(im_size,im_size);
        border_cord_y = NaN(im_size,im_size);
        for i = 1:size(mask_NaN_border,1)
            for j = 1:size(mask_NaN_border,2)
                if ~isnan(mask_NaN_border(i,j))
                    border_cord_x(i,j) = (j-center_x-1)*res;
                    border_cord_y(i,j) = -(i-center_y-1)*res;
                end
            end
        end

%         img_disp = NaN(im_size,im_size);
%         img_cord_x = NaN(im_size,im_size);
%         img_cord_y = NaN(im_size,im_size);
%         for i = 1:size(mask_NaN_area,1)
%             for j = 1:size(mask_NaN_area,2)
%                 if ~isnan(mask_NaN_area(i,j))
%                     img_cord_x(i,j) = (j-center_x-1)*res;
%                     img_cord_y(i,j) = -(i-center_y-1)*res;
%                     [theta,rho] = cart2pol(img_cord_x(i,j),img_cord_y(i,j));
% 
%                     if theta > pi
%                         theta = -theta + 2*pi;
%                     end
% 
%                     if rho > disp.rho(end)
%                         %print('Point outside range')
%                         img_disp(i,j) = 0;
%                     else
%                         [~,idx_rho]=min(abs(disp.rho-rho));
%                         [~,idx_theta]=min(abs(disp.theta-theta));
% 
%                         img_disp(i,j) = disp.map(idx_rho,idx_theta);
%                     end
%                 end
%             end
%         end
% 
%         % figure
%         % surf(img_cord_x,img_cord_y,img_disp,'edgecolor','none')

        border_cord_1D = [];
        border_cord_1D(1,:) = border_cord_x(~isnan(border_cord_x));
        border_cord_1D(2,:) = border_cord_y(~isnan(border_cord_y));
        dist_cord_1D = zeros(size(border_cord_1D));

        for i = 1:size(border_cord_1D,2)

            [theta,rho] = cart2pol(border_cord_1D(1,i),border_cord_1D(2,i));

            RHO = 0:step:rho;
            [X,Y] = pol2cart(repmat(theta,size(RHO)),RHO);

            Z = zeros(1,length(RHO));
            for j = 1:length(Z)
                Z(1,j) = calc_z(theta,RHO(j),eye);
            end

            total_distance = zeros(length(RHO)-1,1);
            for j = 1:(length(RHO)-1)
                distance(1) = X(j+1) - X(j);
                distance(2) = Y(j+1) - Y(j);
                distance(3) = Z(j+1) - Z(j);
                total_distance(j) = norm(distance);
            end
            new_rho = sum(total_distance);
            [dist_cord_1D(1,i),dist_cord_1D(2,i)] = pol2cart(theta,new_rho);
        end

        % figure
        % scatter(border_cord_1D(1,:),border_cord_1D(2,:),4,'r')
        % hold on
        % scatter(dist_cord_1D(1,:),dist_cord_1D(2,:),4,'b')

        dist_mask_border = uint8(zeros(size(mask_border)));
        for i = 1:size(dist_cord_1D,2)
            img_x = round(dist_cord_1D(1,i)/res)+center_x+1;
            img_y = -round(dist_cord_1D(2,i)/res)+center_y+1;
            dist_mask_border(img_y,img_x) = 255;
        end

        % figure 
        % imshow(mask_border)
        % figure 
        % imshow(dist_mask_border)

        dist_mask_border_filter = imdilate(dist_mask_border,strel('disk',1));
        dist_mask_border_filter = imfilter(dist_mask_border_filter,fspecial('average'));
        % figure
        % imshow(dist_mask_border_filter)

        % Add centroid and fovea center
        mask_fill = imfill(dist_mask_border_filter,'holes');
        mask_fill = imbinarize(mask_fill);
        stats = regionprops(mask_fill);
        centroid = stats.Centroid;
        dist_mask_border_filter(round(centroid(2)),round(centroid(1))) = 255;
        dist_mask_border_filter(center_y,center_x) = 255;
   
        % Add centroid
        mask_fill = imfill(mask_border,'holes');
        stats = regionprops(mask_fill);
        centroid = stats.Centroid;
        mask_border(round(centroid(2)),round(centroid(1))) = 255;
        
        orig_masks{m} = mask_border;
        dist_masks{m} = dist_mask_border_filter;
       
    end

    comb_orig_mask = sum(cat(3,orig_masks{:}),3);
    comb_dist_mask = sum(cat(3,dist_masks{:}),3);

    % figure
    % imshow(comb_dist_mask)

    C = imfuse(mask_bw,comb_dist_mask,'ColorChannels',[2 0 1]);
    figure
    imshow(C)

    imwrite(comb_orig_mask, append(name,'_orig.tif'));
    imwrite(comb_dist_mask, append(name,'_disp.tif'));
end
%% Functions
function disp = calc_z(theta,rho,eye)

    %Coefficients of piecewise cubic splines fit to pooled total displacements
    x_nasal = [0.0 0.6243 2.6231 3.9632 5.0]; 
    a_nasal = [-4.3774 1.2022 0.0 0.0];
    b_nasal = [1.1856  -1.5470  0.0 0.0];
    c_nasal = [0.6898 0.5770 -0.1098 0.0];
    d_nasal = [0.0 0.4841 0.1470 0.0];
    
    x_temporal = [0.0 1.2337 2.5360 5.0]; %The x range can be extended to 5
    a_temporal = [-0.1030 1.3537 0.0];
    b_temporal = [-0.7650 -0.8921 0.0];
    c_temporal = [0.9336 -0.0885 -0.0689];
    d_temporal = [0.0 0.5374 0.1639];
    
    theta = abs(theta);

    if rho >= x_nasal(end)
        %print('Point outside range')
        disp = 0;
    else
        for j = 1:(length(x_nasal)-1)
            if rho >= x_nasal(j) && rho < x_nasal(j+1)
                i = j;
            end
        end
        T = rho - x_nasal(i);
        disp_nasal = a_nasal(i)/6*T^3 + b_nasal(i)/2*T^2 + c_nasal(i)*T + d_nasal(i);

        for j = 1:(length(x_temporal)-1)
            if rho >= x_temporal(j) && rho < x_temporal(j+1)
                i = j;
            end
        end
        T = rho - x_temporal(i);
        disp_temporal = a_temporal(i)/6*T^3 + b_temporal(i)/2*T^2 + c_temporal(i)*T + d_temporal(i);
        
        if strcmp(eye,'OD')
            disp = interp1([0 pi],[disp_nasal disp_temporal],theta,'linear');
        elseif strcmp(eye,'OS')
            disp = interp1([0 pi],[disp_temporal disp_nasal],theta,'linear');
        end
        if disp < 0
            disp = 0;
        end 
    end
end

function disp = calc_disp_map(print,eye)

    %Coefficients of piecewise cubic splines fit to pooled total displacements
    x_nasal = [0.0 0.6243 2.6231 3.9632 5.0]; 
    a_nasal = [-4.3774 1.2022 0.0 0.0];
    b_nasal = [1.1856  -1.5470  0.0 0.0];
    c_nasal = [0.6898 0.5770 -0.1098 0.0];
    d_nasal = [0.0 0.4841 0.1470 0.0];
    
    x_temporal = [0.0 1.2337 2.5360 5.0]; %The x range can be extended to 5
    a_temporal = [-0.1030 1.3537 0.0];
    b_temporal = [-0.7650 -0.8921 0.0];
    c_temporal = [0.9336 -0.0885 -0.0689];
    d_temporal = [0.0 0.5374 0.1639];
    
    disp.rho = linspace(0,5,50);
    disp.theta = linspace(0,pi,50);
    disp.map = zeros(length(disp.rho),length(disp.theta));
    for m = 1:length(disp.rho)
        for j = 1:(length(x_nasal)-1)
            if disp.rho(m) >= x_nasal(j) && disp.rho(m) < x_nasal(j+1)
                i = j;
            end
        end
        T = disp.rho(m) - x_nasal(i);
        
        if strcmp(eye,'OD')
            disp.map(m,1) = a_nasal(i)/6*T^3 + b_nasal(i)/2*T^2 + c_nasal(i)*T + d_nasal(i);
        elseif strcmp(eye,'OS')
            disp.map(m,end) = a_nasal(i)/6*T^3 + b_nasal(i)/2*T^2 + c_nasal(i)*T + d_nasal(i);
        end
    end
    for m = 1:length(disp.rho)
        for j = 1:(length(x_temporal)-1)
            if disp.rho(m) >= x_temporal(j) && disp.rho(m) < x_temporal(j+1)
                i = j;
            end
        end
        T = disp.rho(m) - x_temporal(i);

        if strcmp(eye,'OD')
            disp.map(m,end) = a_temporal(i)/6*T^3 + b_temporal(i)/2*T^2 + c_temporal(i)*T + d_temporal(i);
        elseif strcmp(eye,'OS')
            disp.map(m,1) = a_temporal(i)/6*T^3 + b_temporal(i)/2*T^2 + c_temporal(i)*T + d_temporal(i);
        end
    end
    for m = 1:length(disp.rho)
        disp.map(m,:) = interp1([disp.theta(1) disp.theta(end)],[disp.map(m,1) disp.map(m,end)],disp.theta,'linear'); 
    end

    if print == true
        [Theta,Rho] = meshgrid(disp.theta,disp.rho);
        figure
        surf(Rho,Theta,disp.map);
        surf(Rho,Theta,disp.map,'edgecolor','none');
        [X,Y] = pol2cart(Theta,Rho);
        figure
        hold on
        surf(X,Y,disp.map) 
        surf(X,-Y,disp.map)
        surf(X,Y,disp.map,'edgecolor','none');
        surf(X,-Y,disp.map,'edgecolor','none');
    end
end
