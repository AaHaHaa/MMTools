function resized_spatial_field = spatial_field_resizing( spatial_field, ratio )
%SPATIAL_FIELD_RESIZING Resize the spatial field, magnifying or shrinking
%
% spatial_field - (Nx,Nx,num_modes) matrix
% ratio - the resizing ratio. magnification: ratio>1; shrinking: ratio<1

Nx = size(spatial_field,1);
num_modes = size(spatial_field,3);

resized_spatial_field = spatial_field;
for i = 1:num_modes
    if ratio>1
        Nx_crop = floor(Nx/ratio);
        x_crop = floor((Nx-Nx_crop)/2):floor((Nx-Nx_crop)/2)+Nx_crop-1;
        crop_spatial_field = spatial_field(x_crop,x_crop,i);
        x_new = linspace(1,Nx_crop,Nx);
        [xmesh,ymesh] = meshgrid(x_new,x_new);
        resized_spatial_field(:,:,i) = interp2(crop_spatial_field,xmesh,ymesh,'spline');
    elseif ratio<1
        Nx_enlarge = Nx/ratio;
        x_new = linspace(-floor((Nx_enlarge-Nx)/2),-floor((Nx_enlarge-Nx)/2)+Nx_enlarge,Nx);
        [xmesh,ymesh] = meshgrid(x_new,x_new);
        resized_spatial_field(:,:,i) = interp2(spatial_field(:,:,i),xmesh,ymesh,'spline',0);
    end
end

end

