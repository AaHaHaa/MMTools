function check_nummodes(sim, fiber, fields)
%CHECK_NUMMODES It checks the consistency of the number of modes of betas, SR tensors, and fields.

% Get the number of spatial modes from "SR".
num_spatial_modes = size(fiber.SR,1);

% -------------------------------------------------------------------------
% Check the number of spatial modes from betas and SR.
if sim.scalar
    num_spatial_modes_betas = size(fiber.betas,2);
else % polarized fields
    num_modes_betas = size(fiber.betas,2);
    if num_modes_betas == num_spatial_modes
        num_spatial_modes_betas = num_modes_betas;
    else % betas has already included polarization modes
        num_spatial_modes_betas = num_modes_betas/2;
    end
end

if num_spatial_modes_betas ~= num_spatial_modes
    if size(fiber.betas,1)==1 % if the user give betas in row vector for single mode. It should be in column vector.
        error('GMMNLSE_propagate:NumModesError',...
            '"betas" of each mode should be be in the form of "column vector".\nThe number of spatial modes of betas and SR tensors should be the same.');
    else
        error('GMMNLSE_propagate:NumModesError',...
            'The number of spatial modes of betas and SR tensors should be the same.');
    end
end

% -------------------------------------------------------------------------
% Check the number of modes of fields
num_modes_fields = size(fields,2); % spatial modes (+ polarization modes)

field_modes_mismatched = false;
if sim.scalar
    if num_modes_fields ~= num_spatial_modes
        field_modes_mismatched = true;
    end
else
    if num_modes_fields ~= 2*num_spatial_modes
        field_modes_mismatched = true;
    end
end
if field_modes_mismatched
    error('GMMNLSE_propagate:NumModesError',...
        'The number of modes of fields doesn''t match those of betas and SR tensors.\nIf not scalar fields, num_modes(field)=2*num_spatial_modes(SR,betas).');
end

end