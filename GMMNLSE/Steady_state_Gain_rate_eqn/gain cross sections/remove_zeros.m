function input_matrix = remove_zeros( input_matrix )
%REMOVE_ZEROS clear zeros in "input_matrix" and interpolate them with the points nearby
%
%   input_matrix: a matrix (not multidimensional one)

sM = size(input_matrix);
numM = prod(sM(2:end)); % the number of column vectors

% "interp" doesn't support gpuArray.
if isequal(class(input_matrix),'gpuArray')
    gather_field = true;
    input_matrix = gather(input_matrix);
else
    gather_field = false;
end

% Real or Complex number
if isreal(input_matrix(1))
    is_real = true;
else
    is_real = false;
end

% Start removing zeros
for n = 1:numM
    max_number = max(abs(input_matrix(:,n)));
    if max_number == 0 % skip it if all zeros
        continue;
    else
        iszero = input_matrix(:,n)==0;
        zero_idx = find(iszero);
        nonzero_idx = find(~iszero);
        if length(nonzero_idx) > 1
            input_matrix(zero_idx,n) = interp1(nonzero_idx,input_matrix(~iszero,n),zero_idx,'pchip');
        end
        
        if is_real
            input_matrix(input_matrix(:,n)==0,n) = eps(max_number); % if it's zero after the approximation, then put "eps" instead.
        else % complex
            input_matrix(input_matrix(:,n)==0,n) = complex(eps(max_number)); % if it's zero after the approximation, then put "eps" instead.
        end
    end
end

% Send it back to GPU if it's gpuArray before.
if gather_field
    input_matrix = gpuArray(input_matrix);
end

end

