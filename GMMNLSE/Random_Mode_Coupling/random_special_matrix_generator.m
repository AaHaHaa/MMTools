function fun = random_special_matrix_generator(use_gpu)
%RANDOM_SPECIAL_MATRIX_GENERATOR It generates block-diagonal Hermitian and skew-Hermtian matrices
%   Use:
%       (1) First generate the function handle,
%
%           fun = random_special_matrix_generator;
%           fun = random_special_matrix_generator(use_gpu);
%           
%           use_gpu - 0 or 1 (true or false); use gpuArray or not
%           If no "use_gpu", it'll run with CPU.
%
%       (2) Then call the function and generate the output.
%           
%           random_matrix = fun.hermitian(n,n3);
%           random_matrix = fun.skew_hermitian(n,n3);
%           
%           n - an array specifying the dimension (n,n) of the block-diagonal random matrices,
%               where each block is a unitary matrix
%           n3 - the number of random matrices to generate
%           random_matrix (output) - (n,n,n3) matrix with each nxn a random matrix
%

if ~exist('use_gpu','var')
    use_gpu = false;
end

fun.hermitian = @(n,n3) rand_hermitian(n,n3,use_gpu);
fun.skew_hermitian = @(n,n3) rand_skew_hermitian(n,n3,use_gpu);

end

%% RAND_HERMITIAN
function rand_matrix = rand_hermitian(n,n3,use_gpu)
%RAND_HERMITIAN It generates a random hermitian matrix.

% "n" should be in a row, not a column.
if size(n,1) > 1
    n = n';
end

if n == 1
    rand_matrix = matrix_initialization(n,n3,use_gpu,@randn);
else
    
    % Initialization
    rand_matrix = matrix_initialization(n,n3,use_gpu,@zeros);
    
    sn = sum(n);
    
    % Point out the linear indices of independent varibles within upper-half
    % Hermitian.
    previous_n = cumsum(n);
    previous_n = [0 previous_n(1:end-1)];
    num_independent_elements = [0 (n-1).*n/2];
    previous_num_independent_elements = cumsum(num_independent_elements);
    idx = zeros(1,sum(num_independent_elements),'uint32');
    for i = 1:length(n)
        nidx = n(i);
        for m = 1:(nidx-1)
            idx(previous_num_independent_elements(i) + (m-1)*m/2 + (1:m)) = (previous_n(i)+m)*sn + previous_n(i) + (1:m);
        end
    end
    n3idx_increment = uint32(permute(((1:n3)-1)*sn^2,[1 3 2]));
    idx = bsxfun(@plus, idx,n3idx_increment);
    
    % Put random numbers in these slots of independent variables.
    rand_matrix(idx(:)) = (randn(1,sum(num_independent_elements)*n3) + 1i*randn(1,sum(num_independent_elements)*n3))/sqrt(2);
    
    % Calculate the lower-half part according to the upper-half part.
    if use_gpu
        rand_matrix = rand_matrix + pagefun(@ctranspose, rand_matrix);
    else
        rand_matrix = rand_matrix + conj(permute(rand_matrix,[2 1 3]));
    end
    
    % Diagonal part: real numbers
    idx = bsxfun(@plus, uint32( (1:sn) + ((1:sn)-1)*sn ), n3idx_increment);
    rand_matrix(idx(:)) = randn(1,sn*n3);
end

end

%% RAND_SKEW_HERMITIAN
function rand_matrix = rand_skew_hermitian(n,n3,use_gpu)
%RAND_HERMITIAN It generates a random skew-hermitian matrix.

% "n" should be in a row, not a column.
if size(n,1) > 1
    n = n';
end

if n == 1
    rand_matrix = 1i*matrix_initialization(n,n3,use_gpu,@randn);
else
    
    % Initialization
    rand_matrix = matrix_initialization(n,n3,use_gpu,@zeros);
    
    sn = sum(n);
    
    % Point out the linear indices of independent varibles within upper-half
    % Hermitian.
    previous_n = cumsum(n);
    previous_n = [0 previous_n(1:end-1)];
    num_independent_elements = [0 (n-1).*n/2];
    previous_num_independent_elements = cumsum(num_independent_elements);
    idx = zeros(1,sum(num_independent_elements),'uint32');
    for i = 1:length(n)
        nidx = n(i);
        for m = 1:(nidx-1)
            idx(previous_num_independent_elements(i) + (m-1)*m/2 + (1:m)) = (previous_n(i)+m)*sn + previous_n(i) + (1:m);
        end
    end
    n3idx_increment = uint32(permute(((1:n3)-1)*sum(n)^2,[1 3 2]));
    idx = bsxfun(@plus, idx,n3idx_increment);
    
    % Put random numbers in these slots of independent variables.
    rand_matrix(idx(:)) = (randn(1,sum(num_independent_elements)*n3) + 1i*randn(1,sum(num_independent_elements)*n3))/sqrt(2);
    
    % Calculate the lower-half part according to the upper-half part.
    if use_gpu
        rand_matrix = rand_matrix - pagefun(@ctranspose, rand_matrix);
    else
        rand_matrix = rand_matrix - conj(permute(rand_matrix,[2 1 3]));
    end
    
    % Diagonal part: imaginary numbers
    idx = bsxfun(@plus, uint32( (1:sn) + ((1:sn)-1)*sn ), n3idx_increment);
    rand_matrix(idx(:)) = 1i*randn(1,sn*n3);
end

end

%% MATRIX_INITIALIZATION
function output = matrix_initialization(n,n3,use_gpu,func)
%MATRIX_INITIALIZATION It initializes (sum(n),sum(n),n3) matrix according 
%to "use_gpu" and the functions defined in varargin. It also works with 
%block-diagonal initialization for "randn".

single_matrix = cell(1,length(n));
switch func2str(func)
    case 'randn'
        for nidx = 1:length(n)
            if use_gpu
                single_matrix{nidx} = feval(func,n(nidx),n(nidx),n3,'gpuArray');
            else
                single_matrix{nidx} = feval(func,n(nidx),n(nidx),n3);
            end
        end

        if length(n) == 1
            output = single_matrix{1};
        else
            output = matrix_initialization(sum(n),1,use_gpu,@zeros);
            for n3idx = 1:n3
                single_matrix_in_n3idx = cellfun(@(m) m(:,:,n3idx), single_matrix,'UniformOutput',false);
                output(:,:,n3idx) = blkdiag(single_matrix_in_n3idx{:});
            end
        end
    case {'zeros','ones'}
        if use_gpu
            output = feval(func,sum(n),sum(n),n3,'gpuArray');
        else
            output = feval(func,sum(n),sum(n),n3);
        end
    otherwise
        error('Not supported yet.');
end

end