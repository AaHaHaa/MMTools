function func = blkdiag_rand(use_gpu)
%BLKDIAG_RAND It generates block-diagonal random matrices, each with 
%a certain type of special matrix. It can generate unitary,
%close-to-identity unitary, hermitian, and skew-hermitian random matrices.
%
%   func.haar - totally unitary random matrix
%   func.identity_exp - close-to-identity unitary random matrix
%   func.identity_rootn - close-to-identity unitary random matrix
%   func.hermitian - hermitian random matrix
%   func.skew_hermitian - skew-hermitian random matrix
%
%   func.gpu_yes - whether to use "gpu" or not
%
%   For more information, please look into
%   "random_unitary_matrix_generator.m" or
%   "random_special_matrix_generator.m"

% Unitary matrix
fun1 = random_unitary_matrix_generator(use_gpu);
% Hermitian, Skew-Hermitian
fun2 = random_special_matrix_generator(use_gpu);

% status: whether to use gpu or not
status = struct('gpu_yes',use_gpu);

func = catstruct(fun1,fun2,status);

end

