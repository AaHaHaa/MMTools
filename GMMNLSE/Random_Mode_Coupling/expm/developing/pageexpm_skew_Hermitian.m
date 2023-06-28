function expD = pageexpm_skew_Hermitian(D)
%UNTITLED2 Paged matrix exponential for skew-Hermitian matrices

if isequal(class(D),'gpuArray')
    use_gpu = true;
else
    use_gpu = false;
end

[T,Q] = pageTriDiagHouseholder(D);

[V,e] = pageTriDiagEig(T,Q);

if use_gpu
    expD = pagefun(@mrdivide, V.*exp(e),V);
else % CPU
    MATLAB_version = version('-release');
    MATLAB_year = str2double(MATLAB_version(1:4));
    if MATLAB_year < 2023
        expD = complex(zeros(size(D)));
        for i = 1:size(D,3)
            expD(:,:,i) = V(:,:,i).*exp(e(:,:,i))/V(:,:,i);
        end
    else
        expD = pagemrdivide(V.*exp(e),V);
    end
end

end