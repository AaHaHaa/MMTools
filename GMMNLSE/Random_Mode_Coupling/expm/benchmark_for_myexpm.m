% Compare "myexpm_" with MATLAB "expm"

n = 5; % the size of a square matrix
l = 1:10:1000; % the number of square matrices

%reset(gpuDevice);

tmatlab = zeros(1,length(l));
tmyexpm_ = zeros(1,length(l));

i = 1;
for n3 = l
    m = rand(n,n,n3,'gpuArray');
    
    f1 = @()matexpm(m,n3);
    t1 = timeit(f1);
    tmatlab(i) = t1;
    
    f2 = @()myexpm_(m,[],[],1,0,1);
    t2 = gputimeit(f2);
    tmyexpm_(i) = t2;
    
    %f3 = pageexpm(m);
    %t3 = timeit(f3);
    %tpageexpm(i) = t3;
    
    i = i+1;
end

figure;
fp = get(gcf,'position');
set(gcf,'position',[fp(1:2) fp(3:4)*5/4]);
yyaxis right;
h1 = semilogy(l,tmatlab./tmyexpm_);
ylabel('speedup ratio');

yyaxis left;
h2 = semilogy(l,tmyexpm_,'b-');
hold on;
h3 = semilogy(l,tmatlab,'r-');
h4 = semilogy(l,tpageexpm,'r-');
hold off;
xlabel('the number of matrices');
ylabel('time (s)');
title(sprintf('"expm" model comparison\nfor size "%1ux%1u" randn matrices',n,n));
lhandle = legend('myexpm\_','MATLAB expm');
set(lhandle,'location','southeast');
set(h1,'linewidth',3); set(h2,'linewidth',3); set(h3,'linewidth',3);
set(gca,'fontsize',14);

print(gcf,'benchmark_expm', '-depsc');
save('benchmark_expm.mat');

    function matexpm(m,n3)
        for i = 1:n3
            expm(m(:,:,i));
        end
    end