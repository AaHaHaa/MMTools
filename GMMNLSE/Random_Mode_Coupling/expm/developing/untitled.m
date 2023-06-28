clearvars;

addpath('expm');

m=rand(20,20,2000)+1i*rand(20,20,2000); m = 1i*(m+pagectranspose(m))/2;
disp(ishermitian(m(:,:,1),'skew'));

n = gather(m);

for j = 1:20
    
tic;
for i = 1:2000
    expm(n(:,:,i));
end
t1(j)=toc;
a=expm(m(:,:,100)); %disp(a);

tic;
b=pageexpm_skew_Hermitian(m);
t2(j)=toc;
%disp(b(:,:,100));

tic;
c=myexpm_(m,[],[],1,0,0);
t3(j)=toc;

end

disp(mean(t1)/mean(t2))
disp(mean(t1)/mean(t3))