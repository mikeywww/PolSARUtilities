% 对极化协方差/相干矩阵进行极化分解
function decomp = Decompose(method, datatype, input)
arguments
    % 分解方法
    method string
    % 数据类型
    datatype string
    % 输入数据
    input
end
datatype = lower(datatype);
method = lower(method);
decompose_func = str2func(strcat(datatype, '_', method));
decomp = decompose_func(input);
end

function out = ctlr_freeman(in)
G = real(in);
g0 = G(:,:,1);
g1 = G(:,:,2);
g2 = G(:,:,3);
g3 = G(:,:,4);
x1 = g0-sqrt(g1.^2+g2.^2+g3.^2);
x = x1*0.65;

mask1 = g3<0;
pd1 = ((g0+g3-x).*(g0-g3-x)-g1.^2-g2.^2)./(2*(g0-g3+x));
ps1 = ((g0-g3-x).^2+g1.^2+g2.^2)./(2*(g0-g3-x));
mask2 = ~mask1;
pd2 = ((g0+g3-x).^2+g1.^2+g2.^2)./(2*(g0+g3-x));
ps2 = ((g0-g3-x).*(g0+g3-x)-g1.^2-g2.^2)./(2*(g0+g3-x));

Pv = x;
Pd = pd1.*mask1+pd2.*mask2;
Ps = ps1.*mask1+ps2.*mask2;
out = struct(...
    'Ps', Ps, ...
    'Pd', Pd, ...
    'Pv', Pv ...
    );
end

function out = quad_freeman(in)
Fv = squeeze(4*real(in(:,:,2,2)));
A = squeeze(in(:,:,1,3)-Fv/8);
B = squeeze(real(in(:,:,3,3))-3*Fv/8);
C = squeeze(real(in(:,:,1,1))-3*Fv/8);
D = conj(A);

% real(C_Lee(m,n,1,3)-Fv/8)>0
mask1 = real(squeeze(in(:,:,1,3))-Fv/8)>0;
a1 = -1;
Fd1 = real((B.*C-A.*D)./(A+B+C+D));
Fs1 = B-Fd1;
b1 = (A+Fd1)./Fs1;

% real(C_Lee(m,n,1,3)-Fv/8)<=0
mask2 = real(squeeze(in(:,:,1,3))-Fv/8)<=0;
b2 = 1;
Fs2 = real((A.*D-B.*C)./(A+D-B-C));
Fd2 = B-Fs2;
a2 = (A-Fs2)./Fd2;

Fs = Fs1.*mask1+Fs2.*mask2;
Fd = Fd1.*mask1+Fd2.*mask2;
alpha = a1.*mask1+a2.*mask2;
beta = b1.*mask1+b2.*mask2;

Ps = abs(Fs.*(1+abs(beta).^2));
Pd = abs(Fd.*(1+abs(alpha).^2));
Pv = abs(Fv);
out = struct(...
    'Ps', Ps, ...
    'Pd', Pd, ...
    'Pv', Pv ...
    );
end

function out = quad_freeman_advance(in)
[height,width,~,~] = size(in);
B = (in(:,:,2,2)-in(:,:,3,3))/2;
E = real(in(:,:,2,3)+in(:,:,3,2))/2;
cos4t = B./sqrt(B.^2+E.^2);
sin4t = E./sqrt(B.^2+E.^2);
cos2t = sqrt((1+cos4t)/2);
sin2t = sin4t./(2*cos2t);

Ps = zeros(height, width);
Pd = zeros(height, width);
Pv = zeros(height, width);
for m=1:height
    for n=1:width
        q = [1 0 0;0 cos2t(m,n) sin2t(m,n);0 -sin2t(m,n) cos2t(m,n)];
        t = squeeze(in(m,n,:,:));
        t1 = q*t*q';
        if t1(1,1)<=t1(3,3)
            Pv(m,n) = real(3*t1(1,1));
            Ps(m,n) = 0;
            Pd(m,n) = real(t1(2,2)+t1(3,3)-2*t1(1,1));
        else
            Pv(m,n) = 3*real(t1(3,3));
            x11 = real(t1(1,1)-t1(3,3));
            x22 = real(t1(2,2)-t1(3,3));
            if abs(t1(1,2))^2>x11*x22
                if x11>x22
                    Ps(m,n) = x11+x22;
                    Pd(m,n) = 0;
                else
                    Ps(m,n) = 0;
                    Pd(m,n) = x11+x22;
                end
            else
                if x11>x22
                    Ps(m,n) = x11+abs(t1(1,2))^2/x11;
                    Pd(m,n) = x22-abs(t1(1,2))^2/x11;
                else
                    Ps(m,n) = x11-abs(t1(1,2))^2/x22;
                    Pd(m,n) = x22+abs(t1(1,2))^2/x22;
                end
            end
        end
    end
end

out = struct(...
    'Ps', Ps, ...
    'Pd', Pd, ...
    'Pv', Pv ...
    );
end

function out = quad_cloude(in)
[height,width,~,~] = size(in);
H = zeros(height, width);
A = zeros(height, width);
a = zeros(height, width);
for m = 1:height
    for n = 1:width
        t = reshape(in(m,n,:,:), [3,3]);
        [u,l] = eig(t);
        l = sort(sum(abs(l), 1), 'descend');
        p = l/sum(l);
        H(m,n) = -sum(p.*log(p)/log(3));
        a(m,n) = 180*sum(p.*acos(abs(u(1,:))))/pi;
        A(m,n) = (l(2)-l(3))/(l(2)+l(3));
    end
end
out = struct(...
    'H', H, ...
    'A', A, ...
    'a', a ...
    );
end

function out = dual_cloude(in)
[height,width,~,~] = size(in);
H = zeros(height, width);
a = zeros(height, width);
A = zeros(height, width);
for m = 1:height
    for n = 1:width
        t = squeeze(in(m,n,:,:));
        [u,l] = eig(t);
        l = sum(abs(l), 1);
        p = l/sum(l);
        H(m,n) = -sum(p.*log(p)/log(2));
        a(m,n) = 180*sum(p.*acos(abs(u(1,:))))/pi;
        A(m,n) = abs(l(1)-l(2))/sum(l);
    end
end
out = struct(...
    'H', H, ...
    'a', a, ...
    'A', A ...
    );
end
