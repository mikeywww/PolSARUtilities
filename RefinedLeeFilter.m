% 精致Lee滤波的极化降斑
function C_Lee = RefinedLeeFilter(C, nLooks)
arguments
    % 极化SAR数据
    C
    % 多视视数
    nLooks
end
[height,width,d,~] = size(C);
SPAN = zeros(height, width);
for i = 1:d
    SPAN = SPAN+squeeze(abs(C(:,:,i,i)));
end
SPAN_tmp = [zeros(3, width+3*2); % 在SPAN周围补一圈0
    zeros(height, 3), SPAN, zeros(height, 3);
    zeros(3, width+3*2)];

W1 = [-1,0,1;-1,0,1;-1,0,1];
W2 = [0,1,1;-1,0,1;-1,-1,0];
W3 = [1,1,1;0,0,0;-1,-1,-1];
W4 = [1,1,0;1,0,-1;0,-1,-1];
W = cat(3, W1, W2, W3, W4);   % 边缘检测算子

pw1 = repmat([0,0,0,1,1,1,1], [7,1]);
pw2 = [1,1,1,1,1,1,1;
    0,1,1,1,1,1,1;
    0,0,1,1,1,1,1;
    0,0,0,1,1,1,1;
    0,0,0,0,1,1,1;
    0,0,0,0,0,1,1;
    0,0,0,0,0,0,1];
pw3 = repmat([1,1,1,1,0,0,0]', [1,7]);
pw4 = flip(pw2, 2);
pw5 = flip(pw1, 2);
pw6 = flip(pw4, 1);
pw7 = flip(pw3, 1);
pw8 = flip(pw2, 1);
PW = {pw1, pw2, pw3, pw4, pw5, pw6, pw7, pw8}; % prewitt算子

C_Lee = zeros(height, width, d, d);
for m = (1:height)+3
    for n = (1:width)+3
        window = SPAN_tmp((m-3):(m+3),(n-3):(n+3));
        M = clac_m3x3(window);
        i = calc_template(M, W);
        pw = PW{i};
        b = calc_b(window, pw, 1/nLooks);
        c_mean = clac_mean(C, pw, m-3, n-3, height, width, d);
        c_caret = c_mean+b*(C(m-3,n-3,:,:)-c_mean);
        C_Lee(m-3,n-3,:,:) = c_caret;
    end
end

end


% 计算3x3的平均矩阵
function M3x3 = clac_m3x3(window7x7)
M3x3 = zeros(3, 3);
for i = 1:3
    for j = 1:3
        M3x3(i,j) = mean( ...
            window7x7((1:3)+2*(i-1),(1:3)+2*(j-1)), ...
            "all");
    end
end
end


% 计算匹配的模板
function I = calc_template(M3x3, W)
m = M3x3.*W;
[~, I] = max(abs(sum(m, [1,2])));
switch I
    case 1
        delta1 = abs(M3x3(2,1)-M3x3(2,2));
        delta2 = abs(M3x3(2,3)-M3x3(2,2));
        if delta1 > delta2
            I = I+4;
        end
    case 2
        delta1 = abs(M3x3(3,1)-M3x3(2,2));
        delta2 = abs(M3x3(1,3)-M3x3(2,2));
        if delta1 > delta2
            I = I+4;
        end
    case 3
        delta1 = abs(M3x3(3,2)-M3x3(2,2));
        delta2 = abs(M3x3(1,2)-M3x3(2,2));
        if delta1 > delta2
            I = I+4;
        end
    case 4
        delta1 = abs(M3x3(3,3)-M3x3(2,2));
        delta2 = abs(M3x3(1,1)-M3x3(2,2));
        if delta1 > delta2
            I = I+4;
        end
end
end


% 估计b参数
function b = calc_b(window7x7, pw, var_v)
z = window7x7.*pw;
z_mean = mean(z, 'all');
var_z = mean((z-z_mean).^2, "all");
var_x = (var_z - (z_mean^2)*var_v)/(1+var_v);
b = var_x/var_z;
end


% 计算协方差矩阵的统计平均值
function C_mean = clac_mean(C, pw, m, n, height, width, d)
rg = max([m-3, 1]):min([m+3, height]);
if m < 4
    rg_pw = 5-m:7;
elseif m > height - 4
    rg_pw = 1:4+(height-m);
else
    rg_pw = 1:7;
end
az = max([n-3, 1]):min([n+3, width]);
if n < 4
    az_pw = 5-n:7;
elseif n > width - 4
    az_pw = 1:4+(width-n);
else
    az_pw = 1:7;
end

C_masked = zeros(length(rg), length(az), d, d);
for i = 1:d
    for j = 1:d
        C_masked(:,:,i,j) = C(rg,az,i,j).*pw(rg_pw,az_pw);
    end
end
C_mean = mean(C_masked, [1,2]);
end