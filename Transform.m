% 极化协方差矩阵C和极化相干矩阵T的转换
function result = Transform(how, datatype, input)
arguments
    % 转换方法
    how string
    % 数据类型
    datatype string
    % 待转换数据
    input
end
[height,width,d,~] = size(input);
switch datatype
    case 'hhhv'
        U = [sqrt(2), 0; 0, 1+1i]/sqrt(2);
    case 'vvvh'
        U = [sqrt(2), 0; 0, 1+1i]/sqrt(2);
    case 'hhvv'
        U = [1, 1; 1, -1];
    case 'quad'
        U = [1,0,1;1,0,-1;0,sqrt(2),0]/sqrt(2);
end
result = zeros(height, width, d, d);
switch how
    case 'C2T'
        for i=1:height
            for j=1:width
                c = squeeze(input(i,j,:,:));
                result(i,j,:,:) = U*c*U';
            end
        end
    case 'T2C'
        for i=1:height
            for j=1:width
                t = squeeze(input(i,j,:,:));
                result(i,j,:,:) = U'*t*U;
            end
        end
end
end