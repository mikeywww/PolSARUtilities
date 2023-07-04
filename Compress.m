% 将全极化的协方差矩阵压缩成双极化或简缩极化的协方差矩阵
function out = Compress(C, mode, options)
arguments
    % 全极化协方差矩阵
    C
    % 要压缩成的形式
    mode string
    % 是否转为斯托克斯矢量（仅限压缩模式属于简缩极化）
    options.asStokesVector = false
end
mode = lower(mode);
[height,width,~,~] = size(C);
switch mode
    case 'pi/4'
        out = zeros(height, width, 2, 2);
        out(:,:,1,1) = ...
            (C(:,:,1,1)+C(:,:,2,2)/2+sqrt(2)*real(C(:,:,1,2)))/2;
        out(:,:,1,2) = ...
            (C(:,:,1,3)+C(:,:,2,2)/2+(C(:,:,1,2)+C(:,:,2,3)/sqrt(2)))/2;
        out(:,:,2,1) = conj(out(:,:,1,2));
        out(:,:,2,2) = ...
            (C(:,:,3,3)+C(:,:,2,2)/2+sqrt(2)*real(C(:,:,2,3)))/2;
        if options.asStokesVector
            out = to_stokes_vector(out);
        end
    case 'ctlr'
        out = zeros(height, width, 2, 2);
        out(:,:,1,1) = ...
            (C(:,:,1,1)+C(:,:,2,2)/2-sqrt(2)*imag(C(:,:,1,2)))/2;
        out(:,:,1,2) = ...
            (1i*C(:,:,1,3)-1i*C(:,:,2,2)/2+(C(:,:,1,2)+C(:,:,2,3))/...
            sqrt(2))/2;
        out(:,:,2,1) = conj(out(:,:,1,2));
        out(:,:,2,2) = ...
            (C(:,:,3,3)+C(:,:,2,2)/2-sqrt(2)*imag(C(:,:,2,3)))/2;
        if options.asStokesVector
            out = to_stokes_vector(out);
        end
    case 'hhhv'
        out = zeros(height, width, 2, 2);
        out(:,:,1,1) = C(:,:,1,1);
        out(:,:,1,2) = C(:,:,1,2)/sqrt(2);
        out(:,:,2,1) = conj(out(:,:,1,2));
        out(:,:,2,2) = C(:,:,2,2)/sqrt(2);
    case 'hhvv'
        out = zeros(height, width, 2, 2);
        out(:,:,1,1) = C(:,:,1,1);
        out(:,:,1,2) = C(:,:,1,3);
        out(:,:,2,1) = conj(out(:,:,1,2));
        out(:,:,2,2) = C(:,:,3,3);
    otherwise
        error("不支持的操作");
end
end

function g = to_stokes_vector(C2)
[height,width,~,~] = size(C2);
g = zeros(height, width, 4);
g(:,:,1) = C2(:,:,1,1)+C2(:,:,2,2);
g(:,:,2) = C2(:,:,1,1)-C2(:,:,2,2);
g(:,:,3) = 2*real(C2(:,:,1,2));
g(:,:,4) = -2*imag(C2(:,:,1,2));
end