% Wishart迭代分类算法
function cls = WishartClassification(method, C, decomp, options)
arguments
    % 分类算法
    method string
    % 协方差矩阵
    C
    % 极化分解结果
    decomp
    options.MAX_ITER = inf
end

[height,width,~,~] = size(C);

switch method
    case 'H/a'
        [init_classes,num_classes] = ...
            make_init_cls_H_a(decomp.H, decomp.a);
    case 'H/a/A'
        [init_classes,num_classes] = ...
            make_init_cls_H_a_A(decomp.H, decomp.a, decomp.A);
    case 'Freeman'
        [init_classes,num_classes] = ...
            make_init_cls_Freeman(decomp.Ps, decomp.Pd, decomp.Pv, C);
    case 'Dual-pol H/a'
        [init_classes,num_classes] = ...
            make_init_cls_Dual_H_a(decomp.H, decomp.a);
    case 'CTLR Freeman'
        [init_classes,num_classes] = ...
            make_init_cls_Freeman(decomp.Ps, decomp.Pd, decomp.Pv, C);
    otherwise
        error('不支持的操作');
end

transfer_rate = 1;
iter_input = init_classes;
i = 1;
while transfer_rate > 0.1 && i <= options.MAX_ITER
    Vs = calc_centers(num_classes, C, iter_input);
    
    iter_output = zeros(height, width);
    for i = 1:height
        for j = 1:width
            c = squeeze(C(i,j,:,:));
            distances = calc_distances(num_classes, c, Vs);
            [~,new_class] = min(distances);
            iter_output(i,j) = new_class;
        end
    end
    
    transfer_rate = sum(iter_output~=iter_input, 'all')/(height*width);
    
    iter_input = iter_output;
    i = i+1;
end
cls = iter_output;

end


function [init_classes, num_classes] = make_init_cls_H_a(H, a)
[height,width] = size(H);
num_classes = 8;
init_classes = zeros(height, width);
for i = 1:height
    for j = 1:width
        if H(i,j) <= 0.5
            if a(i,j) >= 47.5
                init_classes(i,j) = 1;
            elseif a(i,j) >= 42.5
                init_classes(i,j) = 2;
            else
                init_classes(i,j) = 3;
            end
        elseif H(i,j) <= 0.9
            if a(i,j) >= 50
                init_classes(i,j) = 4;
            elseif a(i,j) >= 40
                init_classes(i,j) = 5;
            else
                init_classes(i,j) = 6;
            end
        else
            if a(i,j) >= 55
                init_classes(i,j) = 7;
            elseif a(i,j) >= 40
                init_classes(i,j) = 8;
            else
                init_classes(i,j) = 9; % 不存在
            end
        end
    end
end
end


function [init_classes, num_classes] = make_init_cls_H_a_A(H, a, A)
[init_classes,~] = make_init_cls_H_a(H, a);
num_classes = 16;
init_classes = init_classes+(8*(A>0.5));
end


function [init_classes, num_classes] = make_init_cls_Dual_H_a(H, a)
[height,width] = size(H);
num_classes = 8;
init_classes = zeros(height, width);
for i = 1:height
    for j = 1:width
        if H(i,j) <= 0.6
            if a(i,j) >= 46
                init_classes(i,j) = 1;
            elseif a(i,j) >= 40
                init_classes(i,j) = 2;
            else
                init_classes(i,j) = 3;
            end
        elseif H(i,j) <= 0.95
            if a(i,j) >= 46
                init_classes(i,j) = 4;
            elseif a(i,j) >= 34
                init_classes(i,j) = 5;
            else
                init_classes(i,j) = 6;
            end
        else
            if a(i,j) >= 46
                init_classes(i,j) = 7;
            elseif a(i,j) >= 33.2
                init_classes(i,j) = 8;
            else
                init_classes(i,j) = 9; % 不存在
            end
        end
    end
end
end


function [init_classes, num_classes] = make_init_cls_Freeman(Ps, Pd, Pv, C)
[height,width] = size(Ps);
% 初始化30个类别
[~,scatter_classes] = max(abs(cat(3, Ps, Pd, Pv)), [], 3);
surface = scatter_classes==1;
double = scatter_classes==2;
volume = scatter_classes==3;
s_borders = divide_vector(sort(Ps(surface)), 10);
d_borders = divide_vector(sort(Pd(double)), 10);
v_borders = divide_vector(sort(Pv(volume)), 10);
init_classes = zeros(height, width);
for m=1:height
    for n=1:width
        switch scatter_classes(m,n)
            case 1
                init_classes(m,n) = 10;
                for idx=1:9
                    if Ps(m,n)<=s_borders(idx)
                        init_classes(m,n) = idx;
                        break;
                    end
                end
            case 2
                init_classes(m,n) = 20;
                for idx=1:9
                    if Pd(m,n)<=d_borders(idx)
                        init_classes(m,n) = 10+idx;
                        break;
                    end
                end
            case 3
                init_classes(m,n) = 30;
                for idx=1:9
                    if Pv(m,n)<=v_borders(idx)
                        init_classes(m,n) = 20+idx;
                        break;
                    end
                end
        end
    end
end
% 合并类别
num_classes = 8;
merged_classes = init_classes;
s_classes = 1:10;
d_classes = 11:20;
v_classes = 21:30;
for cnt=30:-1:(num_classes+1)
    c = calc_centers(cnt, C, merged_classes);
    distances = Inf(cnt);
    for i=1:(cnt-1)
        for j=(i+1):cnt
            d = (log(real(det(c{i})))+log(real(det(c{j})))+...
                trace(c{i}\c{j}+c{j}\c{i}))/2;
            distances(i,j) = d;
        end
    end
    % 确保合并的类别属于一个散射体制，同时总数不会过多
    while 1
        [~,I] = min(distances, [], "all");
        [p,q] = ind2sub(size(distances), I); % 计算出类别p和类别q
        % p和q的散射机制相同
        is_surface = ismember(p, s_classes) && ismember(q, s_classes);
        is_double = ismember(p, d_classes) && ismember(q, d_classes);
        is_volume = ismember(p, v_classes) && ismember(q, v_classes);
        if is_surface || is_double || is_volume
            % p和q合并后数量不过多
            if sum(merged_classes==p, "all")+...
                    sum(merged_classes==q, "all")<=2*height*width/(cnt-1)
                break;
            end
        end
        distances(p,q) = inf;
    end
    % 合并类别
    for m=1:height
        for n=1:width
            if merged_classes(m,n)==q
                merged_classes(m,n) = p;
            elseif merged_classes(m,n)>q
                merged_classes(m,n) = merged_classes(m,n)-1;
            end
        end
    end
    % 更新三种散射机制的类别表
    if is_surface
        s_classes(s_classes>=q) = s_classes(s_classes>=q)-1;
        s_classes = sort(unique(s_classes));
        d_classes = d_classes-1;
        v_classes = v_classes-1;
    elseif is_double
        d_classes(d_classes>=q) = d_classes(d_classes>=q)-1;
        d_classes = sort(unique(d_classes));
        v_classes = v_classes-1;
    elseif is_volume
        v_classes(v_classes>=q) = v_classes(v_classes>=q)-1;
        v_classes = sort(unique(v_classes));
    end
end
init_classes = merged_classes;
end


% 计算类别中心
function Vs = calc_centers(numClasses, C, classes)
[~,~,d,~] = size(C);
Vs = cell(1, numClasses);
for m = 1:numClasses
    Vm = zeros(d, d);
    mask = classes==m;
    count = sum(mask, "all");
    if count ~= 0
        C_masked = C.*mask;
        Vm = squeeze(sum(C_masked, [1 2])/count);
    end
    Vs{m} = Vm;
end
end


% 计算类别间的距离
function ds = calc_distances(numClasses, c, Vs)
ds = zeros(1, numClasses);
for m = 1:numClasses
    Vm = Vs{m};
    if all(Vm==0, "all")
        ds(m) = inf;
    else
        ds(m) = log(real(det(Vm)))+real(trace(Vm\c));
    end
end
end


% 在排好序的向量中等间隔地取N个值
function borders = divide_vector(V, N)
borders = zeros(1,N);
len = length(V);
for i=1:N
    borders(i) = V(round(len*i/N));
end
end