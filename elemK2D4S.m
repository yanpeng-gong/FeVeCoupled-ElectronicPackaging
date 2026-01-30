function [stress_gp, N_gp] = elemK2D4S(EX,mu,x1,y1,x2,y2,x3,y3,x4,y4,reduce,ux,uy,index2)
% 计算单元刚度矩阵
% 四节点四边形单元
% 采用二点高斯积分，全积分高斯点位置+0.577350269189626,-0.577350269189626，权系数为1
% 若采取减缩积分，积分点位置为0，权系数为2
% 可以采用更高阶的高斯积分
% mnode == 4
% 4----3
% |    |
% |    |
% 1----2


D = EX/(1-mu^2)*[1,mu,0;mu,1,0;0,0,(1-mu)/2];

 
if reduce == 0 % 采用完全积分
    nip = 2; % 单维方向的积分点数
    x = [-0.577350269189626,0.577350269189626]; % 积分点
    w = [1,1]; % 权系数
else
    nip = 1;
    x = 0; % 积分点
    w = 2; % 权系数
end
stress_gp = zeros(4, 3); % 存储每个积分点的应力分量[σx, σy, τxy]
N_gp = zeros(4, 4);
gp_count = 0;
for m = 1:nip
    for n = 1:nip
        gp_count = gp_count + 1;
        [B,N] = elemB2D4P(x1,y1,x2,y2,x3,y3,x4,y4,x(n),x(m));
        N_gp(gp_count, :) = N;
        elemU = [ux(index2);uy(index2)];
        eS00 = D * B * elemU;
        
          
        stress_gp(gp_count, :) = eS00;

    end 
end

