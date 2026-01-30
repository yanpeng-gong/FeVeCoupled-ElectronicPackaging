function [stress,mises] = calculateStress(node,elem,elem1,elem2,uh,mat)
ndim = 2;

sumElem = size(elem,1); % the number of element
sumNode = size(node,1);
sumElem1 = size(elem1,1);
sumElem2 =size(elem2,1);


ux = uh(1:sumNode);
uy = uh(sumNode+1:end);

stress1 = zeros(sumNode,3);

nodeUsed1 = zeros(sumNode,1);

for n = 1:sumElem1
    index1 = elem1{n};

    coor1 = node(index1,:);
    [dM,Pis] = calculatePi(coor1);

    A1 = [1,0;0,0;0,1];
    A2 = [0,0;0,1;1,0];
    A = [A1,A2];

    Pis = blkdiag(Pis,Pis);

    E = mat.E; nu =mat.nu;
    % D = E*(1-nu)/((1+nu)*(1-2*nu))*[1,nu/(1-nu),0;nu/(1-nu),1,0;0,0,(1-2*nu)/(2*(1-nu))]; % plane strain
    C = E/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];

    dM1 = blkdiag(dM',dM');
    eS = C*A*dM1*Pis*[ux(index1);uy(index1)];
    stress1(index1,:) = stress1(index1,:)+eS';
    nodeUsed1(index1) = nodeUsed1(index1)+1;
end

if ndim == 2
    nip = 4;
else
    nip = 8;
end
ndimB = ndim*3-3;  % 单元B矩阵的维数，二维为3，三维为6
stress2 = zeros(sumNode,ndimB);
nodeUsed2 = zeros(sumNode,1);
for n = 1:sumElem2
    index2 = elem2{n};
    coor2 = node(index2,:);
    x1 = coor2(1,1);
    y1 = coor2(1,2);
    x2 = coor2(2,1);
    y2 = coor2(2,2);
    x3 = coor2(3,1);
    y3 = coor2(3,2);
    x4 = coor2(4,1);
    y4 = coor2(4,2);
    [stress_gp, N_gp] = elemK2D4S(mat.E,mat.nu,x1,y1,x2,y2,x3,y3,x4,y4,0,ux,uy,index2); 
    STRESS00 =  N_gp'*stress_gp;
    nodeUsed2(index2) = nodeUsed2(index2)+1;
    stress2(index2,:) = stress2(index2,:)+STRESS00;
end

stress=stress1+stress2;
nodeUsed=nodeUsed1+nodeUsed2;
stress = stress./nodeUsed;

smax = (stress(:,1)+stress(:,2))./2+sqrt(((stress(:,1)-stress(:,2))./2).^2+stress(:,3).^2);
smin = (stress(:,1)+stress(:,2))./2-sqrt(((stress(:,1)-stress(:,2))./2).^2+stress(:,3).^2);


mises = sqrt((smax.^2+smin.^2+(smax-smin).^2)./2);
