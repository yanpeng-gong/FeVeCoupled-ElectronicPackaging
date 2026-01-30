function GK = globalK(node,elem,elem1,elem2,mat)

sumElem = size(elem,1); % the number of element
sumNode = size(node,1);
sumElem1 = size(elem1,1);
sumElem2 =size(elem2,1);

elemLen1 = cellfun('length',elem1); 
nnz1 = sum((2*elemLen1).^2);
elemLen2 = cellfun('length',elem2); 
nnz2 = sum((2*elemLen2).^2);
ii = zeros(nnz1,1); jj = zeros(nnz1,1); ss = zeros(nnz1,1); 
ii2 = zeros(nnz2,1); jj2 = zeros(nnz2,1); ss2 = zeros(nnz2,1); 

ia = 0;
for n = 1:sumElem1
    index1 = elem1{n};
    coor1 = node(index1,:);
    AK = elemKVEM(mat.E,mat.nu,coor1);
    AB = reshape(AK',1,[]);

    % --------- assembly index for ellptic projection -----------
    indexDof = [index1, index1+sumNode];
    Ndof = length(indexDof);
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ss(ia+1:ia+Ndof^2) = AB(:);
    ia = ia + Ndof^2;
end
GK1 = sparse(ii,jj,ss,sumNode*2,sumNode*2);

h=1;
    ia = 0;
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
        [K] = elemK2D4(mat.E,mat.nu,h,x1,y1,x2,y2,x3,y3,x4,y4);  % 生成单元刚度阵
        AB = reshape(K',1,[]);
        indexDof = [index2, index2+sumNode];
        Ndof = length(indexDof);
        ii2(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
        jj2(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
        ss2(ia+1:ia+Ndof^2) = AB(:);
        ia = ia + Ndof^2;
    end
GK2 = sparse(ii2,jj2,ss2,sumNode*2,sumNode*2);

GK = GK1+GK2;


