function ff = getForce(nodeNew,sumNode,sumNodeNew,elem,auxT,k,press,direction)


ff = zeros(sumNodeNew*2,1);
% press = load('press.dat');
sumP = size(press,1);

nip = k+2;
[x,w] = gaussInt(nip);

for n = 1:sumP
    p = zeros(k+1,1);
    elemID = press(n,1);
    faceID = press(n,2);
    value = press(n,3);
    index = elem{elemID};
    indexEdge = auxT.elem2edge{elemID};
    Nv = length(index);
    v1 = 1:Nv; v2 = [2:Nv,1]; % loop index for vertices or edges
    
    if k == 1
        elem1 = [v1(:), v2(:)];
        faceNode = elem1(faceID,:);
        faceNodeID = index(faceNode);
        edgeNodeCoor = nodeNew(faceNodeID,:);
    elseif k == 2
        elem1 = [v1(:), v2(:), v1(:)+Nv];
        faceNode = elem1(faceID,:);
        midNode = sumNode+indexEdge(faceID);
        faceNodeID = [index([faceNode(1:2)]),midNode];
        edgeNodeCoor = nodeNew(faceNodeID,:);
    end
    
    L = nodeNew(faceNodeID(1:2),:);
    L = (L(1,:)-L(2,:));
    Normal = [L(2),-L(1)]/norm(L);
    for m = 1:nip
        N = fun(x(m),0,0,1,k+1);
        dN = dfunc(x(m),0,0,1,k+1);
        if k == 2
            N = [N(1);N(3);N(2)];
            dN = [dN(1);dN(3);dN(2)];
        end
        J = sqrt((dN'*edgeNodeCoor(:,1))^2+(dN'*edgeNodeCoor(:,2))^2);
        p = p+w(m)*N*value*J;
    end
    if strcmp(direction,'normal')
        ff(faceNodeID) = ff(faceNodeID)+p*Normal(1);
        ff(faceNodeID+sumNodeNew) = ff(faceNodeID+sumNodeNew)+p*Normal(2);
    elseif strcmp(direction,'x')
        ff(faceNodeID) = ff(faceNodeID)+p;
    elseif strcmp(direction,'y')
        ff(faceNodeID+sumNodeNew) = ff(faceNodeID+sumNodeNew)+p;
    end
end
ff = sparse(ff);