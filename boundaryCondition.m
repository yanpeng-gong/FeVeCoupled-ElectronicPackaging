function[KNew,FNew] = boundaryCondition(K,F,fixDof1,gievDisp1,fixDof2,gievDisp2,m)


KNew = K;
FNew = F;

if m == 1
    FNew(fixDof1) = gievDisp1;
    FNew(fixDof2) = gievDisp2;

    KNew(fixDof1,:) = 0;
    KNew1 = KNew+sparse(fixDof1,fixDof1,ones(size(fixDof1,1),1),size(K,1),size(K,1));
    KNew1(fixDof2,:) = 0;
    KNew2 = KNew1+sparse(fixDof2,fixDof2,ones(size(fixDof2,1),1),size(K,1),size(K,1));
    KNew =KNew2;
else % 乘大数法
    alpha=max(K(:))*1e5;
    for n = 1:length(fixDof)
        dof = fixDof(n);
        KNew(dof,dof) = alpha*K(dof,dof);
        FNew(dof) = alpha*K(dof,dof)*gievDisp(n);
    end
end
