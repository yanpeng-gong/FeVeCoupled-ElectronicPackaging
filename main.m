%% FE-VE Coupled Method
%  Author: Yanpeng Gong/Sishuai Li, Beijing University of Technology
%  Contact: yanpenggong@gmail.com (Feel free to contact for any questions)
%  Reference: Gong et al. (2026), Eng. Anal. Bound. Elem., 184:106640
clear;
load("matlab.mat");
figure;
showmesh(node,elem);

sumNode = size(node,1);
nodeT  = find(node(:,1)>19.999);
nodeX = find(node(:,1)<0.01);
nodeY = find(node(:,2)<0.01);
sX = size(nodeX,1);
sY = size(nodeY,1);
fixNode3 = [nodeX,ones(sX,1),zeros(sX,1)];
fixMes3 = [fixNode3(:,1),fixNode3(:,3)];
fixNode4 = [nodeY,2*ones(sY,1),zeros(sY,1)];
fixMes4 = [fixNode4(:,1)+(fixNode4(:,2)-1)*sumNode,fixNode4(:,3)];

auxT = auxstructure(node,elem);
face = findFace(node,elem,nodeT);
press = [face,ones(length(face),1)*(5)];
nodeForce = getForce(node,sumNode,sumNode,elem,auxT,1,press,'x');

E = 10000; nu = 0.3; 
mat.E = E; mat.nu = nu;

GK = globalK(node,elem,elem1,elem2,mat);
[GK,F] = boundaryCondition(GK,nodeForce,fixMes3(:,1),fixMes3(:,2),fixMes4(:,1),fixMes4(:,2),1);
uh = GK\F;
uh = full(uh);
ux = uh(1:sumNode);
uy = uh(sumNode+1:end);

[stress,mises] = calculateStress(node,elem,elem1,elem2,uh,mat);
figure;
showsolution(node+[ux,uy],elem,mises);
