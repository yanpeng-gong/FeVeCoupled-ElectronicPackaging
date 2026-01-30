function nodeNew = getMap(node)

sumNode = size(node,1);
nodeNew = zeros(sumNode,2);

fourNode = [0,0;48,44;48,60;0,44];

for n = 1:sumNode
    xi = node(n,1);
    eta = node(n,2);
    N = fun(xi,eta,0,2,4);
    nodeNew(n,:) = N'*fourNode;
end