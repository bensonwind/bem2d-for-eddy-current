function [edge] = findboundary(elem)
    nv = size(elem,2);
    allEdge = uint32(sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2));

    ne = nv; % number of edges in one element
    Neall = length(allEdge);
    [~, i2, j] = myunique(allEdge);
    NT = size(elem,1);
    i1(j(Neall:-1:1)) = Neall:-1:1; 
    i1 = i1';
    bdEdgeidx = i1(i1==i2);

    edge = allEdge(bdEdgeidx,:);
end

