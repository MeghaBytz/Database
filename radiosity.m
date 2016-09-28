%function vfactorinfo = radiosity()
clear all
close all
%make these global variables
r=1; % reactor radius %actually .08
L=0.075; % reactor length
l_p=L/2; % half length
area=2*pi*r*(r+L); % total surface area


lambda = .51;
%parametrize reactor into polygon circle
% circle
C=[0 0];
%Define center
cx = 0;
cy = 0;

theta=0:2*pi/360:2*pi; % the angle
center(:,1) = cx*ones(length(theta),1);
center(:,2) = cy*ones(length(theta),1);
m=r*[cos(theta')+C(1) sin(theta')+C(2)]+center; 
L = ones(length(m),1)*lambda;
sheath = [m(:,1) L m(:,2)]; %switch y and z axes
%later change input to function to scaled geometry
R1 = [3,4,-1,1,1,-1,0.5,0.5,-0.75,-0.75]';
R2 = [3,4,-0.2,0.2,.2,-.2,0.5, 0.5, 0,0]';
% xq = [-1,0,.3];
% yq = [-.7,.2,-.4];
xv = [-1,-1,-.2,-.2,.2,.2,1,1];
yv = [-.7,.5,.5,0,0,.5,.5,-.7];
%[in,on] = inpolygon(xq,yq,xv,yv);
%define surface element as element in which 2 nodes are on edge
gm = [R1,R2];
sf = 'R1-R2';
ns = char('R1','R2');
ns = ns';
g = decsg(gm,sf,ns);
figure;
pdegplot(g);
model = createpde;
geometryFromEdges(model,g);
mesh = generateMesh(model);
pdeplot(model);
pg = pdeGeometryFromEdges(g);
test = zeros(length(mesh.Nodes),1);
%compute viewfactors of all mesh triangles from coordinates above
vfactorinfo = zeros(length(mesh.Elements),4);
figure;
scatter3(sheath(:,1),sheath(:,2),sheath(:,3))
hold on
scatter3(mesh.Nodes(1,:),mesh.Nodes(2,:),test','filled')
%detect surface mesh elements
edgeNodes = [];
edgeElems = [];
edgeNodesCoordinates = [];
nodeViewFactor=[]

for i = 1:length(mesh.Nodes)
    xq = mesh.Nodes(1,i);
    yq = mesh.Nodes(2,i);
    [in,on] = inpolygon(xq,yq,xv,yv);
    if on == 1
        edgeNodes = [edgeNodes i];
        edgeNodesCoordinates = [edgeNodesCoordinates;xq yq]
    end
end

for i = 1:length(mesh.Elements)
    elem = mesh.Elements(:,i);
    intersect(elem,edgeNodes)
        if length(intersect(elem,edgeNodes)>=2)
            edgeElems = [edgeElems i];
        end
end
%vfactor info = [vfactor12,vfactor21,area1,area2]

for i=1:length(edgeElems)%get coordinates of nodes attached to each element
    elem = mesh.Elements(:,edgeElems(i));
    nodeCoord1 = [mesh.Nodes(:,elem(1))' 0];
    nodeCoord2 = [mesh.Nodes(:,elem(2))' 0];
    nodeCoord3 = [mesh.Nodes(:,elem(3))' 0];
   [vfactorinfo(edgeElems(i),1),vfactorinfo(edgeElems(i),2),vfactorinfo(edgeElems(i),3),vfactorinfo(edgeElems(i),4)]=viewfactor([nodeCoord1;nodeCoord2;nodeCoord3],[sheath],12);
   if ismember(nodeCoord1(1:2),edgeNodesCoordinates,'rows')
    nodeViewFactor = [nodeViewFactor; mesh.Nodes(:,elem(1))' vfactorinfo(edgeElems(1))];
   end
   if ismember(nodeCoord2(1:2),edgeNodesCoordinates,'rows')
    nodeViewFactor = [nodeViewFactor; mesh.Nodes(:,elem(2))' vfactorinfo(edgeElems(2))];
   end
   if ismember(nodeCoord2(1:2),edgeNodesCoordinates,'rows')
    nodeViewFactor = [nodeViewFactor; mesh.Nodes(:,elem(3))' vfactorinfo(edgeElems(3))];
   end
end
nodeViewFactorSorted = sortrows(nodeViewFactor,3);
figure;
pdeplot(model)
hold on
a=50;
views = unique(nodeViewFactorSorted(:,3));
colorScheme = linspace(1,10,length(views));
xColorTable = [];
for i=1:length(nodeViewFactor(:,3))
   [u,ia,ib] = intersect(nodeViewFactorSorted(i,3),views)
    xColorTable = [xColorTable colorScheme(ib)]
end
c = xColorTable;

scatter(nodeViewFactorSorted(:,1),nodeViewFactorSorted(:,2),a,c,'filled')
