% Voronoi tesselation for cell culture simulation
% Anwesha Barua 10_09_24

function [M, seedpoint, vertexcounter, cellarea, cellperimeter, numNeighbour, cvn] = voronoi_tesselation(x,y)

% clc
% clear all
% clf

% Generate random points for x and y
% s = 300;
% M0 = 100;
% x = s*rand(M0,1);
% y = s*rand(M0,1);
% x = [119.1964; 78.7692;  23.5916; 66.0359; 145.4955;  154.5889; 251.7028;  123.9456; 120.6527;   285.8061];
% y = [100.5295; 40.1251; 208.6309; 276.8320; 146.2574;  151.7077; 1.3877;  56.5406; 233.8918; 183.2478];
xmin = min(x); xmax = max(x); ymin = min(y); ymax = max(y);
M0 = size(x,1); % no. of points
AreaTotal = (xmax-xmin)*(ymax-ymin);
AreaAve = AreaTotal/M0;
AreaLimit = 4;

% Plot Voronoi diagram using the randomly generated points
h = voronoi(x,y);
hold on
axis equal
set(h, 'linewidth', 1, 'markersize', 3);

% Truncate the display area with a dashed line box
axis([xmin xmax ymin ymax]);
plot([xmin xmax xmax xmin xmin],[ymin ymin ymax ymax ymin], 'k--','linewidth',2);

% Delaunay triangulation - geometric dual of Voronoi diagram. 
% V contains the coordinates of all the vertices - the first row always
% represents an infinite vertex.
% Each row of R contains the indices of the vertices making up that
% particular voronoi region. 
dt = delaunayTriangulation(x(:),y(:));
[V, R] = voronoiDiagram(dt);
k =0;
% Plot all the vertices - red are the ones lying outside the display box
% while white represents the vertices inside the box
for i = 1:size(V,1)
    if((V(i,1) > xmin && V(i,1) <xmax) && V(i,2) > ymin && V(i,2) < ymax)
        plot(V(i,1),V(i,2), 'w*')
    else
        plot(V(i,1),V(i,2),'r*')
    end
end

% Get the array containing the coordinates of the vertices making up each
% Voronoi region
% for loop that will run for the length of R (number of rows in R, hence 1)
for i = 1:size(R,1)
    temp = R{i};                % stores that particular row of R. This value changes with every iteration where the row number changes.
    X = zeros(size(temp,2),1);  % create column zero vectors of X and Y of size equal to the number of elements in that particular row of R
    Y = zeros(size(temp,2),1);
    % for loop to go through each vertex index in each row of R
    for j = 1:size(temp,2)
        temp1 = temp(j);
        X(j) = V(temp1,1); % X and Y coordinates of that particular vertex stored as jth element
        Y(j) = V(temp1,2);
    end
    Area(i) = polyarea(X,Y);
    flag = 0;
    if ~isnan(Area(i))
        if Area(i) < AreaLimit*AreaAve
            flag = 1;
            validcell = isnan(Area(i));
            Vcell = [X Y];
            polyin = polyshape(Vcell);
            k = k+1;
            % Find the number of vertices associated with each voronoi cell
            vertexcounter(k) = size(X,1);
            % Find the area and perimeter of the voronoi cells
            cellarea(k,:) = area(polyin);
            cellperimeter(k,:) = perimeter(polyin);
            plot(polyin);
            % Find which seed point lies inside kth cell
            in = inpolygon(x,y,X,Y);
            seedpoint(k) = find(in==1);
        end
    end
end

meanArea = mean(Area, "omitnan");

% Number of Voronoi cells with finite area less than 3*average area
M = size(seedpoint,2);

% Find neighbouring cells
n = size(x,1);
CellN = zeros(n,n);
for i = 1:n
    for j = i+1:n
        s = size(intersect(R{i}, R{j})); %intersect returns the data common to both
        if (s(2) > 1)
            CellN(i,j) = 1;
            CellN(j,i) = 1;
        end
    end
    if ~isnan(Area(i)) && Area(i) < AreaLimit*AreaAve
        numNeighbour(i,:) = sum(CellN(i,:)); %number of neighbours
    end
end
numNeighbour(numNeighbour == 0) = [];

% Find coordinates of the common edges with the neighbours for cells with
% closed polygons
cvn = zeros(M0,5);
for i = 1:M
    cn = seedpoint(i);
    nN(cn,1) = sum(CellN(cn,:)); % Number of neighbours
    cv = zeros(M0,5);
    for j = 1:M0
        if CellN(cn,j) == 1
            cvtemp = intersect(R{cn},R{j});
            cv(j,1) = cn;
            cv(j,2) = j;
            cv(j,3) = cvtemp(1,1);
            cv(j,4) = cvtemp(1,2);
            vx1 = V(cvtemp(1,1),1);
            vx2 = V(cvtemp(1,2),1);
            vy1 = V(cvtemp(1,1),2);
            vy2 = V(cvtemp(1,2),2);
            cv(j,5) = sqrt((vx2-vx1)^2+(vy2-vy1)^2); % Find length of common edge for every pair of neighbouring cells
        end
    end
    cvn(:,:,i) = cv;
   
end
nN(nN==0) = [];

% Consolidated table containing seedpoint number, x and y
% coordinates, number of vertices associated with the voronoi cell, area,
% perimeter and number of neighbours
Cellproperties = zeros(M,7);
colNames = {'seedpoint', 'x coordinate of seedpoint', 'y coordinate of seedpoint', 'number of vertices per cell', 'cell area', 'cell perimeter', 'number of neighbours'};
for p = 1:M
    Cellproperties(p,:) = [seedpoint(p), x(p), y(p), vertexcounter(p), cellarea(p), cellperimeter(p), numNeighbour(p)];
end
output = array2table(Cellproperties, 'VariableNames', colNames);


                  

