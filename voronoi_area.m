clear all
clc 
clf
x = 200.*rand(500,1);
y = 200.*rand(500,1);
h = voronoi(x,y);
set(h, 'LineWidth', 4, 'MarkerSize', 10)
dt = delaunayTriangulation(x(:),y(:));    %creating column vectors of x and y and then performing delaunay triangulation. Delaunay triangulation and voronoi diagrams are geometric duals - you can get one from the other. [V,R] = voronoiDiagram(dt) %Here V gives the coordinates of the vertices and R gives the voronoi regions of each point. It contains the indices of the vertices present in each voronoi region.
[V,R] = voronoiDiagram(dt)    %V represents the coordinates of the vertices in the voronoi diagram. R represents the voronoi regions of each point. Each row contains the indices of the vertices that are present in that particular voronoi region. 
for i = 1:size(R,1);          %for loop that will run for the length of R (number of rows in R, hence 1)
    temp = R{i};        %stores that particular row of R. This value changes with every iteration where the row number changes. 
    X = zeros(size(temp,2),1);        %create column zero vectors of X and Y of size equal to the number of elements in that particular row of R
    Y = zeros(size(temp,2),1);
    for j = 1:size(temp,2);        %for loop to go through the different elements present in that particular row. The size of j depends on the length of temp which indicate the number of elements in that particular row. 
        temp1 = temp(j);
        temp2 = V(temp1,:);
        X(j) = V(temp1,1);       %X coordinate of that particular vertex stored as jth element
        Y(j) = V(temp1,2);       %Y coordinate of that particular vertex stored as jth element
        j = j+1;           %subsequent iterations to go to the next element in that row
    end 
    X             %display the X coordinates of all the vertices of a particular voronoi region
    Y             %display the Y coordinates of all the vertices of a partciular voronoi region
   Area(i) = polyarea(X,Y)      %display the area of that particular voronoi region and store it as ith element.
   i = i+1; %subsequent iterations to go to the next row
end 
Area = Area(:)   %displays all the areas of the voronoi regions. 
meanArea = mean(Area, "omitnan") %mean area is calculated, excluding the NaN values otherwise meanarea is coming NaN
