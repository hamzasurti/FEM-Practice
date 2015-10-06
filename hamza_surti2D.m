function DemoFEM2DLaplaceEssentialBoundaryCondition
%
% This function demonstrates how the 2D Laplace equation with prescribed 
% essential boundary conditions is solved by the finite element method.  
%

clc, clear, close all;

% geometric parameters of a rectangule domain

Length = 90;
Height = 30;
ElementSize = 1;

% temperature at the boundary

TemperatureAtTop    = 300;
TemperatureAtBottom = 260;
TemperatureAtLeft   = 260;
TemperatureAtRight  = 260;

% number of elements

NumberOfElementInX = 2*(floor(0.5*Length/ElementSize)+1); %these two linesdetermine how many nodes there will be based on the size of the problem
NumberOfElementInY = ceil(Height/ElementSize);

% Element size

% ElementSizeInX = Length/NumberOfElementInX;
% ElementSizeInY = Height/NumberOfElementInY;

% construct mesh data

NumberOfNode    = (NumberOfElementInX+1)*(NumberOfElementInY+1); %sets up nodes
NumberOfElement = 2*NumberOfElementInX*NumberOfElementInY;

% node coordinates
% tic;
% X = zeros(NumberOfNode,1);
% Y = zeros(NumberOfNode,1);
% 
% for j = 1 : NumberOfElementInY+1
%     for i = 1 : NumberOfElementInX+1
%         n = (j-1)*(NumberOfElementInX+1)+i;
%         X(n) = (i-1)*ElementSizeInX;
%         Y(n) = (j-1)*ElementSizeInY;
%     end
% end
% toc;

tic;
x = linspace(0, Length, NumberOfElementInX+1); %These lines set up the array
y = linspace(0, Height, NumberOfElementInY+1);
[X, Y] = meshgrid(x,y);
X = X';
Y = Y';
X=X(:);
Y=Y(:);
toc;

% I = randperm(length(X));
% X=X(I);
% Y=Y(I);

% element node connectivity

% ElementNode = zeros(NumberOfElement,3);
% 
% for j = 1 : NumberOfElementInY
%     for i = 1 : NumberOfElementInX
%         n = (j-1)*NumberOfElementInX+i;
%         if i <= NumberOfElementInX/2
%             ElementNode(2*n-1, 1) = (j-1)*(NumberOfElementInX+1)+i;
%             ElementNode(2*n-1, 2) = j*(NumberOfElementInX+1)+i+1;
%             ElementNode(2*n-1, 3) = j*(NumberOfElementInX+1)+i;
%             ElementNode(2*n, 1)   = (j-1)*(NumberOfElementInX+1)+i;
%             ElementNode(2*n, 2)   = (j-1)*(NumberOfElementInX+1)+i+1;
%             ElementNode(2*n, 3)   = j*(NumberOfElementInX+1)+i+1;
%         else
%             ElementNode(2*n-1, 1) = (j-1)*(NumberOfElementInX+1)+i;
%             ElementNode(2*n-1, 2) = (j-1)*(NumberOfElementInX+1)+i+1;
%             ElementNode(2*n-1, 3) = j*(NumberOfElementInX+1)+i;
%             ElementNode(2*n, 1)   = (j-1)*(NumberOfElementInX+1)+i+1;
%             ElementNode(2*n, 2)   = j*(NumberOfElementInX+1)+i+1;
%             ElementNode(2*n, 3)   = j*(NumberOfElementInX+1)+i;
%         end
%     end
% end

ElementNode = delaunay(X,Y-0.001*ElementSize*sin(pi*X/Length));
% h=trimesh(ElementNode,X,Y, 'Color', 'k', 'LineWidth', 3);
% get(h);
% axis equal

% boundary condition

Temperature  = zeros(NumberOfNode,1); %set up empty array
NodeEquation = ones(NumberOfNode,1);

for j = 1 : NumberOfElementInY+1
    for i = 1 : NumberOfElementInX+1
        n = (j-1)*(NumberOfElementInX+1)+i;
        if j == 1
            Temperature(n) = TemperatureAtBottom;
            NodeEquation(n) = 0;
        end
        if j == NumberOfElementInY+1
            Temperature(n) = TemperatureAtTop;
            NodeEquation(n) = 0;
        end
        if i == 1
            Temperature(n) = TemperatureAtLeft;
            NodeEquation(n) = 0;
        end
        if i == NumberOfElementInX+1
            Temperature(n) = TemperatureAtRight;
            NodeEquation(n) = 0;
        end
    end
end

% number of equations and node corresponding equation

NumberOfEquation = length(NodeEquation(NodeEquation==1));
NodeEquation(NodeEquation==1) = 1 : NumberOfEquation;

% solve temperature using the finite element method

Temperature = FEM2DLaplaceEssentialBoundaryCondition(X, Y, Temperature, NodeEquation, ElementNode);

% analytical solution

ExactTemperature =zeros(size(Temperature));
n = 0;
Error = 1;
while Error > 1e-10
    n = n + 1;
    Term = (sin((2*n-1)*pi*X/Length).*sinh((2*n-1)*pi*Y/Length)/((2*n-1)*sinh((2*n-1)*pi*Height/Length)));
    ExactTemperature = ExactTemperature + Term;
    Error = max(Term)/max(ExactTemperature);
end
ExactTemperature = ExactTemperature*4*TemperatureAtTop/pi;

% plot finte element and analytical results for comparison

ScreenSize = get(0,'ScreenSize');
figure('Name',' DemoFEM2DLaplaceEssentialBoundaryCondition', ...
       'NumberTitle','off', ...
       'OuterPosition',[5, 1/3*ScreenSize(4) 4/3*ScreenSize(4) 2/3*ScreenSize(4)]);

subplot(1,2,1);
patch('Faces',ElementNode,'Vertices',[X Y],'FaceColor','interp','FaceVertexCData',Temperature,'EdgeColor','none'); 
axis equal  
colorbar('location','westoutside')
xlabel('X'); 
ylabel('Y');  
title('FEM Solution');

subplot(1,2,2);
patch('Faces',ElementNode,'Vertices',[X Y],'FaceColor','interp','FaceVertexCData',ExactTemperature,'EdgeColor','none'); 
axis equal  
colorbar('location','Eastoutside')
xlabel('X'); 
ylabel('Y');  
title('Exact Solution');


function T = FEM2DLaplaceEssentialBoundaryCondition(X, Y, T, NodeEquation, ElementNode)

% solve the following Laplace equation using FEM
% T,xx + T,yy = 0
% T = g on boundary

% location matrix

LocationMatrix = NodeEquation(ElementNode);

% number of equations and elements

NumberOfEquation = max(NodeEquation);
NumberOfElement  = size(ElementNode,1);

% assemble system matrix and vector

SystemMatrix = zeros(NumberOfEquation);
SystemVector = zeros(NumberOfEquation,1);

for n = 1 : NumberOfElement
    
    % element coordinates
    
    XE = X(ElementNode(n,:));
    YE = Y(ElementNode(n,:));
    
    % derivatives of shape functions
    
    Area = 0.5*(XE(2)-XE(1))*(YE(3)-YE(1))-(YE(2)-YE(1))*(XE(3)-XE(1));
    
    NX = [YE(2)-YE(3); YE(3)-YE(1); YE(1)-YE(2)]/(2*Area);
    NY = [XE(3)-XE(2); XE(1)-XE(3); XE(2)-XE(1)]/(2*Area);
    
    % element matrix
    
    ElementMatrix = (NX*NX'+NY*NY')*Area;
    
    % system matrix and vector
    
    for i = 1 : 3
        ii = LocationMatrix(n,i);
        if ii~=0
            for j = 1 : 3
                jj = LocationMatrix(n,j);
                if jj~=0
                    SystemMatrix(ii,jj) = SystemMatrix(ii,jj) + ElementMatrix(i,j);
                else
                    SystemVector(ii)    = SystemVector(ii) - ElementMatrix(i,j)*T(ElementNode(n,j));
                end
            end
        end
    end
    
end

% solve system matrix equations
tic
SystemVariable = sparse(SystemMatrix)\SystemVector;
toc
B=SystemMatrix~=0;
alpha = length(B(:))/length(SystemVector)^2

% assign the solution to T

T(NodeEquation>0) = SystemVariable;

% In the beginning, you set up the problem based on the parameters given.
% We had to change the boundary conditions to those that matched the problem
% that we were doing. We also set up our conditions so that matlab could figure
% out how many elements we have in the x and y directions. From that we were
% able to get our total number of nodes and our total number of elements.
% Next we created a linearly spaced vector and Professor Xu added elements to
% see how quickly the computer could do that calculations. Then we used the
% deluanay triangulation function which created the triangles we need to solve
% the problem. Then the temperature vector is initialized and loops are used
% to insert the boundary conditions into the vectors. The number of equations
% and the corresponding total equation is then calculated and temperature is
% solved for using the finite element methods that we set up. The next while
% loop helps us calculate the the exact solutions so that we can compare it
% to what we got using FEM. After that everything is just plotted on a heat
% map. The next function sets up the location matrix and solves it by 
% initializing and filling new arrays. The final for loop does most of the
% matrix calculations in order for us to get our system variable. This allows
% us to solve our final matrix equations and get our final answer.
