close all
clear all
clc
%%
% *choosing parameters for genetic algorithm* 

p = 100; % population size base on genoms
cross_num = 20;
mut_num = 20;
selection_num = 60;

%%
% *creating cities coordinates* 

x = 5; y = 5; z = 5;
% pick 20 of random points
x = [0 x*rand(1,19)]; %x-coordinate of a point
y = [0 y*rand(1,19)]; %y-coordinate of a point
z = [0 z*rand(1,19)]; %z-coordinate of a point

%inspect cities locatins
figure('Name','Cities');
labels = ('A' : 'T')';
scatter3(x, y, z, 'MarkerFaceColor','b','MarkerEdgeColor','b');
text(x', y', z',labels);
title('Cities locations');

%%
% *assign initial random trail*

labels = [labels;'A'];
x_plot = [x';0];
y_plot = [y';0];
z_plot = [z';0]; 
cities =  [x' y' z'];
figure('Name','path');
plot3(x_plot, y_plot, z_plot);
text(x_plot, y_plot, z_plot,labels);
title('first trail before evolution');

%%
% *creating genetic algorithm loop*

for gen = [1 100 500]    
    pop = population(p,20); %initialize population
    cost = zeros(1,gen);    
    for i = 1: gen     
        %compute cost and evaluation
        E = evaluation(pop,cities);
        cost(1,i) = min(E);        
        prev_pop = pop; %save previous population
        
        %make new generation
        pop(1:cross_num,:) = crossover(prev_pop,cross_num);
        pop(cross_num+1:cross_num + mut_num ,:) = mutation(prev_pop,mut_num);
        pop(cross_num + mut_num+1:p,:) = selection(prev_pop,E,selection_num); 
    end

    %plot cost
    figure('Name','Cost');
    plot(cost);
    cost_label = sprintf('cost after generation %d = %f',gen,cost(1,gen));
    title(cost_label);
    xlabel('Generation');
    ylabel('Fitness funcion - distance of cities');

    %plot optimal path
    [x_opt y_opt z_opt]= changeCordinates(pop,cities);
    figure('Name','path');
    plot3(x_opt, y_opt, z_opt);
    text(x_plot, y_plot, z_plot,labels);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    path_label = sprintf('path after generation %d',gen);
    title(path_label);
end

% population func
function Y=population(n,nc)
% nc = number of cities
% n = pop size
B=zeros(n,nc);
for j=1:n
    A=1:1:nc;   
    [x y]=size(A);
    for i = 1:nc
        r = randi(y);
        B(j,i)=A(r);
        A(:,r) = [];
        [x y]=size(A);
    end
end
Y = B;
end

% evaluation function
function YY = evaluation(P,cities)
% P = population
[n ,nc] = size(P);
Y = zeros(n,1);
for i = 1:n
    A = P(i,:); %get chromosome
    B = zeros(size(A));
    for j = 1:nc-1
        x_diff = (cities(A(1,j),1)-cities(A(1,j+1),1))^2;
        y_diff = (cities(A(1,j),2)-cities(A(1,j+1),2))^2;
        z_diff = (cities(A(1,j),3)-cities(A(1,j+1),3))^2;
        B(1,j) = sqrt(x_diff + y_diff + z_diff);
    end
        x_diff = (cities(A(1,1),1)-cities(A(1,nc),1))^2;
        y_diff = (cities(A(1,1),2)-cities(A(1,nc),2))^2;
        z_diff = (cities(A(1,1),3)-cities(A(1,nc),3))^2;
        B(1,j+1) = sqrt(x_diff + y_diff + z_diff);    
    Y(i,1) = sum(B);
end
YY = Y;
end

% Crossover function
function YY=crossover(X,n)
% X = population
% n = number of chromosomes to be crossed
[x1 y1] = size(X);
Y = zeros(n,y1);
for z = 1:n
    B = X(randi(x1),:); % select parent chromosome
    r1 = 1 + randi(y1-1);
    C = B(1,1:r1);
    B(:,1:r1) = [];  % cut
    [x3 y3] = size(B);
    B(1,y3+1:y1) = C;
    Y(z,:) = B;
end
YY = Y;
end

% mutation function
function YY = mutation(X,n)
% X = population
% n = number of chromosomes to be mutated
[x1 y1]=size(X);
Y=zeros(n,y1);
for z=1:n
    A=X(randi(x1),:); % select parent chromosome
    r1=1+randi(y1-1,1,2);
    while r1(1)==r1(2) 
            r1=1+randi(y1-1,1,2);
    end
    B=A(1,r1(1));
    A(1,r1(1))=A(1,r1(2));
    A(1,r1(2))=B;
    Y(z,:)=A;
end
YY = Y;
end

% Selection function
function YY = selection(P,E,selection_num)
% P - popolation
A = zeros(selection_num,20);
copy_E = E;
for i = 1:selection_num
    [row c] = find(copy_E == min(copy_E));
    row = row(1,1);

    A(i,:) = P(row,:);
    copy_E(row,1) = 0;
end
 YY = A;
end

% change coordination function
function [x,y,z] = changeCordinates(pop,cities)
E = evaluation(pop,cities);
[row c] = find(E == min(E));
row = row(1,1);
A = pop(row,:);
B = zeros(size(cities));
for i = 1:20
    B(i,:) = cities(A(1,i),:);
end
x = B(:,1);
y = B(:,2);
z = B(:,3);

x = [x;x(1)];
y = [y;y(1)];
z = [z;z(1)];
end
