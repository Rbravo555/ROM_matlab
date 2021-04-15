%script for running isomap 
clear all;
close all;
clc;

%load data and plot
load swiss_roll_data
scatter3(X_data(1,:),X_data(2,:),X_data(3,:),[],Y_data'*ones(size(Y_data,1),1) , 'filled');
colormap jet(256);
view(-15,20); axis('equal'); axis('off');

%calculate the L2 distance and run isomap
D = L2_distance(X_data(:,1:1000), X_data(:,1:1000), 1);
options.dims = 1;
[Y, R, E] = Isomap(D, 'k', 7, options);
%scatter(Y.coords{1}(1,:), Y.coords{1}(2,:),[],Y_data(:,1:1000)'*ones(size(Y_data(:,1:1000),1),1), 'filled');




