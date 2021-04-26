function VisualizeSelectedElements(elements, X)
%naive way of visualizing selected elements in a 1D mesh
    %the selected elements will be shown using the array of elements
    X(  setdiff([1:length(X)],elements) ) = [];

    
    figure(88)
    hold on
    plot(X,zeros(size(X)), 'x', 'LineWidth',2)
    title('Selected elements in 1D mesh')

end