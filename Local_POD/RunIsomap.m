function embedding = RunIsomap( HighDimensionalRepresentation, dimensions_of_embedding, number_of_neighbors);

D = L2_distance(HighDimensionalRepresentation,HighDimensionalRepresentation,1);
options.dims = dimensions_of_embedding;
options.display = 0;
options.overlay = 0;
options.verbose = 0;
[Y, R] = Isomap(D, 'k', number_of_neighbors, options);
embedding = Y.coords{1};
end