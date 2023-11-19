# julia_educational_examples
Some simple scripts to show the trickier parts of Julia

Run the script in VSCode with the Julia extension.

The reason it is not a notebook is to more clearly demonstrate global and local scope.

To load the sparse array into Julia encoding the Fluent mesh in "halfpipe_tri.msh", use:
cell_sparse = load("halfpipe_tri.jld2","cell_sparse")
