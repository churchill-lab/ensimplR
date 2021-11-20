## ensimplR

A simple R library that makes use of [ensimpl](https://churchilllab.jax.org/ensimpl).

## Installation

    devtools::install_github('churchill-lab/ensimplR', force=T)


## Functions

    releases(debug=FALSE)
> Information about the supported releases.

	searchGenes(term, species, release, debug=FALSE)
> Search for a gene.

    getGene(id, species, release, detail=FALSE, debug=FALSE)
> Get information about a particular gene.

    getGeneHistory(id, species, start_release, end_release, debug=FALSE)
> Get information about a particular gene from start_release to end_release.

    batchGenes(ids, species, release, details=FALSE, debug=FALSE)
> Get information for a list of genes.
