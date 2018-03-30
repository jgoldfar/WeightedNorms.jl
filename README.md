# WeightedNorms

[![Build Status](https://travis-ci.org/jgoldfar/WeightedNorms.jl.svg?branch=master)](https://travis-ci.org/jgoldfar/WeightedNorms.jl)
[![Coverage Status](https://coveralls.io/repos/jgoldfar/WeightedNorms.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/jgoldfar/WeightedNorms.jl?branch=master)
[![codecov.io](http://codecov.io/github/jgoldfar/WeightedNorms.jl/coverage.svg?branch=master)](http://codecov.io/github/jgoldfar/WeightedNorms.jl?branch=master)

### What is this repository for? ###

* For some applications, the natural way of measuring distance between two vectors or matrices is intrinsically dependent on the grid on which they are defined. In particular, precise estimates for the convergence of a particular method may be stated most simply using such a norm. Of course, these norms are equivalent to the "standard" Euclidean norm for vectors of a fixed length, but a grid-dependent norm will have different stability properties as the number of elements in the vector changes.

### How do I get set up? ###

    Pkg.clone("git@github.com:jgoldfar/WeightedNorms.jl.git")

### Usage ###

    using WeightedNorms

will make Lebesgue- and Sobolev-type discrete norms available in the current scope. More usage instructions TBD.

### Contribution guidelines ###

Contributions welcome! Submit an issue or pull request!

### Who do I talk to? ###

* Jonathan Goldfarb <jgoldfar@gmail.com>
