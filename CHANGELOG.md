# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased](unreleased)
- Removed: Support for Julia v0.5

## [0.2.0](2015-09-26)
- Added: Initial support for efficient discrete difference calculation, currently just l_2 norms
- Added: honesty mode, a flag which can be set per-repo to enable more-or-less accurate approximations to Lebesgue integrals (true) or regular discrete norms (false).
- Added: Support for Julia v0.4, and simplified/faster calculation with LinSpace grids

## [0.1.0](2015-04-18)
- Added: Initial implementations of these vector norms, based on ISP.jl/PDE.jl versions
- Changed: Performance of implemented routines is now improved due to devectorization of the relevant routines
- Added: Completely split-up implementations of discrete Sobolev norms
- Added: Weighted matrix norms
- Added: Unit tests for relatively obvious sanity checks
