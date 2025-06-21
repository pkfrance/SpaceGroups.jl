# SpaceGroups.jl

**SpaceGroups.jl** is a lightweight Julia package for handling crystallographic symmetry operations in any spatial dimension.

---

## Design Principles

**SpaceGroups.jl** is built on the following principles:

* **Lightweight and Self-Contained:** Depends only on Julia's standard library.
* **Dimension-Agnostic:** All algorithms are designed to work generically in any dimension.
* **Performance-Focused:** Uses `StaticArrays` under the hood.
* **Exact Representation:** Symmetry operations are represented precisely using Julia's `Rational` type for fractional components, avoiding the inaccuracies of floating-point numbers.

## What SpaceGroups.jl is Not

* **It is not a computer algebra system.** For tasks like enumerating all non-equivalent space groups in a given dimension, established tools like GAP are more suitable.
* **It is not a space group database.** If you need a comprehensive database of crystallographic space groups, consider such excellent packages as [Crystalline.jl](https://github.com/thchr/Crystalline.jl) or [Spglib.jl](https://github.com/singularitti/Spglib.jl).

## Core Concepts

### Space Group Quotient

The central type in the package is `SpaceGroupQuotient`, which represents the quotient of a space group by its translational subgroup. While this group is always isomorphic to the point group of the structure, the `SpaceGroupQuotient` retains sufficient information to reconstruct the full space group.

* **Construction:** A `SpaceGroupQuotient` is instantiated from its generators.
* **Eager Instantiation:** The full set of group elements is computed and stored upon creation.
    * ⚠️ **Note:** This approach may be memory-intensive in high-dimensional settings.

### Wyckoff Positions

The package provides a dedicated `WyckoffPosition` type to represent symmetry-equivalent classes of atomic sites within a crystal structure.

Key functionalities include:
* Applying space group operations to sites.
* Computing site-symmetry groups.
* Validating that a Wyckoff position is compatible with a given space group.

### Bragg Peaks

For any given point in the reciprocal lattice, the package can analyze its symmetry orbit to determine whether it corresponds to:

* A **Bragg peak** with an arbitrary phase.
* A **Bragg peak** with a phase restricted to a sign difference.
* An **extinct peak** (forbidden by symmetry).

---

## Basic Usage

Here is an example of defining a non-symmorphic primitive octagonal space group in 4D and checking for systematic extinctions:

```julia-repl
julia> using SpaceGroups

# Define the generators for the space group
julia> g1 = @SGE([0 0 0 -1; 1 0 0 0; 0 1 0 0; 0 0 1 0]);
julia> g2 = @SGE([0 1 0 0; 1 0 0 0; 0 0 0 -1; 0 0 -1 0], [1//2, 1//2, 1//2, 1//2]);

# Create the space group quotient
julia> G = SpaceGroupQuotient([g1, g2])
SpaceGroupQuotient (dimension 4, order 16)

# A reflection with restricted phase
julia> make_orbit([1, 1, 0, 0], G)
RealOrbit with 8 elements

# An extinct reflection
julia> make_orbit([1, 0, 0, 0], G)
ExtinctOrbit with 8 elements
```

## Installation

This package is not yet in the official Julia General registry. You can install it using the `Pkg` REPL (press `]` in the Julia REPL to enter `Pkg` mode) with:

```julia-repl
julia> ]
(@v1.10) pkg> add https://github.com/pkfrance/SpaceGroups.jl(https://github.com/pkfrance/SpaceGroups.jl)
```

## Documentation
 - [DEV](https://pkfrance.github.io/SpaceGroups.jl/dev/): documentation for the development version.
 - STABLE: *work in progress", waiting for the first release.

## Contributing
Contributions, bug reports, and feature suggestions are welcome! Feel free to open an issue or a pull request.