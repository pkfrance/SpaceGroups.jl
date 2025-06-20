# SpaceGroups.jl

**SpaceGroups.jl** is a lightweight Julia package for working with **crystallographic space groups in arbitrary dimensions**. It provides efficient data structures for handling crystallographic symmetry operations, emphasizing performance and minimal dependencies.

---

## Features

- **Arbitrary-dimensional space groups**: Work with space groups in any dimension.
- **Performance-focused**: Uses static arrays internally for speed and memory efficiency.
- **Lightweight and self-contained**: Depends only on Julia's standard library.
- **Numerical focus**: Optimized for numerical computations where performance and dimensional generality are essential. *Note*: This is not a symbolic algebra tool and does not aim to match the scope of packages like GAP.
- **Custom integer types**: All types are parameterized by the underlying integer type (default: `Int`, but supports arbitrary-precision integers and other integer types).

---

## Core Concepts

### Space Group Quotient

The central type in the package is `SpaceGroupQuotient`, representing the **quotient of a space group by its translational subgroup**. While this group is always isomorphic to the point group of the structure, the `SpaceGroupQuotient` retains sufficient information to **reconstruct the full space group** as a semidirect product.

- Constructed from a list of generators.
- **Eager instantiation**: Computes and stores the full set of group elements upon creation.
  - ⚠️ *Note*: May be memory-intensive in high-dimensional settings.

### Wyckoff Position

The package provides a dedicated type, `WyckoffPosition`, representing symmetry classes of atomic sites in a crystal structure.  

Key functionalities include:
- Applying space group operations to these sites.
- Computing **site-symmetry groups**.
- **Validating** Wyckoff positions against the space group.

### Bragg Peaks

For a given point in the reciprocal lattice, the package determines whether it corresponds to:
- A **Bragg peak** with an arbitrary phase.
- A **Bragg peak** with a phase fixed up to a sign.
- An **extinct peak** (forbidden by symmetry).

---

## Installation

Until the package is officially registered, install it via the Julia `Pkg` REPL (press `]` in the Julia REPL to enter `Pkg` mode) with:

```julia-repl
julia> ]
(@v1.10) pkg> add SpaceGroups
```

## Documentation
*Work in progress*. For now, please refer to the source code and docstrings.

## Contributing
Contributions, bug reports, and feature suggestions are welcome! Feel free to open an issue or a pull request.