# About

License, acknowledgements, and contributors for Q2D.  This section mirrors the MOSE `about/` documentation 1:1.

## License

Q2D is free and open-source software released under the **GNU General Public License v3.0**. See `LICENSE` in the repository root for the full license text.

## Acknowledgements

Q2D shares the heritage of the [`open-hydra`](https://github.com/open-hydra) release with [MOSE](https://github.com/open-hydra/MOSE), the companion 3D solver. Q2D was forked from an early version of MOSE and evolved independently with focus on three geometric regimes (2D planar, 2D axisymmetric, Q2D area-variation), Q2D-specific numerical features, and applications to hybrid rocket motors (paraffin grain modeling) and Borda-flow choked nozzles (Fievisohn 2016/2017).

The structural refactoring to mirror MOSE was completed in 2026, restoring architectural parity while preserving Q2D's post-fork additions.

## Companion Libraries

Q2D builds on the same companion library stack as MOSE:

- [FLINT](https://github.com/MarcoGrossi92/FLINT) — Thermodynamic and chemical kinetics database.
- [ORION](https://github.com/MarcoGrossi92/ORION) — Multi-format I/O (Tecplot, VTK, Plot3D).
- [OSLO](https://github.com/MarcoGrossi92/OSLO) — ODE solver library (stiff chemistry integration).
- [FiNeR](https://github.com/szaghi/FiNeR) — INI configuration file parser.

## Contributors

The Q2D core team and contributors are listed in the `AUTHORS` file at the repository root.

For citation: see `CITATION.cff` at the repository root.
