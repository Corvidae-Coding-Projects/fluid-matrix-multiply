# Fluidic matrix multiplier

Explore a 4×4 fluidic matrix multiplier through computed three-dimensional
pressure and velocity fields.

**Live interactive demo:** https://corvidae-coding-projects.github.io/fluid-matrix-multiply/

Rotate and zoom the device, switch between computed CFD runs, inspect pressure
or speed, show velocity arrows, and slice the geometry by height. The viewer
shows bundled simulation results; choosing a case does not run a new simulation.

The complete source, validation results, and volume data are in
[`fluidic-cfd/`](fluidic-cfd/). See the [method and reproduction
instructions](fluidic-cfd/README.md) for the model, limitations, numerical checks,
and how to run a new experiment locally.

To verify the bundled source and results:

```sh
cd fluidic-cfd
python3 verify-package.py
```

You can also open `fluidic-cfd/fluidic-matrix-multiplier-3d.html` directly in a
modern browser. It bundles its viewer library and computed data for offline use.
WebGL and the standard DecompressionStream API are required.

GitHub Pages automatically publishes the bundled standalone viewer whenever
`main` changes. The deployment checks the package manifest before publishing.
