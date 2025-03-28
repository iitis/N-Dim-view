# N-Dim-view plugin for dpVision

**N-Dim-view** is a plugin for the open-source platform [dpVision](https://github.com/pojdulos/dpVision) that provides an interactive 3D visualisation of multidimensional data using expressive avatars. The plugin combines spatial projection with glyph-based visual metaphors to help users explore complex datasets in an intuitive and visually engaging way.

## 🧠 Key Features

- **Chernoff-inspired avatars** – Each data point is visualised as a small figure ("avatar") whose facial traits represent selected data dimensions.
- **3D spatial projection** – Additional dimensions are mapped to position in a 3D coordinate system, allowing for interactive navigation and inspection.
- **Custom slab slicing** – Users can control the thickness of the visible data layer in nD space, focusing on cross-sections of interest.
- **Manual or automatic assignment** – Flexible mapping of data dimensions to avatar appearance or spatial coordinates.
- **Support for real and synthetic data** – Tested with a 15D wine quality dataset and an artificial 4D tesseract-shaped dataset.

## 📦 Installation

The plugin is designed to work with `dpVision` version X.X or higher.

1. Clone this repository:
    ```bash
    git clone https://github.com/iitis/n-dim-view.git
    ```

2. Copy the plugin folder to the `dpVision/plugins/` directory.

3. Reconfigure your dpVision project with CMake or add a plugin project to dpVision solution in MS Visual Studio

4. Build dpVision project 

5. Launch `dpVision` and select the plugin in the plugin manager.

6. Use plugin to load your dataset.

<!-- 

## 📊 Input Format

The plugin expects a dataset in the form of a matrix:
- Rows: dimensions (features)
- Columns: data points

A separate configuration file (optional) can be used to:
- assign dimensions to facial traits (e.g., eye shape, mouth curvature)
- assign dimensions to spatial coordinates (e.g., X, Y, Z)
- define viewing parameters, PCA settings, or label colours

-->

## 🚀 Example Use Cases

- Visualising clusters in high-dimensional data
- Exploring correlations between human-perceived qualities and technical parameters
- Educational visualisation of abstract geometries (e.g., 4D rotation)
- Dimensionality reduction with interpretability

<!--
## 📷 Screenshots


![Swarm example](docs/img/swarm_example.png)
*A swarm of avatars based on wine quality data.*

![Tesseract rotation](docs/img/tesseract_rotating.gif)
*A rotating 4D anisotropic hypercube projected into 3D.*

## 📘 Related Publications

If you use this plugin in your research, please cite:

> Leszek Luchowski, Dariusz Pojda.  
> *Visualisation of a Multidimensional Point Cloud as a 3D Swarm of Avatars – A Plugin for dpVision.*  
> IEEE Transactions on Visualization and Computer Graphics, 2025.
-->

## 🛠️ Roadmap

Planned features include:
- Support for alternative glyph styles (trees, neutral geometric forms)
- Heuristic or ML-based automatic dimension-to-feature mapping
- Scalability enhancements for large datasets
- Export to interactive WebGL viewers
- User feedback and usability testing

## 🧑‍💻 Authors

- **Leszek Luchowski** – Concept, algorithm design, geometry
- **Dariusz Pojda** – Implementation, system integration, optimisation

## 📄 License

This plugin is released under the LGPL-2.1 License. See [LICENSE](LICENSE) for details.
