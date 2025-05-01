<div align="center">

<h1>Changes of cerebellar vascularization in a murine model of apnea of prematurity</h1>
 
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-blue.svg)](https://creativecommons.org/licenses/by/4.0/)
<!-- [![DOI](https://zenodo.org/badge/873064553.svg)](https://doi.org/???/zenodo.???) -->

</div>

> **Note**  
> This repository contains the data and R code for the _"Innovative 3D-image analysis of cerebellar vascularization highlights transcriptomic changes in a murine model of apnea of prematurity."_ paper

## ❔ Requirements

* R version 4.3 or newer
* R Studio version 2022.07 or newer

## 💻 Repository structure

* `analysis/`: R Markdown files containing the data analysis pipeline.
  * The first code chunk of any `.Rmd` file typically installs and loads all required packages based on the `renv.lock`.
* `data/`: Contains the raw and processed data used in the analyses.
* `src/`: R scripts declaring custom functions called within the analysis files (e.g., for plotting, data processing).
* `fig/`: Stores figures generated during the analysis.
* `_configs.yml`: Lists paths to external files (e.g., data) and potentially other configuration settings used within the code.
* `_dependencies.yml`: Lists the R packages required for this project (managed by `renv`, you don't need to edit this file).
* `renv.lock`: Records the exact versions of R packages used, ensuring reproducibility (managed by the `renv` package).
* `renv/`: Directory used by the `renv` package to manage project-specific libraries. You don't need to edit this folder.

## 📜 License

This project is licensed under the Creative Commons Attribution 4.0 International License - see the [LICENSE](LICENSE) file for details, or visit [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## 💬 Citation

* **Paper:** _TBA_
* **Code:** Agalic Rodriguez-Duboc, & Marc-Aurèle Rivière. (2025). agalic-rd/Vasc-AoP: Initial release (v0.1). Zenodo.

## ✨ Contributors

* **Agalic Rodriguez-Duboc**:  
[![ORCID](https://img.shields.io/badge/ORCID-A6CE39?style=flat-square&labelColor=white&logo=orcid&logoColor=A6CE39)][ORCID_ARD]
[![Research Gate](https://img.shields.io/badge/ResearchGate-00CCBB?style=flat-square&labelColor=white&logo=researchgate&logoColor=00CCBB)][RG_ARD]

* **Marc-Aurèle Rivière**:  
[![ORCID](https://img.shields.io/badge/ORCID-A6CE39?style=flat-square&labelColor=white&logo=orcid&logoColor=A6CE39)][ORCID_MAR]
[![Research Gate](https://img.shields.io/badge/ResearchGate-00CCBB?style=flat-square&labelColor=white&logo=researchgate&logoColor=00CCBB)][RG_MAR]

## 📫 Contact

For questions regarding the analyses, please contact: [**Agalic Rodriguez-Duboc**](mailto:agalic.rodriguez.duboc@ntnu.no?subject=Vasc%20AoP%20project)

For any other questions, please contact the corresponding author: [**Delphine Burel**](mailto:delphine.burel@univ-rouen.fr?subject=Vasc%20AoP%20project)

<!----------------------------------->

[RG_MAR]: https://www.researchgate.net/profile/Marc_Aurele_Riviere2
[ORCID_MAR]: https://orcid.org/0000-0002-5108-3382
[RG_ARD]: https://www.researchgate.net/profile/Agalic-Rodriguez-Duboc
[ORCID_ARD]: https://orcid.org/0000-0002-2084-3780
