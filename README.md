<div id="top"></div>

<!-- PROJECT SHIELDS -->  
<div align="center">
  
  [![Documentation Status](https://readthedocs.org/projects/tensortrack/badge/?version=latest)](https://tensortrack.readthedocs.io/en/latest/index.html)
  [![CI](https://github.com/quantumghent/TensorTrack/actions/workflows/CI.yml/badge.svg)](https://github.com/quantumghent/TensorTrack/actions/workflows/CI.yml)
  [![Codecov](https://codecov.io/gh/quantumghent/TensorTrack/branch/main/graph/badge.svg?token=1I0XEB69TQ)](https://codecov.io/gh/quantumghent/TensorTrack)
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
  [![DOI](https://zenodo.org/badge/481924253.svg)](https://zenodo.org/badge/latestdoi/481924253)

</div>

<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/quantumghent/TensorTrack">
    <img src="docs/src/img/logo.png" alt="Logo" width="800">
  </a>

  <p align="center">
    An open-source tensor network library for MATLAB.
    <br />
    <a href="https://tensortrack.readthedocs.io/en/latest/"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://tensortrack.readthedocs.io/en/latest/examples/examples.html">View Demo</a>
    ·
    <a href="https://github.com/quantumghent/TensorTrack/issues">Report Bug</a>
    ·
    <a href="https://github.com/quantumghent/TensorTrack/issues">Request Feature</a>
  </p>
</div>

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#requirements">Requirements</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>


<!-- ABOUT THE PROJECT -->
## About The Project

<!-- [![Product Name Screen Shot][product-screenshot]](https://example.com) -->
This is a package which aims to efficiently implement the various elementary algorithms that arise in the context of tensor networks. Currently, this includes:

* Various basic and utility tensor operations, such as creation routines, various linear algebra routines (norms, traces, overlaps, ...), index permutations, ...
* Tensor contraction routines, both pairwise as well as through a network contraction routine.
* Tensor factorizations, such as QR, LQ and polar decompositions, eigendecompositions and singular value decompositions.
* Solver algorithms for eigen systems and linear systems.
 
Additionally, these tensors support a general global symmetries, in which case both memory and CPU usage are optimized. The framework is able to support both Abelian and non-Abelian symmetries, as well as symmetry groups with multiplicities, which can have bosonic or fermionic braiding rules.
The design of the algorithms is chosen such that the inclusion of symmetries should not alter the code after the creation of the tensors. Currently, the following symmetries are implemented:

* Z2
* U1
* SU2
* O2
* Direct product groups

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- GETTING STARTED -->
## Getting Started

### Requirements

This project depends on the following:
- MATLAB version R2020b or newer
- [Parallel Computing Toolbox](https://de.mathworks.com/products/parallel-computing.html)
- A C++ compiler compatible with your MATLAB version for MEX-file compilation.

### Installation
1. Clone the repo into a local folder.
  ```sh
  git clone https://github.com/quantumghent/TensorTrack.git mylocalfolder
  ```
  
2. Add the folder and subfolders to your MATLAB path.
  - Via the MATLAB UI:
      Home > Environment > Set Path > Add with Subfolders > mylocalfolder/src
  - Via the MATLAB Command Window:
    ```matlabsession
    addpath(genpath('mylocalfolder/src'))
    ```
  By default, the path is reset every time you close the application. You can permanently add this package to the path by calling ```savepath```.
  
3. Precompile the necessary mex files.
  Within the MATLAB Command Window, call:
  ```matlabsession
  GetMD5
  uninit
  ```
  
<p align="right">(<a href="#top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->
## Usage

_For examples, please refer to the [Documentation](https://tensortrack.readthedocs.io/en/latest/examples/examples.html)_

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- CONTRIBUTING -->
## Contributing

Contributions as well as feature requests are greatly appreciated.

If you have a suggestion that would make this project better, please do not hesitate to open an issue with the tag "enhancement". Alternatively, you could also fork the repo and create a pull request:

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

Don't forget to give the project a star! Thanks again!

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- CONTACT -->
## Contact
Lukas Devos - [lkdvos](https://github.com/lkdvos) - Lukas.Devos@ugent.be

Lander Burgelman - [leburgel](https://github.com/leburgel) - Lander.Burgelman@ugent.be

Project Link: [https://github.com/quantumghent/TensorTrack](https://github.com/quantumghent/TensorTrack)

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- ACKNOWLEDGMENTS -->
## Acknowledgments
This project has been made possible through the work of the following people:
* [Lukas Devos](https://orcid.org/0000-0002-0256-4200)
* [Lander Burgelman](https://orcid.org/0000-0003-1724-5330)
* [Bram Vanhecke](https://orcid.org/0000-0001-9557-1591)
* [Jutho Haegeman](https://orcid.org/0000-0002-0858-291X)
* [Frank Verstraete](https://orcid.org/0000-0003-0270-5592)
* [Laurens Vanderstraeten](https://orcid.org/0000-0002-3227-9822)
* ...


<p align="right">(<a href="#top">back to top</a>)</p>


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/github_username/repo_name.svg?style=for-the-badge
[contributors-url]: https://github.com/github_username/repo_name/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/github_username/repo_name.svg?style=for-the-badge
[forks-url]: https://github.com/github_username/repo_name/network/members
[stars-shield]: https://img.shields.io/github/stars/github_username/repo_name.svg?style=for-the-badge
[stars-url]: https://github.com/github_username/repo_name/stargazers
[issues-shield]: https://img.shields.io/github/issues/github_username/repo_name.svg?style=for-the-badge
[issues-url]: https://github.com/github_username/repo_name/issues
[license-shield]: https://img.shields.io/github/license/github_username/repo_name.svg?style=for-the-badge
[license-url]: https://github.com/github_username/repo_name/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/linkedin_username
[product-screenshot]: images/screenshot.png
