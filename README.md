# Structure Function (SF) Repository

## Description

This repository contains tools and algorithms for calculating the Structure Function (SF) to investigate turbulent properties of two-dimensional images, such as brightness temperature or column density images. The SF is utilized for understanding the spatial variability within these images, overcoming limitations seen in Fourier-transform methodologies, including edge effects. This repository provides a comprehensive solution for computing SF, along with adequate documentation and examples.

## Table of Contents

- [Background](#background)
- [Methodology](#methodology)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Background

The SF is an advantageous tool for analyzing the turbulence properties in two-dimensional images, offering an in-depth investigation without image size and shape constraints. Despite its limited spatial scale coverage compared to the Spectral Power Spectrum (SPS), the SF presents a robust analysis by working directly in the image domain. An in-depth analysis of the practical use of the SF can be found in [Mapping Spatial Variations of HI Turbulent Properties in the Small and Large Magellanic Cloud](https://iopscience.iop.org/article/10.3847/1538-4357/ab53df) found in The Astrophysical Journal.

## Methodology

1. **Image Integration and Normalization:**
   - Integrate the image over all velocity channels.
   - Normalize the integrated image by its maximum pixel value.
   
2. **Calculation of SF:**
   - Use the given Equation (5) to calculate the SF at each discrete pixel separation `r`, where `r′` represents any arbitrary pixel in the image and `I` represents either the normalized intensity or normalized column density of the pixel.
   - Calculate SF `SFr` for all realizations to estimate the mean `SFr` and the standard deviation over different realizations.
   
3. **Binning:**
   - Bin `r` values by 0.05 pixels in log space.
   - Utilize estimated uncertainties for each `r` to calculate the weighted average and standard deviation of a given bin.

4. **Special Cases:**
   - For specific cases, adjust the thickness of velocity channels (`Δv`), as detailed in the repository documentation.

## Installation

Ensure that the required software for running the SF algorithms is installed on your computer. Follow the instructions provided in the documentation to clone and set up the repository.

```bash
# Clone this repository
$ git clone https://github.com/user/SF-Repository

# Go into the repository
$ cd SF-Repository

# Install dependencies as listed in the documentation
```

## Usage

After installation, refer to the README files in each subdirectory for instructions on running the SF algorithms and examples. Ensure to adhere to the guidelines for selecting and adjusting parameters for accurate SF calculations.

## Contributing

Community contributions are welcomed to enhance the functionalities of the SF Repository:

- Fork the repository.
- Create a feature branch (`git checkout -b feature_branch`).
- Commit the changes (`git commit -am 'Add some feature'`).
- Push to the branch (`git push origin feature_branch`).
- Submit a Pull Request.

Before contributing, please read our [Contributing Guidelines](CONTRIBUTING.md).

## License

This project is under the MIT License. Refer to the [LICENSE.md](LICENSE.md) file for more details.

## Contact

For inquiries or concerns, please open an issue or contact the maintainers directly at email@example.com. Explore the SF Repository for insightful turbulent property analysis in two-dimensional images!
