# CALIN: Composable Asynchronous Logic Integrase Networks

CALIN permits automatization of design of integrase based logic gates, either for Boolean and History-dependent programs. Designs follow the framework corresponding to the pre-print paper: https://doi.org/10.1101/150987. 
These scripts permits the generation of the genetic layout of the biological design (using dnaplotlib, [Der et al, ACS Synbio 2016](http://pubs.acs.org/doi/abs/10.1021/acssynbio.6b00252) and the DNA sequence.

A web-interface using as source code this code is available at http://synbio.cbs.cnrs.fr/calin/.

## Installing

CALIN code requires Python 2.7 and installation of DNAplotlib, biopython and matplotlib 1.2 or newer. 
For biopython, please download the module at http://biopython.org/wiki/Download and to install: pip install biopython.
For DNAplotlib, please download the python code at https://github.com/VoigtLab/dnaplotlib and add the file to your calin directory.

To use CALIN, download all python files and place them to a single python diectory (with DNAplotlib).

## Getting Started

For use of the full CALIN property, [main.py](https://github.com/sguiz/calin/blob/master/main.py) should be excute as arguments: 
- the type of logic, for Boolean logic: 'comb' and for history-dependent logic: 'seq'.
- the logic function to implement.
- the output directory.
- the output folder.

### Example for Boolean logic


### Example for history-dependent logic


## Authors

* **Sarah Guiziou**

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* Jerome Bonnet for adivces
* Ashley Nord for tips on Python

