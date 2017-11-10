# CALIN: Composable Asynchronous Logic Integrase Networks

CALIN permits automatization of design of integrase based logic gates, either for Boolean and History-dependent programs. Designs follow the framework corresponding to the pre-print paper: https://doi.org/10.1101/150987. 
These scripts permits the generation of the genetic layout of the biological design (using dnaplotlib, [Der et al, ACS Synbio 2016](http://pubs.acs.org/doi/abs/10.1021/acssynbio.6b00252) and the DNA sequence.

A web-interface using as source code this code is available at http://synbio.cbs.cnrs.fr/calin/.

## Installing

CALIN code requires Python 2.7 and installation of DNAplotlib, biopython and matplotlib 1.2 or newer. 
For biopython, please download the module at http://biopython.org/wiki/Download and to install: 'pip install biopython'.
For DNAplotlib, please download the python code at https://github.com/VoigtLab/dnaplotlib and add the file to your calin directory.

To use CALIN, download all python files and place them to a single python diectory (with DNAplotlib).

## Getting Started

For execution of the Python codon such as in the CALIN web-interface, [main.py](https://github.com/sguiz/calin/blob/master/main.py) should be excuted with the following arguments: 
- the type of logic, for Boolean logic: 'comb' and for history-dependent logic: 'seq'.
- the number of inputs.
- the logic function to implement.
- the output folder.
Please first create a results folder in your CALIN python directory. The results are saved in ../results/output-folder.

### Example for Boolean logic

To obtain the biological design of the boolean function: f=A.B.C+not(A).not(B).not(C), please execute the code as folowing:

'python main.py 'comb' '2' '10001' 'XOR-gate' '

A XOR-gate folder will be create in your results directory with 2 genbank files, and 2 image files.

### Example for history-dependent logic


## Authors

* **Sarah Guiziou**

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* Jerome Bonnet for adivces
* Ashley Nord for tips on Python

