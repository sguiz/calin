# CALIN: Composable Asynchronous Logic Integrase Networks

CALIN permits automatization of design of integrase based logic gates, either for Boolean and History-dependent programs. Designs follow the framework corresponding to the pre-print paper, [Guiziou et al., BiorXiv 2017](https://doi.org/10.1101/150987). The design of Boolean programs have been published on ACS Synth. Bio. : [Guiziou et al., ACS Synth. Bio.](https://doi.org/10.1021/acssynbio.8b00016)
These scripts permits the generation of the genetic layout of the biological design (using dnaplotlib, [Der et al, ACS Synbio 2016](http://pubs.acs.org/doi/abs/10.1021/acssynbio.6b00252) and the DNA sequence.

A web-interface using as source code this code is available at http://synbio.cbs.cnrs.fr/calin/.

## Installing

CALIN code requires Python 2.7 and installation of DNAplotlib, qm, biopython and matplotlib 1.2 or newer. 
- For biopython, please download the module [here](http://biopython.org/wiki/Download) and to install: `pip install biopython`.
- For DNAplotlib, please download the python code [here](https://github.com/VoigtLab/dnaplotlib) and add the file to your calin script directory.
- For qm, please download the python code [here](https://pypi.python.org/pypi/qm) and add the file to your calin script directory.

To use CALIN, download all python files and place them to the same diectory than DNAplotlib and qm.

## Getting Started

For execution of the Python codon such as in the CALIN web-interface, [main.py](https://github.com/sguiz/calin/blob/master/main.py) should be excuted with the following arguments: 
- the type of logic, for Boolean logic: 'comb' and for history-dependent logic: 'seq'.
- the number of inputs.
- the logic function to implement.
- the output folder.
Please first create a results folder in your CALIN python directory which should contain as well your script directory. The results are saved in ../results/output-folder.

### Example for Boolean logic

To obtain the biological design of the boolean function: f=A.B.C+not(A).not(B).not(C), 
please execute the code in the terminal as following: `python directory/main.py 'comb' '2' '10001' 'example_boolean'`

A example_boolean folder will be create in your results directory with 2 genbank files, and 2 image files (see below). All results are available [here](https://github.com/sguiz/calin/tree/master/results/example_boolean).

<img src="https://github.com/sguiz/calin/blob/master/results/example_boolean/example_boolean_Strain1.png" width="400"> <img src="https://github.com/sguiz/calin/blob/master/results/example_boolean/example_boolean_Strain2.png" width="400">

### Example for history-dependent logic

To obtain the biological design of the following history-dependent function:     <img src="https://github.com/sguiz/calin/blob/master/results/example_seq/seq_function.png" width="80">

Please execute the code in the terminal as following:
`python directory/main.py 'seq' '2' '32001' 'example_seq'`

A example_seq folder will be create in your results directory with 2 genbank files, and 2 image files (see below). All results are available [here](https://github.com/sguiz/calin/tree/master/results/example_seq).

<img src="https://github.com/sguiz/calin/blob/master/results/example_seq/example_seq_Strain1.png" width="400"> <img src="https://github.com/sguiz/calin/blob/master/results/example_seq/example_seq_Strain2.png" width="400">

## Author

**Sarah Guiziou**

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

* **Jerome Bonnet** for adivces
* Ashley Nord for tips on Python

