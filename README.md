# Thesis methods

The main goal of this repository is to show the persons that read my PhD dissertation how the computational analyses were done. The scripts have commented code inside them, so everyone can follow what the scripts are doing. Also, the Materials and Methods section explains what every script does.



***



## Installation
The use of the scripts requires several programs in order to run properly and reproduce the results. In order to save the user the tedious installation processes and the searching of the programs manually, we have generated conda recipes. These conda recipes can be used directly by users in order to generate ready-to-analyse conda environments. The only thing the user has to change from these recipes are the prefix variables.

If you don't have Conda installed, we recommend using Anaconda (or Miniconda) to start. Check [this link](https://docs.conda.io/en/latest/miniconda.html) to install conda from Miniconda. 

In each folder of the different analyses, there are YML files that contain these recipes.

Creating the environment for the ATAC-seq would be like this 
```
conda env create -f ATAC_env.yml
```
Once you have created the environment, do not forget to activate the environment with:

```
conda activate ATAC_env
```

With these steps, you have loaded all the necessary programs to reproduce all our ATAC-seq analyses.


## Authors and acknowledgement
The main author of this repository is Alejandro Gil (agilgal AT upo DOT es agilgal@upo.es)
Some scripts and code come from diverse persons from my lab, like Juan J. Tena (ATAC_pipe.pl) or Alberto Perez Posadas.


## License
MIT License

Copyright (c) 2023 Alejandro Gil

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


