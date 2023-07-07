# fontanka

Fountains are novel features of Hi-C maps, resembling "flares" [[1](#ref1)], "plumes" [[2](#ref2)] and "jets" [[3](#ref3)]. They can be detected as inverted triangular structures along the main diagonal of Hi-C maps.
Fountains are also similar to "hairpins" [[4](#ref4)]. 

Fontanka aims to provide a flexible Python API and specialized CLI for calling fountains from cool Hi-C files. 

Note that Chromosight [[4](#ref4)] can be modified to do fountains calling with fountain mask. 

## Installation

It is highly recommended to work in conda environment, which can be found in environment.yml:
```bash
conda env create -f environment.yml
conda activate fontanka
```

Optionally, you may want to activate IPython kernel to work with API:
```bash
conda install ipykernel
python -m ipykernel install --user --name fontanka --display-name "fontanka"
 ```

```bash
git clone https://github.com/agalitsyna/fontanka.git
cd fontanka
pip install -e ./
```

## Example usage

Call for the binary fountain mask:

```bash
fontanka call-fountains test.mcool::resolutions/10000 \
         output.fountains.tsv \
         -A 0.7854 \
         -W 200_000 \
         -p 20 \
         --store-snips output.snips.fountains.npy \
         --regions chromosome.regions.txt 
```

For running fontanka for Galitsyna et al. 2023, see the notebook in the "examples/".

### References: 

1. <a name="ref1" href="https://www.genome.org/cgi/doi/10.1101/gr.269860.120">Wike, C. L., Guo, Y., Tan, M., Nakamura, R., Shaw, D. K., Díaz, N., ... & Cairns, B. R. (2021). 
   Chromatin architecture transitions from zebrafish sperm through early embryogenesis. 
   Genome Research, 31(6), 981-994.</a>

2. <a name="ref2" href="https://doi.org/10.1101/2021.08.27.457977">Liu, N. Q., Magnitov, M., Schijns, M., van Schaik, T., van der Weide, R. H., Teunissen, H., ... & de Wit, E. (2021). 
   Rapid depletion of CTCF and cohesin proteins reveals dynamic features of chromosome architecture.
   </a>

3. <a name="ref3" href="https://doi.org/10.1016/j.molcel.2022.09.003">Guo, Y., Al-Jibury, E., Garcia-Millan, R., Ntagiantas, K., King, J.W., Nash, A.J., Galjart, N., Lenhard, B., Rueckert, D., Fisher, A.G. and Pruessner, G., 2022. Chromatin jets define the properties of cohesin-driven in vivo loop extrusion. Molecular Cell, 82(20), pp.3769-3780.</a>


4. <a name="ref4" href="https://doi.org/10.1038/s41467-020-19562-7">Matthey-Doret, C., Baudry, L., Breuer, A., Montagne, R., Guiglielmoni, N., Scolari, V., ... & Cournac, A. (2020). 
   Computer vision for pattern detection in chromosome contact maps. 
   Nature Communications, 11(1), 1-11.</a>