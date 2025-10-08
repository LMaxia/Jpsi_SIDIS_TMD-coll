# J/psi production at intermediate qT in SIDIS
The code generates a qT distribution of J/psi inclusive production in SIDIS, employing the NRQCD approach.
Here, qT corresponds to the photon transverse momenta, which is directly related to the PT of the Jpsi via the Jpsi energy fraction zh: |qT| = |PT|/zh.
Results are obtained for the lowest non-zero order in alpha_s, in both TMD and collinear factorizations.

Input variables are selected in the input.yaml file.
TMD.py and Collinear.py generate a .tsv file of the distributions in the TMD and collinear factorization, respectively.

The output includes variations if the option is activated in the input file.
In the header of the output file, only the parameters that deviate from their central values are explicitly shown.

If not included, PDF sets must be downloaded from the LHAPDF database to successfully run the code.

The factorization scale is chosen as mu^2 = M^2 + Q^2, where M is the J/psi mass and Q the photon virtuality.


Tips:
1.  Different regions of qT might require more points to produce high-quality figures. In the necessary case, modify the qT array in the YAML file, and change the name of the output file to avoid overwriting.
