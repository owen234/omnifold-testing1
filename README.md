## Sandbox for exploring OmniFold

Various tests.

Some things were taken or adapted from Ben Nachman's OmniFold tutorial.  https://github.com/hep-lbdl/OmniFold

Here's a few descriptions of things that are in here.



Most important stuff

- **omnifold6b.py** : Code for Omnifold.  Adapted from omnifold.py in https://github.com/hep-lbdl/OmniFold
  
- **toy_study6b-2d-7b.ipynb** : Complete example of a 2D multi-variate Gaussian study using a NN in Omnifold.  Uses omnifold6b.py.

- **toy_study6b-4d-pdcov.ipynb** : Complete example of a 4D multi-variate Gaussian study using a NN in Omnifold.  Uses omnifold6b.py.

- **bootstrap_toy_study6b-2d.ipynb** : Runs bootstraps for a 2D multi-variate Gaussian using a NN in Omnifold.  Uses omnifold6b.py.

- **analyze-bootstrap6b.ipynb** :  Analyzes the event weight correlations of the bootstrap study.  Saves the results for plotting by plot-bootstrap6b.ipynb.

- **plot-bootstrap6b.ipynb** :  Make plots of event weight correlations as a function of distance in feature space.  Uses output of analyze-bootstrap6b.ipynb.

- **RooMultiVarGaussian2e.cxx,h** :  RooFit class for a multivariate Gaussian PDF.  Used for fitting bootstrap samples.  Comple the cxx file before running the Jupyter Notebook by starting root and then typing `.L RooMultiVarGaussian2e.cxx++`

- **bootstrap-nd-fitting6.ipynb** :  Performs unbinned maximum likelihood fits of bootstrap samples.  Uses RooMultiVarGaussian2e.cxx.

- **inspect-bootstraps.ipynb** :  Inspects the bootstrap samples produced by bootstrap_toy_study*.ipynb.  Does simple calculations of model parameters from the sample.  Compares these to the results of unbinned maximum likelihood fits of the samples from bootstrap-nd-fitting6.ipynb, if they are available.




Less important stuff

- **multivariate_normal.ipynb** :  Simple code to generate random numbers from a multivariate normal distribution.

- **roofit-testing-rf611.ipynb** :  Jupyter notebook adaptation of $ROOTSYS/tutorials/roofit/rf611_weightedfits.C

- **roofit-testing-rf101.ipynb** :  Jupyter notebook adaptation of $ROOTSYS/tutorials/roofit/rf611_basics.C

- **jnb-mvn2d-test1.ipynb** : Runs a test of the RooMvn2d.cxx code which is a class for a 2D multivariate normal PDF in RooFit.  You need to compile the cxx file before running the Jupyter Notebook by starting root and then typing `.L RooMvn2d.cxx++`



- 
