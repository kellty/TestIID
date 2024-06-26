# Implementation of Tests for Data Being IID

## Data
The MNIST dataset from [http://yann.lecun.com/exdb/mnist/](http://yann.lecun.com/exdb/mnist/), two network datasets from [https://github.com/matteoiacopini/ZIL-T-MS](https://github.com/matteoiacopini/ZIL-T-MS) and a dataset of air pollutants from [https://doi.org/10.24432/C5RK5G](https://doi.org/10.24432/C5RK5G) were used for this work. The files "***mnist.RData***", "***Data_email.mat***", "***Data_financial.mat***" and "***PRSA_Data_Aotizhongxin_20130301-20170228.csv***" are corresponding copies.

## Code
- "***func.R***": Build functions to implement the IID tests.
- "***simu.R***": Perform IID tests on simulated data in various settings.
- "***real.R***": Perform IID tests on the aforementioned real data (with certain preprocessing for the MNIST dataset) and their bootstraped samples.

## Workflow
- Simulation: Make sure that "***func.R***" is in the working directory, and *then* run "***simu.R***".
- Real data analysis: Make sure that "***func.R***", "***mnist.RData***", "***Data_email.mat***" and "***Data_financial.mat***" are in the working directory, and *then* run "***real.R***".
