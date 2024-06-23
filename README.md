# Implementation of Tests for Data Being IID

## Data
The MNIST dataset from [http://yann.lecun.com/exdb/mnist/](http://yann.lecun.com/exdb/mnist/) and two network datasets from [https://github.com/matteoiacopini/ZIL-T-MS](https://github.com/matteoiacopini/ZIL-T-MS) were used for this work. The files "***mnist.RData***", "***Data_email.mat***" and "***Data_financial.mat***" are corresponding copies.

## Code
- "***func.R***": Build functions to implement the IID tests.
- "***simu.R***": Perform IID tests on simulated data in various settings.
- "***real.R***": Perform IID tests on the aforementioned real data (with certain preprocessing for the MNIST dataset) and their bootstraped samples.

## Workflow
- Simulation: Make sure that "***func.R***" is in the working directory, and *then* run "***simu.R***".
- Real data analysis: Make sure that "***func.R***", "***mnist.RData***", "***Data_email.mat***" and "***Data_financial.mat***" are in the working directory, and *then* run "***real.R***".
