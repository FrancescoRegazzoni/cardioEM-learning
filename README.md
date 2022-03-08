# A Machine Learning method for real-time numerical simulations of cardiac electromechanics

This repository contains the code accompanying the paper [1], in which we propose a machine learning method to build a system of differential equations that approximates the dynamics of 3D electromechanical models for the human heart, accounting for the dependence on a set of parameters.

The code is based on the library [**<kbd>model-learning</kbd>**](https://github.com/FrancescoRegazzoni/model-learning) and the methods proposed in [2].

## Installation

1. Clone the library [**<kbd>model-learning</kbd>**](https://github.com/FrancescoRegazzoni/model-learning). Then, follow the [instructions](https://model-learning.readthedocs.io/en/latest/installation.html) therein contained to set-up the library. If you use `ssh` (you can also use `https` or simply download the files), execute:
```bash
git clone git@github.com:FrancescoRegazzoni/model-learning.git
```

2. Move into the `examples` folder.
```bash
cd model-learning/examples
```

3. Clone this repository under the name `app_cardioEM-learning`.
```bash
git clone git@github.com:FrancescoRegazzoni/cardioEM-learning.git app_cardioEM-learning
```

## Documentation

The code is divided into four steps. The corresponding scripts are identified with the prefix `"step_n_..."`.

In this repository we consider two use cases. In the first one, we consider variability with respect to a single parameter of the electromechanical model, while in the second one we vary four parameters simultaneously (see [1] for more details). In the scripts, these use cases are respectively identified with labels:
```matlab
use_case = 'one_param';
use_case = 'all_params';
```

### Step 1: Dataset generation

First, you need to generate the training data and convert them into the format of the library[**<kbd>model-learning</kbd>**](https://github.com/FrancescoRegazzoni/model-learning). For demonstration purposes, in the folder `data` of this repository you can find the results of simulations that we have generated using [**<kbd>life<sup>x</sup></kbd>**](https://lifex.gitlab.io/). For more details, see [1]. We recall that, thanks to the non-intrusive nature of the method proposed in [1], you can use data generated with any cardiac mechanics model, provided that you are able to export the pressure and volume transients.

### Step 2: Dataset visualization

Now that you have converted the data into the correct format, it is a good idea to look at what they look like.

### Step 3: Training

You can now train a reduced-order model (ROM). To do this, simply use the command:
```matlab
model_learn('path/to/setting/file.ini')
```
You can find two pre-compiled setting files in the folder `training_options`.
Training of the models will take a few hours. 

### Step 4: Using the trained models

The [**<kbd>model-learning</kbd>**](https://github.com/FrancescoRegazzoni/model-learning) library saves the trained models in a standard path. It is possible to load the models later and to deploy them, e.g. in combination with other applications. Below we show how the models can be loaded in either **<kbd>MATLAB</kbd>** or **<kbd>Python</kbd>**. If you are interested in loading the models in **<kbd>C++</kbd>** or other languages, feel free to contact us.

**Loading the model in Matlab**. To open the trained model, use for example the commands
```matlab
problem = problem_get('app_cardioEM-learning', 'problems/EM_one_param.ini');
ANNmod = read_model_fromfile(problem, model_dir);
```
where `model_dir` is the model folder (it is printed on the screen during Step 3, make a note of that!). Now, the right-hand side and the initial state of the trained models are respectively available as `ANNmod.f` and `ANNmod.x0`.

**Loading the model in Python**. To load the model in **<kbd>Python</kbd>**, you can use the wrapper [`pyModelLearning`](https://model-learning.readthedocs.io/en/latest/installation.html#python-wrapper).
```python
from pyModelLearning import ANNmodel
ANNmod = ANNmodel('app_cardioEM-learning/networks_EM_all_params/' + model_dir + '/compact.mat')
```
Now, the right-hand side and the initial state of the trained models are respectively available as `ANNmod.rhs` and `ANNmod.initial_state`.

In this repository we show an example of the use of ROM, implemented both in **<kbd>MATLAB</kbd>** and in **<kbd>Python</kbd>**. To demonstrate the versatility of our method, we test the ROM by coupling it with a circulation model different than the one used to generate the training dataset. In particular, we use the windkessel model described in [this paper](https://doi.org/10.1016/j.cma.2020.113268). By running the script, plots comparing the PV loops obtained with the full-order model (FOM) and the reduced-order model (ROM) will be shown, and errors will be printed.

## References

[1] F. Regazzoni, M. Salvador, L. Dede', A. Quarteroni [A Machine Learning method for real-time numerical simulations of cardiac electromechanics](https://arxiv.org/abs/2110.13212), *Computer Methods in Applied Mechanics and Engineering* (2022).

[2] F. Regazzoni, L. Dede', A. Quarteroni [Machine learning for fast and reliable solution of time-dependent differential equations](https://doi.org/10.1016/j.jcp.2019.07.050), *Journal of Computational Physics* (2019).

