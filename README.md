# CoaTran: Coalescent tree simulation along a transmission network
CoaTran is a tool that, given a transmission network and sample times, will simulate a coalescent phylogeny constrained by the transmission network. CoaTran is similar in functionality to the [VirusTreeSimulator](https://github.com/PangeaHIV/VirusTreeSimulator) component of [PANGEA.HIV.sim](https://github.com/olli0601/PANGEA.HIV.sim), but CoaTran is consistently ~100x faster.

# Installation
To compile CoaTran, simply download the latest release or clone this repository, and then compile using `make`:

```bash
$ make
g++ -Wall -pedantic -std=c++11 -O3 -o coatran_constant main.cpp common.cpp common.h coalescent.cpp coalescent.h
g++ -Wall -pedantic -std=c++11 -O3 -DEXPGROWTH -o coatran_expgrowth main.cpp common.cpp common.h coalescent.cpp coalescent.h
g++ -Wall -pedantic -std=c++11 -O3 -DTRANSTREE -o coatran_transtree main.cpp common.cpp common.h coalescent.cpp coalescent.h
g++ -Wall -pedantic -std=c++11 -O3 -DINFTIME -o coatran_inftime main.cpp common.cpp common.h coalescent.cpp coalescent.h
```

If you want to debug/benchmark, you can compile the debug executables using `make debug`:

```bash
$ make debug
g++ -Wall -pedantic -std=c++11 -O0 -g -pg -o coatran_constant_debug main.cpp common.cpp common.h coalescent.cpp coalescent.h
g++ -Wall -pedantic -std=c++11 -O0 -g -pg -DEXPGROWTH -o coatran_expgrowth_debug main.cpp common.cpp common.h coalescent.cpp coalescent.h
g++ -Wall -pedantic -std=c++11 -O0 -g -pg -DTRANSTREE -o coatran_transtree_debug main.cpp common.cpp common.h coalescent.cpp coalescent.h
g++ -Wall -pedantic -std=c++11 -O0 -g -pg -DINFTIME -o coatran_inftime_debug main.cpp common.cpp common.h coalescent.cpp coalescent.h
```

# Usage
When compiled, CoaTran will produce different executables depending on the model of effective population size you choose to use. All modes have at least the following two parameters:

* **`<trans_network>`:** The transmission network, in the [FAVITES format](https://github.com/niemasd/FAVITES/wiki/File-Formats#transmission-network-file-format)
* **`<sample_times>`:** The sample times, in the [FAVITES format](https://github.com/niemasd/FAVITES/wiki/File-Formats#sample-time-file-format)

In all modes, you can specify a constant random number generator seed (e.g. for reproducibility) by setting the `COATRAN_RNG_SEED` environment variable:

```bash
export COATRAN_RNG_SEED=42
```

The Newick trees output by CoaTran have unifurcations (i.e., an internal node with a single child) at the times of infection, which may be useful information. However, if you want to suppress unifurcations (i.e., merge the branches above and below the unifurcating node), you can do so easily with tools like [TreeSwift](https://github.com/niemasd/TreeSwift) or [DendroPy](https://dendropy.org/):

```python3
from treeswift import read_tree_newick
coatran_output_file = "my_tree.nwk"
tree = read_tree_newick(coatran_output_file)
tree.suppress_unifurcations()
print(tree.newick())
```

## Constant Effective Population Size (`coatran_constant`)
You can use `coatran_constant` to simulate phylogenies under coalescence with constant effective population size:

```bash
coatran_constant <trans_network> <sample_times> <eff_pop_size>
```

* **`<eff_pop_size>`:** The effective population size, which remains constant

## ~Exponential Effective Population Size Growth~
**THIS MODE DOES NOT WORK YET!!!**

~You can use `coatran_expgrowth` to simulate phylogenies under coalescence with exponential effective population size growth from the time of infection:~

```bash
coatran_expgrowth <trans_network> <sample_times> <init_eff_pop_size> <eff_pop_growth>
```

* **`<init_eff_pop_size>`:** The initial effective population size at the time of infection (N0)
* **`<eff_pop_growth>`:** The growth rate of the effective population size

## Transmission Tree
You can use `coatran_transtree` to simulate phylogenies that are equivalent to the transmission tree. In other words, if *u* infected *v*, coalescence of their lineages happens as late in time as possible: the time at which *u* infected *v*.

```bash
coatran_transtree <trans_network> <sample_times>
```

## Infection Time
You can use `coatran_inftime` to simulate phylogenies such that coalescence happens at the time of infection. In other words, if *u* infected *v*, coalescence of their lineages happens as early in time as possible: the time at which *u* was infected.

```bash
coatran_inftime <trans_network> <sample_times>
```

# Citing CoaTran
If you use CoaTran in your work, please cite:

> Moshiri N (2020). "CoaTran: Coalescent tree simulation along a transmission network." *bioRxiv*. [doi:10.1101/2020.11.10.377499](https://doi.org/10.1101/2020.11.10.377499)
