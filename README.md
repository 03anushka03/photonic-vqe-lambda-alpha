# Photonic Variational Quantum Eigensolver for the Λ–α System

A variational quantum eigensolver (VQE) implementation for estimating the ground-state energy of a Λ–α system using photonic quantum circuits simulated with [Perceval](https://perceval.quandela.net/).

## Project overview

Hyperon–nucleon interactions are relevant to understanding the properties of dense nuclear matter, including matter inside neutron stars. The appearance of hyperons in neutron-star matter is connected to the so-called *hyperon puzzle*.

This project investigates the use of a variational quantum eigensolver to estimate the ground-state energy of a Λ–α system. The Hamiltonian is mapped to qubit operators and implemented using photonic quantum circuits in Perceval.

The project explores:

- Mapping the Λ–α Hamiltonian to a qubit representation.
- Constructing variational photonic circuits using Perceval.
- Estimating Hamiltonian expectation values through Pauli measurements.
- Optimising variational parameters using classical optimisation.
- Comparing VQE estimates with exact diagonalisation.
- Exploring larger circuit implementations and parallel evaluation of Hamiltonian terms.
- Testing a two-qubit implementation on Perceval's remote `sim:ascella` simulator.

**Important:** The remote experiment used a simulator, not a physical photonic quantum processor.

## Methodology

### 1. Hamiltonian construction

The Λ–α system is described using a model Hamiltonian constructed from the kinetic-energy and Λ–α interaction terms.

The matrix elements are evaluated in a chosen basis:

\[
h_{ij}=\langle i|T+U_{\Lambda\alpha}|j\rangle.
\]

The resulting Hamiltonian is mapped to qubit operators using the Jordan–Wigner transformation.

The resulting qubit Hamiltonian is expressed as a sum of Pauli strings:

\[
H=\sum_j c_jP_j,
\]

where \(c_j\) are the Hamiltonian coefficients and \(P_j\) are tensor products of Pauli operators.

### 2. Variational quantum eigensolver

The VQE method uses a parameterised quantum circuit to prepare a trial state:

\[
|\psi(\boldsymbol{\theta})\rangle
=U(\boldsymbol{\theta})|\psi_0\rangle.
\]

The energy expectation value is

\[
E(\boldsymbol{\theta})
=
\langle\psi(\boldsymbol{\theta})|
H
|\psi(\boldsymbol{\theta})\rangle
=
\sum_j c_j
\langle P_j\rangle_{\boldsymbol{\theta}}.
\]

A classical optimiser updates the circuit parameters to minimise the energy:

\[
E_0 \approx \min_{\boldsymbol{\theta}}E(\boldsymbol{\theta}).
\]

The optimisation is performed using SciPy's COBYLA implementation.

### 3. Photonic circuit implementation

The variational circuits are constructed using Perceval. Dual-rail photonic encoding is used, with each logical qubit represented by two optical modes.

The implementation uses:

- Perceval circuit and processor objects.
- Photonic \(R_y\) and \(X\) operations.
- Heralded CNOT gates.
- Measurement-basis rotations for Pauli expectation values.
- Postselection on the required photon-detection patterns.

For each Pauli string, the circuit is configured with the appropriate measurement-basis rotations. The output distribution is then used to calculate the corresponding expectation value.

For example, the two-qubit Hamiltonian is

\[
\begin{aligned}
H_2={}&3.89655I+0.75187Z_0-4.64842Z_1\\
&-0.629235(X_0X_1+Y_0Y_1).
\end{aligned}
\]

Its energy is evaluated by measuring the four required Pauli terms:

\[
\langle Z_0\rangle,\quad
\langle Z_1\rangle,\quad
\langle X_0X_1\rangle,\quad
\langle Y_0Y_1\rangle.
\]

### 4. Simulation and execution

The project contains multiple implementations:

| Implementation | Execution method |
|---|---|
| Two-qubit VQE | Local Perceval SLOS simulation |
| Three-qubit VQE | Local Perceval SLOS simulation |
| Two-qubit remote experiment | Perceval remote `sim:ascella` simulator |
| Four-qubit VQE | Local simulation with multiprocessing |

The local implementations use Perceval's SLOS backend. The remote two-qubit experiment uses a finite-shot sampling workflow.

The four-qubit implementation evaluates the Hamiltonian's Pauli terms concurrently using Python multiprocessing. This parallelises the measurement simulations, rather than the quantum circuit itself.

## Results

### Local simulation benchmarks

The local VQE results were compared with independently calculated exact eigenvalues of the corresponding qubit Hamiltonians.

| Qubit implementation | VQE energy | Exact ground-state energy | Absolute difference |
|---|---:|---:|---:|
| 2 qubits | -1.64843682 | -1.64843682 | Approximately \(0\) |
| 3 qubits | -2.66207725 | -2.66615292 | 0.00407567 |

The two-qubit VQE result agrees with exact diagonalisation to numerical precision.

The three-qubit VQE result is slightly higher than the exact ground-state energy, consistent with the variational upper-bound property of VQE.

### Remote two-qubit experiment

A two-qubit implementation was also submitted to Perceval's remote `sim:ascella` simulator using finite-shot sampling.

The experiment used a sampling workflow with postselection on valid photonic outcomes. The reported diagnostic run at \(\theta=0.2289515\) produced the following partial measurement output:

| Measurement | Reported expectation value |
|---|---:|
| \( \langle Z_1\rangle \) (`ZI`) | 0.99410029 |
| \( \langle X_1X_0\rangle \) (`XX`) | -0.07917889 |
| \( \langle Y_1Y_0\rangle \) (`YY`) | 0.07647059 |

The reported energy from this diagnostic run was \(-0.24053920\).

This value is **not presented as a validated ground-state estimate**. The available diagnostic output is incomplete, and the remote result has not yet been reconciled with the local exact-diagonalisation benchmark. The remote experiment is therefore included as a preliminary execution test.

### Four-qubit implementation

A four-qubit circuit and multiprocessing-based measurement workflow were developed to explore a larger Hamiltonian.

This implementation is experimental. Its controlled-rotation sequence and optimisation convergence require validation before a numerical result can be reported as a benchmark.

## Repository structure

```text
photonic-vqe-lambda-alpha/
├── README.md
├── requirements.txt
├── .gitignore
├── LICENSE
├── src/
│   ├── vqe_2qubit.py
│   ├── vqe_3qubit.py
│   ├── remote_2qubit.py
│   └── experimental/
│       └── vqe_4qubit_parallel.py
├── benchmarks/
│   ├── exact_diagonalization.py
│   └── benchmark_results.csv
├── results/
├── figures/
├── poster/
│   └── project_poster.pdf
└── references/
    └── references.md
```

## Installation

Create and activate a Python environment, then install the dependencies.

```bash
python -m venv .venv
```

On Linux/macOS:

```bash
source .venv/bin/activate
```

On Windows:

```powershell
.venv\Scripts\activate
```

Install the dependencies:

```bash
pip install -r requirements.txt
```

The exact package versions should be recorded from the environment in which the scripts were tested.

## Running the implementations

Run the local two-qubit implementation:

```bash
python src/vqe_2qubit.py
```

Run the local three-qubit implementation:

```bash
python src/vqe_3qubit.py
```

The four-qubit multiprocessing implementation is experimental:

```bash
python src/experimental/vqe_4qubit_parallel.py
```

The remote two-qubit script requires Perceval remote-simulator access and a valid authentication token. It should not be run without configuring the token as an environment variable.

## Limitations and future work

- Validate the remote two-qubit sampling results against the local simulation and exact diagonalisation.
- Investigate finite-shot statistical uncertainty and postselection success probability.
- Validate the measurement-basis rotation conventions independently.
- Correct and verify the four-qubit controlled-rotation sequence.
- Study optimisation convergence for the larger circuit.
- Compare serial and multiprocessing execution times.
- Investigate the computational and sampling overhead as the number of qubits increases.
- Explore execution on physical photonic quantum hardware, subject to availability.

## Acknowledgements

This project was carried out as an academic research project at IIT Mandi.

The implementation uses [Perceval](https://perceval.quandela.net/) for photonic quantum circuit construction and simulation, and SciPy for classical optimisation.

## References

Add the original nuclear-model reference, the source of the Λ–α interaction parameters, the Jordan–Wigner mapping reference, the Perceval documentation, and relevant VQE literature here.
