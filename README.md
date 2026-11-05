# Photonic Variational Quantum Eigensolver for the Λ–α System

A Variational Quantum Eigensolver (VQE) implementation for estimating the ground-state energy of a Λ–α system using photonic quantum circuits and the Perceval framework.

## 1. Project Overview

Hyperon–nucleon interactions are important for understanding dense nuclear matter, including matter inside neutron stars. The appearance of hyperons in neutron-star matter is connected to the so-called hyperon puzzle.

This project investigates the application of VQE to estimate the ground-state energy of a Λ–α system.

The Hamiltonian is mapped to a qubit representation and implemented using photonic quantum circuits in Perceval.

The project explores:

- Mapping the Λ–α Hamiltonian to a qubit representation.
- Constructing variational photonic circuits using Perceval.
- Estimating Hamiltonian expectation values through Pauli measurements.
- Optimising variational parameters using classical optimisation.
- Comparing VQE estimates with exact diagonalisation.
- Exploring larger circuit implementations and parallel evaluation of Hamiltonian terms.
- Testing a two-qubit implementation on Perceval's remote sim:ascella simulator.

**Note:** The remote experiment used a simulator, not a physical photonic quantum processor.

## 2. Methodology

### 2.1 Hamiltonian construction

The Λ–α system is described using a model Hamiltonian:

H = T + U_Λα

The Hamiltonian matrix elements in the chosen basis are:

h_ij = <i|T + U_Λα|j>

where T is the kinetic-energy operator and U_Λα is the Λ–α interaction potential.

The resulting Hamiltonian is mapped to qubit operators using the Jordan–Wigner transformation.

The qubit Hamiltonian is expressed as:

H = Σ_j c_j P_j

where c_j are Hamiltonian coefficients and P_j are tensor products of Pauli operators.

### 2.2 Variational Quantum Eigensolver

A parameterised quantum circuit prepares a trial state:

|ψ(θ)> = U(θ)|ψ₀>

The energy expectation value is:

E(θ) = <ψ(θ)|H|ψ(θ)>

Using the Pauli decomposition:

E(θ) = Σ_j c_j <P_j>_θ

The classical optimiser updates the variational parameters to minimise the energy:

E₀ ≈ min_θ E(θ)

The optimisation is performed using the COBYLA algorithm from SciPy.

### 2.3 Photonic circuit implementation

The variational circuits are constructed using Perceval.

Dual-rail photonic encoding is used, where each logical qubit is represented by two optical modes.

The implementation uses:

- Perceval circuit and processor objects.
- Photonic X and Ry operations.
- Heralded CNOT gates.
- Measurement-basis rotations for Pauli expectation values.
- Postselection on required photon-detection patterns.

For each Pauli string, the circuit is configured with appropriate measurement-basis rotations. The output distribution is used to calculate the corresponding expectation value.

### 2.4 Two-qubit Hamiltonian

H₂ = 3.89655 I + 0.75187 Z₀ - 4.64842 Z₁
     - 0.629235 (X₀X₁ + Y₀Y₁)

The energy is evaluated using:

- <Z₀>
- <Z₁>
- <X₀X₁>
- <Y₀Y₁>

### 2.5 Three-qubit Hamiltonian

H₃ = H₂ + 9.75975(I - Z₂)
     + 2.27815(X₁X₂ + Y₁Y₂)
     - 2.53445(X₀Z₁X₂ + Y₀Z₁Y₂)

### 2.6 Four-qubit Hamiltonian

H₄ = 3.870065 I + 0.795235 Z₀ - 4.6653 Z₁
     - 0.626698(X₀X₁ + Y₀Y₁)
     + 9.79092(I - Z₂)
     + 2.30925(X₁X₂ + Y₁Y₂)
     - 2.49918(X₀Z₁X₂ + Y₀Z₁Y₂)
     + 14.6668(I - Z₃)
     + 4.99733(X₂X₃ + Y₂Y₃)
     - 1.77978(X₀Z₁Z₂X₃ + Y₀Z₁Z₂Y₃)
     - 2.19287(X₁Z₂X₃ + Y₁Z₂Y₃)

## 3. Simulation and Execution

| Implementation | Execution method |
|---|---|
| Two-qubit VQE | Local Perceval SLOS simulation |
| Three-qubit VQE | Local Perceval SLOS simulation |
| Two-qubit remote experiment | Perceval remote sim:ascella simulator |
| Four-qubit VQE | Local simulation with Python multiprocessing |

### Local simulation

The two- and three-qubit implementations use Perceval's SLOS backend.

### Remote simulation

The two-qubit implementation was also tested using Perceval's remote sim:ascella simulator.

The remote execution uses finite-shot sampling and postselection on valid photonic outcomes.

The available remote diagnostic run evaluates the circuit at a fixed variational parameter. It should be treated as a preliminary remote execution experiment rather than a completed remote VQE optimisation.

### Parallel implementation

The four-qubit implementation uses Python multiprocessing to distribute the evaluation of independent Pauli terms across CPU processes.

This parallelises the measurement simulations, rather than the quantum circuit or classical optimiser itself.

## 4. Results

### 4.1 Local simulation benchmarks

The local VQE results were compared with independently calculated exact eigenvalues of the corresponding qubit Hamiltonians.

| Qubits | VQE energy | Exact ground-state energy | Absolute difference |
|---|---:|---:|---:|
| 2 | -1.64843682 | -1.64843682 | Approximately 0 |
| 3 | -2.66207725 | -2.66615292 | 0.00407567 |

The two-qubit VQE result agrees with exact diagonalisation to numerical precision.

The three-qubit VQE result is slightly higher than the exact ground-state energy, consistent with the variational upper-bound property of VQE.

### 4.2 Remote two-qubit experiment

A diagnostic run at:

θ = 0.2289515

produced the following partial measurement output:

| Measurement | Expectation value |
|---|---:|
| ZI | 0.99410029 |
| XX | -0.07917889 |
| YY | 0.07647059 |

The reported diagnostic energy was:

E = -0.24053920

This value is not presented as a validated ground-state estimate. The available diagnostic output is incomplete, and the remote result has not yet been reconciled with the local exact-diagonalisation benchmark.

The remote experiment is therefore included as a preliminary execution test.

### 4.3 Four-qubit implementation

A four-qubit circuit and multiprocessing-based measurement workflow were developed to explore a larger Hamiltonian.

This implementation remains experimental. The controlled-rotation sequence and optimisation convergence require further validation before a numerical result can be reported as a benchmark.

## 5. Repository Structure

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
