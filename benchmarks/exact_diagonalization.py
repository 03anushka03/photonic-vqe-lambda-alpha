import numpy as np


I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def pauli_string(operators):
    """
    Construct a tensor product of Pauli operators.

    operators are listed in the order:
    q_(n-1), ..., q_1, q_0
    """
    result = np.array([[1]], dtype=complex)

    for operator in operators:
        result = np.kron(result, operator)

    return result


def H2():
    identity = pauli_string([I, I])

    return (
        3.89655 * identity
        + 0.75187 * pauli_string([I, Z])
        - 4.64842 * pauli_string([Z, I])
        - 0.629235 * (
            pauli_string([X, X])
            + pauli_string([Y, Y])
        )
    )


def H3():
    identity = pauli_string([I, I, I])

    return (
        3.89655 * identity
        + 0.75187 * pauli_string([I, I, Z])
        - 4.64842 * pauli_string([I, Z, I])
        - 0.629235 * (
            pauli_string([I, X, X])
            + pauli_string([I, Y, Y])
        )
        + 9.75975 * (
            identity - pauli_string([Z, I, I])
        )
        + 2.27815 * (
            pauli_string([I, X, X])
            + pauli_string([I, Y, Y])
        )
        - 2.53445 * (
            pauli_string([X, Z, X])
            + pauli_string([Y, Z, Y])
        )
    )


def H4():
    identity = pauli_string([I, I, I, I])

    return (
        3.870065 * identity
        + 0.795235 * pauli_string([I, I, I, Z])
        - 4.6653 * pauli_string([I, I, Z, I])
        - 0.626698 * (
            pauli_string([I, I, X, X])
            + pauli_string([I, I, Y, Y])
        )
        + 9.79092 * (
            identity - pauli_string([I, Z, I, I])
        )
        + 2.30925 * (
            pauli_string([I, I, X, X])
            + pauli_string([I, I, Y, Y])
        )
        - 2.49918 * (
            pauli_string([X, Z, I, X])
            + pauli_string([Y, Z, I, Y])
        )
        + 14.6668 * (
            identity - pauli_string([Z, I, I, I])
        )
        + 4.99733 * (
            pauli_string([Z, X, X, I])
            + pauli_string([Z, Y, Y, I])
        )
        - 1.77978 * (
            pauli_string([X, Z, Z, X])
            + pauli_string([Y, Z, Z, Y])
        )
        - 2.19287 * (
            pauli_string([I, X, Z, X])
            + pauli_string([I, Y, Z, Y])
        )
    )


def calculate_eigenvalues(H, name):
    eigenvalues = np.linalg.eigvalsh(H)

    print(f"\n{name}")
    print("-" * 40)

    for i, energy in enumerate(eigenvalues):
        print(f"Eigenvalue {i}: {energy:.10f}")

    print(f"\nGround-state energy: {eigenvalues[0]:.10f}")

    return eigenvalues


if __name__ == "__main__":
    calculate_eigenvalues(H2(), "Two-qubit Hamiltonian")
    calculate_eigenvalues(H3(), "Three-qubit Hamiltonian")
    calculate_eigenvalues(H4(), "Four-qubit Hamiltonian")
