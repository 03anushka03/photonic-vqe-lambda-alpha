import perceval as pcvl
from perceval import catalog
import numpy as np
from scipy.optimize import minimize
from multiprocessing import Pool
import os

# ---------------------------------------------------------------------------
# How many CPU cores to use.
# On your HPC cluster, set this to however many cores your job requests.
# e.g. if your SLURM script has --cpus-per-task=16, set N_WORKERS=16
# ---------------------------------------------------------------------------
N_WORKERS = int(os.environ.get("VQE_WORKERS", os.cpu_count() or 4))
print(f"[VQE] Using {N_WORKERS} parallel workers")

# ---------------------------------------------------------------------------
# Circuit builder  (unchanged from your original)
# ---------------------------------------------------------------------------

def build_proc(t0, t1, t2):
    proc = pcvl.Processor("Naive", 8)

    proc.add([4, 5], catalog['ry'].build_circuit(theta=t0))
    proc.add([6, 7], catalog['x'].build_circuit())
    proc.add([4, 5, 6, 7], catalog['heralded cnot'].build_processor())

    c_ry_pos1 = pcvl.Circuit(2)
    c_ry_pos1.add(0, pcvl.BS.Ry(theta=t1 / 2))
    proc.add([2, 3], c_ry_pos1)

    proc.add([4, 5, 2, 3], catalog['heralded cnot'].build_processor())

    c_ry_neg1 = pcvl.Circuit(2)
    c_ry_neg1.add(0, pcvl.BS.Ry(theta=-t1 / 2))
    proc.add([2, 3], c_ry_neg1)

    proc.add([4, 5, 2, 3], catalog['heralded cnot'].build_processor())
    proc.add([2, 3, 4, 5], catalog['heralded cnot'].build_processor())

    c_ry_pos2 = pcvl.Circuit(2)
    c_ry_pos2.add(0, pcvl.BS.Ry(theta=t2 / 2))
    proc.add([0, 1], c_ry_pos2)

    proc.add([2, 3, 0, 1], catalog['heralded cnot'].build_processor())

    c_ry_neg2 = pcvl.Circuit(2)
    c_ry_neg2.add(0, pcvl.BS.Ry(theta=-t1 / 2))
    proc.add([0, 1], c_ry_neg2)

    proc.add([2, 3, 0, 1], catalog['heralded cnot'].build_processor())
    proc.add([0, 1, 2, 3], catalog['heralded cnot'].build_processor())

    return proc


# ---------------------------------------------------------------------------
# Single measurement  (now a plain top-level function so it can be pickled
# and sent to worker processes)
# ---------------------------------------------------------------------------

def measure_expectation(args):
    """
    args = (theta, pauli_string)

    Kept as a single-argument function so multiprocessing.Pool.map() works
    cleanly (it only passes one argument per call).
    """
    theta, pauli_string = args

    if pauli_string == "IIII":
        return 1.0

    t0, t1, t2 = float(theta[0]), float(theta[1]), float(theta[2])

    # Each worker builds its own processor — no shared state, no race conditions
    proc = build_proc(t0, t1, t2)

    for i, p in enumerate(pauli_string):
        if p == "X":
            proc.add([2 * i, 2 * i + 1], catalog['ry'].build_circuit(theta=np.pi / 2))
        elif p == "Y":
            proc.add([2 * i, 2 * i + 1], catalog['rx'].build_circuit(theta=-np.pi / 2))

    proc.with_input(pcvl.BasicState([1, 0, 1, 0, 1, 0, 1, 0]))
    raw_res = proc.probs()

    if "results" in raw_res:
        dist = raw_res["results"]
    elif "probs" in raw_res:
        dist = raw_res["probs"]
    else:
        dist = raw_res

    total_valid = 0
    exp = 0
    for state, prob in dist.items():
        s = list(state)
        if (s[0] + s[1] != 1) or (s[2] + s[3] != 1) or \
           (s[4] + s[5] != 1) or (s[6] + s[7] != 1):
            continue
        total_valid += prob
        bits = [s[7], s[5], s[3], s[1]]
        parity = 1
        for i, p in enumerate(pauli_string):
            if p != "I" and bits[i] == 1:
                parity *= -1
        exp += parity * prob

    return 0.0 if total_valid == 0 else exp / total_valid


# ---------------------------------------------------------------------------
# Energy function  — all 16 measurements fired in parallel
# ---------------------------------------------------------------------------

PAULI_TERMS = [
    "IIIZ", "IIZI", "IZII", "ZIII",
    "IIXX", "IIYY",
    "IXXI", "IYYI",
    "XXII", "YYII",
    "IXZX", "IYZY",
    "XZXI", "YZYI",
    "XZZX", "YZZY",
]

def energy(theta):
    # Build the argument list once — each element is (theta, pauli_string)
    args = [(theta, p) for p in PAULI_TERMS]

    # Fire all 16 simulations in parallel across N_WORKERS processes
    with Pool(processes=N_WORKERS) as pool:
        results = pool.map(measure_expectation, args)

    # Unpack in the same order as PAULI_TERMS
    (Z0, Z1, Z2, Z3,
     XX01, YY01,
     XX12, YY12,
     XX23, YY23,
     XZX_012, YZY_012,
     XZX_123, YZY_123,
     XZZX_0123, YZZY_0123) = results

    E = (
        3.870065
        + 0.795235  * Z0
        - 4.6653    * Z1
        - 0.626698  * (XX01  + YY01)
        + 9.79092   * (1 - Z2)
        + 2.30925   * (XX12  + YY12)
        - 2.49918   * (XZX_012 + YZY_012)
        + 14.6668   * (1 - Z3)
        + 4.99733   * (XX23  + YY23)
        - 1.77978   * (XZZX_0123 + YZZY_0123)
        - 2.19287   * (XZX_123 + YZY_123)
    )
    return E


# ---------------------------------------------------------------------------
# Main  (guard required for multiprocessing on Linux/HPC)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    initial_theta = np.array([0.1, 0.1, 0.1])

    res = minimize(
        energy,
        initial_theta,
        method="COBYLA",
        tol=1e-6,
        options={"maxiter": 20},
    )

    print("Ground state energy:", res.fun)
    print("Optimal thetas (deg):", np.degrees(res.x))
