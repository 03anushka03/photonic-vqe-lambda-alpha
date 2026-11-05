import perceval as pcvl
from perceval import catalog
import numpy as np
from scipy.optimize import minimize
from perceval.rendering import DisplayConfig, SymbSkin
theta0 = pcvl.Parameter("theta_0")
theta1 = pcvl.Parameter("theta_1")

def build_proc(t0, t1):
    """
    Exact circuit from the image:
    q[0]: ──────────────────── CRy(θ1) ──●──── measure
    q[1]: ── Ry(θ0) ──●──────────●───────⊕─── measure
    q[2]: ── X ────────⊕──────────────────────  measure

    Gates in order:
    1. Ry(θ0)  on q[1]
    2. X       on q[2]
    3. CNOT    control=q[1], target=q[2]
    4. CRy(θ1) control=q[1], target=q[0]
       decomposed as: Ry(+t1/2) · CNOT(q1->q0) · Ry(-t1/2) · CNOT(q1->q0)
    5. CNOT    control=q[0], target=q[1]

    Mode layout:
      q[0] -> modes 0,1
      q[1] -> modes 2,3
      q[2] -> modes 4,5
      herald pairs: 6,7 | 8,9 | 10,11 | 12,13
    """
    proc = pcvl.Processor("SLOS", 6)

    # 1. Ry(θ0) on q[1] (modes 2,3)
    proc.add([2, 3], catalog['ry'].build_circuit(theta=theta0))

    # 2. X on q[2] (modes 4,5)
    proc.add([4, 5], catalog['x'].build_circuit())

    # 3. CNOT: control=q[1], target=q[2]  (heralded cnot manages its own heralds)
    proc.add([2, 3, 4, 5], catalog['heralded cnot'].build_processor())

    # 4. CRy(θ1): control=q[1], target=q[0]
    #    Decomposition: Ry(+t1/2) · CNOT(q1->q0) · Ry(-t1/2) · CNOT(q1->q0)

    # 4a. Ry(+t1/2) on q[0] (modes 0,1)
    c_ry_pos = pcvl.Circuit(2)
    c_ry_pos.add(0, pcvl.BS.Ry(theta=t1/2))
    proc.add([0, 1], c_ry_pos)

    # 4b. CNOT: control=q[1], target=q[0]
    proc.add([2, 3, 0, 1], catalog['heralded cnot'].build_processor())

    # 4c. Ry(-t1/2) on q[0] (modes 0,1)
    c_ry_neg = pcvl.Circuit(2)
    c_ry_neg.add(0, pcvl.BS.Ry(theta=-t1/2))
    proc.add([0, 1], c_ry_neg)

    # 4d. CNOT: control=q[1], target=q[0]
    proc.add([2, 3, 0, 1], catalog['heralded cnot'].build_processor())

    # 5. CNOT: control=q[0], target=q[1]
    proc.add([0, 1, 2, 3], catalog['heralded cnot'].build_processor())

    return proc
def measure_expectation(theta, pauli_string, display1=False):
    if pauli_string == "III":
        return 1.0

    t0, t1 = float(theta[0]), float(theta[1])
    theta0.set_value(t0)
    theta1.set_value(t1)

    proc = build_proc(t0, t1)
    circuit=proc.linear_circuit()
    if(display1):
       #DisplayConfig.select_skin(SymbSkin)
       display(pcvl.pdisplay(circuit, recursive=False))
       return None
       
    for i, p in enumerate(pauli_string):
        if p == "X":
            proc.add([2*i, 2*i+1], catalog['ry'].build_circuit(theta=np.pi/2))
        elif p == "Y":
            proc.add([2*i, 2*i+1], catalog['rx'].build_circuit(theta=-np.pi/2))

    proc.with_input(pcvl.BasicState([1,0, 1,0, 1,0]))
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
        # Heralded cnot enforces post-selection internally, just check qubit modes
        if (s[0]+s[1] != 1) or (s[2]+s[3] != 1) or (s[4]+s[5] != 1):
            continue
        total_valid += prob
        bits = [s[5], s[3], s[1]]
        parity = 1
        for i, p in enumerate(pauli_string):
            if p != "I" and bits[i] == 1:
                parity *= -1
        exp += parity * prob

    if total_valid == 0:
        return 0
    return exp / total_valid
measure_expectation([0.1,0.1], "ZII", display1=True)
def energy(theta):
    Z0   = measure_expectation(theta, "IIZ")
    Z1   = measure_expectation(theta, "IZI")
    Z2   = measure_expectation(theta, "ZII")
    XX01 = measure_expectation(theta, "IXX")
    YY01 = measure_expectation(theta, "IYY")
    XX12 = measure_expectation(theta, "XXI")
    YY12 = measure_expectation(theta, "YYI")
    XZX  = measure_expectation(theta, "XZX")
    YZY  = measure_expectation(theta, "YZY")
    E = (
        3.89655
        + 0.75187  * Z0
        - 4.64842  * Z1
        - 0.629235 * (XX01 + YY01)
        + 9.75975  * (1 - Z2)
        + 2.27815  * (XX12 + YY12)
        - 2.53445  * (XZX + YZY)
    )
    return E
initial_theta = np.array([0.1, 0.1])
res = minimize(
    energy,
    initial_theta,
    method="COBYLA",
    tol=1e-6,
    options={"maxiter": 300}
)
print("Ground state energy:", res.fun)
print("Optimal thetas (deg):", np.degrees(res.x))
