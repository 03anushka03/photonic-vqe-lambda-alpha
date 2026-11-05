import perceval as pcvl
from perceval import catalog
from perceval.algorithm import Sampler
import numpy as np
import time
from scipy.optimize import minimize
TOKEN = "_T_eyJhbGciOiJIUzUxMiIsInR5cCI6IkpXVCJ9.eyJpZCI6MjYxMiwiZXhwIjoxNzc2NDkzMzgwLjExNDE4MDN9.Xj_cRpk1lwcBian4tfYmoywofx9kLVYbRezb7Dc5fv4DQYrfxrEVx8OLSMB6rlooPWPjQJXrBuudcFr8dvpk3w"  
BACKEND = "sim:ascella"

def wait_for_job(job, pauli_string, poll_interval=2):
    while not job.is_complete:
        progress = job.status.progress
        print(f"  [{pauli_string}] progress: {progress:.0f}%", end="\r")
        time.sleep(poll_interval)
    print(f"  [{pauli_string}] done")
def measure_expectation(t, pauli_string, shots=20000):
    if pauli_string == "II":
        return 1.0

    proc = pcvl.Processor("SLOS", 4)
    proc.add([0, 1], catalog['x'].build_circuit())
    proc.add([2, 3], catalog['ry'].build_circuit(theta=float(t)))
    proc.add([2, 3, 0, 1], catalog['heralded cnot'].build_processor())

    for i, p in enumerate(pauli_string):
        if p == "X":
            proc.add([2*i, 2*i+1], catalog['ry'].build_circuit(theta=np.pi/2))
        elif p == "Y":
            proc.add([2*i, 2*i+1], catalog['rx'].build_circuit(theta=-np.pi/2))

    lc = proc.linear_circuit()
    print("Circuit modes: " + str(lc.m)) 

    remote_proc = pcvl.RemoteProcessor(BACKEND, token=TOKEN, m=lc.m)
    remote_proc.set_circuit(lc)

    remote_proc.with_input(pcvl.BasicState([1, 0, 1, 0, 0, 0]))
    remote_proc.min_detected_photons_filter(2)

    remote_proc.set_parameter("final_mode_number", 4)

    sampler = Sampler(remote_proc, max_shots_per_call=shots)
    job = sampler.sample_count.execute_async(shots)
    print("Job submitted for " + pauli_string)
    wait_for_job(job, pauli_string)

    results = job.get_results()
    
    if results is None:
       raise RuntimeError("Job finished but returned no results")
    dist = results["results"]

    # Debug
    print("Distinct states: " + str(len(dist)))
    for state, count in list(dist.items())[:10]:
        print("  " + str(list(state)) + " : " + str(count))

    total_valid = 0
    exp = 0
    for state, count in dist.items():
        s = list(state)
        if len(s) < 4:
            continue
        if (s[0] + s[1] != 1) or (s[2] + s[3] != 1):
            continue
        total_valid += count
        bits = [s[3], s[1]]
        parity = 1
        for i, p in enumerate(pauli_string):
            if p != "I" and bits[i] == 1:
                parity *= -1
        exp += parity * count

    print("total_valid=" + str(total_valid) + "  exp=" + str(exp))
    if total_valid == 0:
        return 0.0
    return exp / total_valid
def energy(theta):

    theta = theta[0]  
    shots = 1000 
    Z0 = measure_expectation(theta, "IZ", shots)
    Z1 = measure_expectation(theta, "ZI", shots)
    XX = measure_expectation(theta, "XX", shots)
    YY = measure_expectation(theta, "YY", shots)

    E = 3.89655 + 0.75187*Z0 - 4.64842*Z1 - 0.629235*(XX + YY)

    print(f"Theta: {theta:.4f}, Energy: {E:.6f}")

    return E

# Run this diagnostic cell first — don't run measure_expectation yet

import perceval as pcvl
from perceval import catalog
from pprint import pprint


BACKEND = "sim:ascella"

# ── 1. Check platform specs ───────────────────────────────────────────────
remote = pcvl.RemoteProcessor(BACKEND, token=TOKEN)
specs = remote.specs
print("Platform specs:")
pprint(specs)

# ── 2. Check how many modes the heralded CNOT actually uses ──────────────
cnot = catalog['heralded cnot'].build_processor()
print("\nHeralded CNOT modes: " + str(cnot.m))
print("Heralded CNOT circuit size: " + str(cnot.circuit_size))

# ── 3. Build your full local processor and check its size ────────────────
proc = pcvl.Processor("SLOS", 4)
proc.add([0, 1], catalog['x'].build_circuit())
proc.add([2, 3], catalog['ry'].build_circuit(theta=0.1))
proc.add([2, 3, 0, 1], catalog['heralded cnot'].build_processor())
print("\nFull processor logical modes: " + str(proc.m))
print("Full processor circuit size:  " + str(proc.circuit_size))

# ── 4. Check what linear_circuit looks like ──────────────────────────────
lc = proc.linear_circuit()
print("Linear circuit modes: " + str(lc.m))
optimal_theta = 0.2289515  # your known value in radians (or whatever units you used)

Z0 = measure_expectation(optimal_theta, "IZ")
Z1 = measure_expectation(optimal_theta, "ZI")
XX = measure_expectation(optimal_theta, "XX")
YY = measure_expectation(optimal_theta, "YY")

E = 3.89655 + 0.75187*Z0 - 4.64842*Z1 - 0.629235*(XX + YY)

print("IZ  = " + str(Z0))
print("ZI  = " + str(Z1))
print("XX  = " + str(XX))
print("YY  = " + str(YY))
print("Ground state energy = " + str(E))

