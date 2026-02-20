import os, re, glob, csv, statistics
import numpy as np
import matplotlib.pyplot as plt

# Parse all tf_np*_*.out files, keep only those with timings, and take median per NP
pat = re.compile(r"^tf_np(\d+)_.*\.out$")
groups = {}     # np -> list of elapsed times
ntasks_map = {} # np -> ntasks

for f in glob.glob("tf_np*_*.out"):
    base = os.path.basename(f)
    m = pat.match(base)
    if not m:
        continue
    np_ranks = int(m.group(1))
    txt = open(f, "r", errors="ignore").read()

    mN = re.search(r"Number of settings:\s*([0-9]+)", txt)
    mT = re.search(r"Elapsed time\s*:\s*([0-9]*\.?[0-9]+)", txt)
    if not (mN and mT):
        continue

    ntasks = int(mN.group(1))
    elapsed = float(mT.group(1))

    groups.setdefault(np_ranks, []).append(elapsed)
    ntasks_map[np_ranks] = ntasks  # same for all runs of same np in your case

if not groups:
    raise SystemExit("No valid tf_np*_*.out files found (missing 'Number of settings'/'Elapsed time').")

rows = []
for np_ranks, times in sorted(groups.items()):
    nworkers = np_ranks - 1
    ntasks = ntasks_map[np_ranks]
    elapsed_med = statistics.median(times)
    ttask = elapsed_med * nworkers / ntasks
    rows.append((np_ranks, nworkers, ntasks, elapsed_med, ttask, len(times)))

# Write CSV
with open("scaling_sorted.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["np","nworkers","ntasks","elapsed_s","ttask_s","nruns"])
    for r in rows:
        w.writerow(r)

# Arrays for plots
N = np.array([r[1] for r in rows], dtype=float)      # workers
T = np.array([r[3] for r in rows], dtype=float)      # wall time (median)
Ttask = np.array([r[4] for r in rows], dtype=float)

# Baseline: 1 worker (np=2)
i1 = np.where(N == 1)[0][0]
T1 = T[i1]
speedup = T1 / T
eff = speedup / N

# Amdahl fit: T/T1 = f + (1-f)/N
y = T / T1
fs_grid = np.linspace(0.0, 0.5, 5001)
mse = []
for fs in fs_grid:
    yhat = fs + (1.0 - fs) / N
    mse.append(np.mean((y - yhat) ** 2))
fs = fs_grid[int(np.argmin(mse))]
fp = 1.0 - fs

Nfine = np.logspace(np.log2(N.min()), np.log2(N.max()), 300, base=2)
speed_theory = 1.0 / (fs + (1.0 - fs) / Nfine)

print(f"Estimated serial fraction f_s ≈ {fs:.4f}")
print(f"Estimated parallel fraction ≈ {fp:.4f}")
print("np  workers  elapsed_s(median)  nruns")
for np_ranks, nworkers, _, elapsed_med, _, nruns in rows:
    print(f"{np_ranks:2d}  {nworkers:7d}  {elapsed_med:14.4f}  {nruns:5d}")

# Plots
plt.figure()
plt.plot(N, Ttask * 1e6, marker="o")
plt.xscale("log", base=2)
plt.xlabel("N_workers")
plt.ylabel("T_task [microseconds]")
plt.title("Strong scaling: per-task time")
plt.grid(True)
plt.savefig("scaling_ttask.png", dpi=200)

plt.figure()
plt.plot(N, speedup, marker="o", label="measured (median)")
plt.plot(Nfine, speed_theory, label=f"Amdahl fit (f_s={fs:.3f})")
plt.xscale("log", base=2)
plt.xlabel("N_workers")
plt.ylabel("Speedup")
plt.title("Strong scaling: speedup with Amdahl fit")
plt.grid(True)
plt.legend()
plt.savefig("scaling_speedup_amdahl.png", dpi=200)

plt.figure()
plt.plot(N, eff, marker="o")
plt.xscale("log", base=2)
plt.xlabel("N_workers")
plt.ylabel("Parallel efficiency")
plt.title("Strong scaling: efficiency")
plt.grid(True)
plt.savefig("scaling_efficiency.png", dpi=200)
