"""Analysis of the reweighting test (see run_reweighting_test.sh).

Compares direct Nf=4 averages with Nf=2 averages reweighted by det(D)^2,
    <O>_{Nf=4} = <O det(D)^2>_{Nf=2} / <det(D)^2>_{Nf=2},
for the average plaquette and log|det D|, and checks <exp(-dH)> = 1 for both chains.

Usage: python3 rw_analysis.py <data directory>
The directory contains files nf2_<tag>_c<i>.dat and nf4_<tag>_c<i>.dat written by rw_chain, with columns
    plaquette    dH    Re log det(D)    Im log det(D) (mod 2 pi)
The exit code is 0 if all checks pass and 1 otherwise.
"""
import glob
import os
import re
import sys

import numpy as np

NBLK = 100        # jackknife blocks
MAX_DEV = 3.0     # maximum allowed deviation in units of the error
MIN_SEP = 5.0     # minimum separation between Nf=2 and Nf=4, otherwise the test has no power


def load(pattern):
    files = sorted(glob.glob(pattern))
    return np.concatenate([np.loadtxt(f, ndmin=2) for f in files]), len(files)


def jackknife(estimator, columns):
    """Block jackknife estimate and error of estimator(*columns)."""
    n = len(columns[0]) // NBLK * NBLK
    columns = [c[:n] for c in columns]
    blk = n // NBLK
    values = []
    for b in range(NBLK):
        keep = np.ones(n, bool)
        keep[b * blk:(b + 1) * blk] = False
        values.append(estimator(*[c[keep] for c in columns]))
    values = np.array(values)
    error = np.sqrt((NBLK - 1) / NBLK * np.sum((values - values.mean()) ** 2))
    return estimator(*columns), error


def mean(o):
    return o.mean()


def mean_exp_minus_dH(dH):
    return np.exp(-dH).mean()


def reweighted_mean(o, logdet):
    w = np.exp(2.0 * (logdet - logdet.max()))  # det(D)^2 = |det D|^2, since det D is real
    return np.sum(o * w) / np.sum(w)


def fmt(value, error):
    return f"{value:.6f} +/- {error:.6f}"


def main(directory):
    tags = sorted({re.sub(r"^nf2_(.*)_c\d+\.dat$", r"\1", os.path.basename(f))
                   for f in glob.glob(os.path.join(directory, "nf2_*_c*.dat"))})
    if not tags:
        print(f"No data files in {directory}")
        return 1
    passed = True
    for tag in tags:
        nf2, n2 = load(os.path.join(directory, f"nf2_{tag}_c*.dat"))
        nf4, n4 = load(os.path.join(directory, f"nf4_{tag}_c*.dat"))
        print(f"=== {tag}: {len(nf2)} trajectories for Nf=2 ({n2} chains), {len(nf4)} for Nf=4 ({n4} chains)")

        for name, data in [("Nf=2", nf2), ("Nf=4", nf4)]:
            value, error = jackknife(mean_exp_minus_dH, [data[:, 1]])
            dev = (value - 1.0) / error
            ok = abs(dev) < MAX_DEV
            passed &= ok
            print(f"  {name}: <exp(-dH)> = {fmt(value, error)}  ({dev:+.2f} sigma from 1)  "
                  f"acceptance = {np.minimum(1, np.exp(-data[:, 1])).mean():.3f}  {'ok' if ok else 'FAIL'}")
            # det D must be real: Im log det is 0 or pi
            im = np.abs(data[:, 3])
            max_im = np.max(np.minimum(im, np.abs(im - np.pi)))
            ok = max_im < 1e-8
            passed &= ok
            print(f"        det D real: max |Im log det D| (mod pi) = {max_im:.1e}  {'ok' if ok else 'FAIL'}")

        logdet = nf2[:, 2]
        w = np.exp(2.0 * (logdet - logdet.max()))
        print(f"  reweighting effective sample size: {w.sum() ** 2 / (w ** 2).sum():.0f} of {len(w)}")

        for obs, col in [("plaquette", 0), ("log|det D|", 2)]:
            direct2 = jackknife(mean, [nf2[:, col]])
            direct4 = jackknife(mean, [nf4[:, col]])
            reweighted = jackknife(reweighted_mean, [nf2[:, col], logdet])
            dev = (direct4[0] - reweighted[0]) / np.hypot(direct4[1], reweighted[1])
            sep = (direct4[0] - direct2[0]) / np.hypot(direct4[1], direct2[1])
            ok = abs(dev) < MAX_DEV and abs(sep) > MIN_SEP
            passed &= ok
            print(f"  {obs}:")
            print(f"      Nf=2 direct              {fmt(*direct2)}")
            print(f"      Nf=4 direct              {fmt(*direct4)}")
            print(f"      Nf=2 reweighted to Nf=4  {fmt(*reweighted)}")
            print(f"      (direct - reweighted)/sigma = {dev:+.2f},  (Nf=4 - Nf=2)/sigma = {sep:+.1f}  "
                  f"{'ok' if ok else 'FAIL'}")

    print("PASSED" if passed else "FAILED")
    return 0 if passed else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    sys.exit(main(sys.argv[1]))
