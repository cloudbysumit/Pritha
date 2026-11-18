#!/usr/bin/env python
"""
Automatic Lipid Droplet Circles (circle search) - Python version of the v4 ImageJ macro.

Finds thin bright RINGS (lipid droplets), fits a circle to each, and measures the mean
intensity of the ORIGINAL image inside every circle.

Requires: numpy, scipy, scikit-image, pandas, matplotlib, tifffile
          (czifile only if you open .czi files:  pip install czifile)

Usage:
    python lipid_droplet_circle_finder.py                       # opens a file chooser
    python lipid_droplet_circle_finder.py image.tif
    python lipid_droplet_circle_finder.py image.czi --channel 0 --z 2 --t 50
    python lipid_droplet_circle_finder.py image.tif --min-contrast 4 --min-brightness 4
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import io as skio
from skimage.morphology import h_maxima
from skimage.measure import label, regionprops


# ----------------------------------------------------------------------
# Loading
# ----------------------------------------------------------------------
def load_plane(path, channel=0, z=0, t=0):
    """Load one 2-D plane (channel / z / t) from a TIFF, CZI or ordinary image file."""
    ext = os.path.splitext(path)[1].lower()

    if ext == ".czi":
        try:
            import czifile
        except ImportError:
            sys.exit("To open .czi files run:  pip install czifile")
        with czifile.CziFile(path) as czi:
            data = czi.asarray()
            axes = czi.axes            # e.g. 'BTCZYX0' or 'TZCYX0'
        return _pick_plane(data, axes, channel, z, t)

    if ext in (".tif", ".tiff"):
        import tifffile
        with tifffile.TiffFile(path) as tf:
            data = tf.asarray()
            axes = tf.series[0].axes   # e.g. 'TZCYX'
        return _pick_plane(data, axes, channel, z, t)

    img = skio.imread(path)
    if img.ndim == 3:                  # RGB / multi-channel picture -> take chosen channel
        img = img[..., min(channel, img.shape[-1] - 1)]
    return img


def _pick_plane(data, axes, channel, z, t):
    """Reduce an N-D array to a 2-D plane using axis letters (T, Z, C, Y, X)."""
    index = []
    for ax, n in zip(axes, data.shape):
        if ax == "C":
            index.append(min(channel, n - 1))
        elif ax in ("Z", "S"):
            index.append(min(z, n - 1))
        elif ax == "T":
            index.append(min(t, n - 1))
        elif ax in ("Y", "X"):
            index.append(slice(None))
        else:                          # any other axis (e.g. sample/RGB, singleton)
            index.append(0)
    plane = data[tuple(index)]
    return np.squeeze(plane)


# ----------------------------------------------------------------------
# Ring detector kernel (same as buildKernel() in the macro)
# ----------------------------------------------------------------------
def build_kernel(r, ring_half, gap, flank_w):
    half = int(np.ceil(r + ring_half + gap + flank_w))
    y, x = np.mgrid[-half:half + 1, -half:half + 1]
    d = np.hypot(x, y)

    inner_out = r - ring_half - gap
    outer_in = r + ring_half + gap

    ring = np.abs(d - r) <= ring_half
    flank = (~ring) & (
        ((d >= inner_out - flank_w) & (d < inner_out))
        | ((d > outer_in) & (d <= outer_in + flank_w))
    )

    k = np.zeros_like(d)
    k[ring] = 1.0 / ring.sum()
    k[flank] = -1.0 / flank.sum()
    return k                                      # sums to zero


# ----------------------------------------------------------------------
# Ring test (same as ringTest() in the macro)
# ----------------------------------------------------------------------
def _sample(img, xs, ys):
    return ndi.map_coordinates(img, [ys, xs], order=1, mode="nearest")


def ring_test(sm, x0, y0, rr, off, min_gap, n_dir=24):
    """Return (fraction of directions where ring beats BOTH sides by min_gap, mean ring value)."""
    ang = np.arange(n_dir) * 2 * np.pi / n_dir
    ca, sa = np.cos(ang), np.sin(ang)

    ring_v = _sample(sm, x0 + rr * ca, y0 + rr * sa)
    v2 = _sample(sm, x0 + (rr - 1) * ca, y0 + (rr - 1) * sa)
    v3 = _sample(sm, x0 + (rr + 1) * ca, y0 + (rr + 1) * sa)
    ring_v = np.maximum(ring_v, np.maximum(v2, v3))

    in_v = _sample(sm, x0 + (rr - off) * ca, y0 + (rr - off) * sa)
    out_v = _sample(sm, x0 + (rr + off) * ca, y0 + (rr + off) * sa)
    side_max = np.maximum(in_v, out_v)

    good = np.count_nonzero(ring_v > side_max + min_gap)
    return good / n_dir, ring_v.mean()


# ----------------------------------------------------------------------
# Main detection
# ----------------------------------------------------------------------
def detect_droplets(
    img,
    r_min=4, r_max=16, r_step=1,
    thickness=3, sigma=1.0,
    thr_pct=15, min_comp=0.75, min_space=0.70,
    contrast_pct=8, bright_pct=12, margin=0.5,
    verbose=True,
):
    ring_half = thickness / 2.0
    gap, flank_w = 1, 3
    off = ring_half + 2

    img = img.astype(np.float32)
    sm = ndi.gaussian_filter(img, sigma) if sigma > 0 else img.copy()
    H, W = sm.shape

    # reference "bright level" = 99.5th percentile of the smoothed image
    ref = np.percentile(sm, 99.5)
    min_contrast_abs = contrast_pct / 100.0 * ref
    min_bright_abs = bright_pct / 100.0 * ref

    # ring response for every radius
    radii = np.arange(r_min, r_max + 1, r_step)
    resp = np.empty((len(radii), H, W), dtype=np.float32)
    for k, r in enumerate(radii):
        if verbose:
            print(f"  searching radius {r} px ...", end="\r")
        resp[k] = ndi.correlate(sm, build_kernel(r, ring_half, gap, flank_w), mode="nearest")
    if verbose:
        print()

    max_resp_map = resp.max(axis=0)
    max_resp = float(max_resp_map.max())
    thr = max(max_resp * thr_pct / 100.0, min_contrast_abs)

    # candidate centres = local maxima with prominence >= thr/2
    cand = []
    if max_resp > 0:
        peaks = h_maxima(max_resp_map, thr / 2.0)
        lab = label(peaks)
        for p in regionprops(lab):
            cy, cx = p.centroid
            cand.append((int(round(cx)), int(round(cy))))

    # best radius and score for each candidate
    cands = []
    for cx, cy in cand:
        profile = resp[:, cy, cx]
        k = int(profile.argmax())
        cands.append(dict(x=cx, y=cy, r=float(radii[k]), score=float(profile[k])))

    # first ring test
    kept = []
    for c in cands:
        if c["score"] < thr:
            continue
        x0, y0, rr = c["x"], c["y"], c["r"]
        if (x0 - rr - ring_half < 0 or y0 - rr - ring_half < 0
                or x0 + rr + ring_half >= W or y0 + rr + ring_half >= H):
            continue
        frac, ring_mean = ring_test(sm, x0, y0, rr, off, margin * c["score"])
        if frac >= min_comp and ring_mean >= min_bright_abs:
            kept.append(c)

    # remove duplicates (strongest first)
    kept.sort(key=lambda c: c["score"], reverse=True)
    accepted = []
    for c in kept:
        ok = True
        for a in accepted:
            d = np.hypot(c["x"] - a["x"], c["y"] - a["y"])
            if d < min_space * (c["r"] + a["r"]):
                ok = False
                break
        if ok:
            accepted.append(c)

    # refine each circle onto the bright ring (centre +/-2 px, radius +/-2 px), then test again
    n_ang = 36
    ang = np.arange(n_ang) * 2 * np.pi / n_ang
    ca, sa = np.cos(ang), np.sin(ang)
    final = []
    for c in accepted:
        best = (-1e30, c["x"], c["y"], c["r"])
        for dy in range(-2, 3):
            for dx in range(-2, 3):
                for rc in np.arange(c["r"] - 2, c["r"] + 2.01, 0.5):
                    if rc < 2:
                        continue
                    px, py = c["x"] + dx, c["y"] + dy
                    if (px - rc - ring_half < 0 or py - rc - ring_half < 0
                            or px + rc + ring_half >= W or py + rc + ring_half >= H):
                        continue
                    ps = _sample(sm, px + rc * ca, py + rc * sa).sum()
                    if ps > best[0]:
                        best = (ps, px, py, rc)
        _, px0, py0, rc0 = best
        frac, ring_mean = ring_test(sm, px0, py0, rc0, off, margin * c["score"])
        if frac >= min_comp and ring_mean >= min_bright_abs:
            final.append(dict(x=px0, y=py0, r=rc0, score=c["score"]))

    if verbose:
        print(f"  candidates: {len(cands)} | after first test: {len(accepted)} | "
              f"final droplets: {len(final)}")
    return final


# ----------------------------------------------------------------------
# Measurement on the ORIGINAL image
# ----------------------------------------------------------------------
def measure(img, droplets, roi_offset=0.0, pixel_size=None):
    H, W = img.shape
    yy, xx = np.mgrid[:H, :W]
    rows = []
    for i, d in enumerate(droplets, start=1):
        rr = max(d["r"] + roi_offset, 1.0)
        # ImageJ draws the circle at (x + 0.5, y + 0.5) and tests pixel centres at +0.5,
        # so the pixel-index distance to (x, y) is the right test.
        mask = (xx - d["x"]) ** 2 + (yy - d["y"]) ** 2 <= rr * rr
        vals = img[mask].astype(np.float64)
        area_px = int(mask.sum())
        row = {
            "ROI": f"LD_{i}",
            "X_px": d["x"] + 0.5,
            "Y_px": d["y"] + 0.5,
            "Radius_px": rr,
            "Area_px": area_px,
            "Mean": vals.mean() if area_px else np.nan,
            "Min": vals.min() if area_px else np.nan,
            "Max": vals.max() if area_px else np.nan,
            "IntDen": vals.sum(),
            "Perimeter_px": 2 * np.pi * rr,
        }
        if pixel_size:
            row["Area_um2"] = area_px * pixel_size ** 2
            row["Radius_um"] = rr * pixel_size
        rows.append(row)
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# Overlay picture
# ----------------------------------------------------------------------
def save_overlay(img, droplets, roi_offset, out_png, show=False):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lo, hi = np.percentile(img, (1, 99.8))
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(img, cmap="afmhot", vmin=lo, vmax=hi)
    for i, d in enumerate(droplets, start=1):
        rr = max(d["r"] + roi_offset, 1.0)
        ax.add_patch(plt.Circle((d["x"], d["y"]), rr, fill=False, ec="yellow", lw=1))
        ax.text(d["x"], d["y"], str(i), color="white", fontsize=6, ha="center", va="center")
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    if show:
        plt.show()
    plt.close(fig)


# ----------------------------------------------------------------------
# Command line
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Lipid droplet circle finder (ring based)")
    ap.add_argument("image", nargs="?", help="image file (.tif, .czi, .png ...)")
    ap.add_argument("--channel", type=int, default=0, help="channel index (0-based)")
    ap.add_argument("--z", type=int, default=0, help="z index (0-based)")
    ap.add_argument("--t", type=int, default=0, help="time index (0-based)")
    ap.add_argument("--pixel-size", type=float, default=None, help="um per pixel (adds um columns)")

    ap.add_argument("--min-radius", type=float, default=4)
    ap.add_argument("--max-radius", type=float, default=16)
    ap.add_argument("--radius-step", type=int, default=1)
    ap.add_argument("--ring-thickness", type=float, default=3)
    ap.add_argument("--sigma", type=float, default=1.0)
    ap.add_argument("--detection-strength", type=float, default=15)
    ap.add_argument("--completeness", type=float, default=0.75)
    ap.add_argument("--min-spacing", type=float, default=0.70)
    ap.add_argument("--min-contrast", type=float, default=8)
    ap.add_argument("--min-brightness", type=float, default=12)
    ap.add_argument("--margin", type=float, default=0.5)
    ap.add_argument("--roi-offset", type=float, default=0.0,
                    help="pixels added to the ring centre-line radius (negative = smaller circle)")
    ap.add_argument("--show", action="store_true", help="show the overlay window")
    args = ap.parse_args()

    path = args.image
    if not path:
        import tkinter as tk
        from tkinter import filedialog
        tk.Tk().withdraw()
        path = filedialog.askopenfilename(title="Select the lipid droplet image")
        if not path:
            sys.exit("No image selected.")

    print(f"Loading {path}")
    img = load_plane(path, args.channel, args.z, args.t)
    print(f"  plane size: {img.shape[1]} x {img.shape[0]} px, dtype {img.dtype}")

    droplets = detect_droplets(
        img,
        r_min=args.min_radius, r_max=args.max_radius, r_step=args.radius_step,
        thickness=args.ring_thickness, sigma=args.sigma,
        thr_pct=args.detection_strength, min_comp=args.completeness,
        min_space=args.min_spacing, contrast_pct=args.min_contrast,
        bright_pct=args.min_brightness, margin=args.margin,
    )

    if not droplets:
        sys.exit("No droplets detected. Lower --min-contrast, --min-brightness or "
                 "--completeness, or check the radius range.")

    df = measure(img, droplets, args.roi_offset, args.pixel_size)

    base = os.path.splitext(path)[0]
    csv_path, png_path = base + "_droplets.csv", base + "_droplets.png"
    df.to_csv(csv_path, index=False)
    save_overlay(img, droplets, args.roi_offset, png_path, show=args.show)

    print(f"\nDroplets detected: {len(df)}")
    print(f"Mean of droplet means: {df['Mean'].mean():.3f}")
    print(f"Results table : {csv_path}")
    print(f"Overlay image : {png_path}")


if __name__ == "__main__":
    main()
