#!/usr/bin/env python3
"""
Finite-Size Scaling Analysis for Spin Glass Phase Transition
"""

import struct
import sys
import os
import glob
import math

# Add Matplotlib for plotting
import matplotlib.pyplot as plt


def read_overlap_binary(filename, nbetas=16):
    """Read binary overlap file (β, <q^2>, <q^4>, F) of length 4*nbetas doubles."""
    with open(filename, 'rb') as f:
        fmt = f'{4 * nbetas}d'
        data_bytes = f.read(8 * 4 * nbetas)
        
        # If the file is shorter/longer than expected, adjust nbetas
        if len(data_bytes) != 8 * 4 * nbetas:
            actual_nbetas = len(data_bytes) // (8 * 4)
            nbetas = actual_nbetas
            fmt = f'{4 * nbetas}d'
        
        all_values = struct.unpack(fmt, data_bytes)
        
        betas = all_values[0:nbetas]
        Q2 = all_values[nbetas:2*nbetas]
        Q4 = all_values[2*nbetas:3*nbetas]
        F = all_values[3*nbetas:4*nbetas]
        
        return {
            'betas': list(betas),
            'Q2': list(Q2),
            'Q4': list(Q4),
            'F': list(F),
            'nbetas': nbetas
        }


def calculate_binder_cumulant(Q2, Q4):
    """
    Given lists Q2[i] = <q^2> and Q4[i] = <q^4> for i=0..nbetas-1,
    return a list g[i] = Q4[i] / (Q2[i])^2.
    """
    return [q4 / (q2 * q2) for q2, q4 in zip(Q2, Q4)]


def average_disorder_realizations(data_list):
    """
    data_list = [dict_i, dict_j, ...], each dict has keys:
      'betas', 'Q2', 'Q4', 'F', 'nbetas'.
    Returns:
      avg_data: { 'betas': [...], 'Q2': [...], 'Q4': [...], 'F': [...], 'nbetas': nbetas }
      errors:  { 'Q2': [...], 'Q4': [...], 'binder': [...] }
    where each list is length nbetas.  'errors["binder"][i]' is σ_g at β_i.
    """
    nbetas = data_list[0]['nbetas']
    n_samples = len(data_list)
    
    # Initialize summed arrays
    avg_data = {
        'betas': data_list[0]['betas'].copy(),
        'Q2': [0.0] * nbetas,
        'Q4': [0.0] * nbetas,
        'F': [0.0] * nbetas,
        'nbetas': nbetas
    }
    
    # Sum over all disorder realizations
    for data in data_list:
        for i in range(nbetas):
            avg_data['Q2'][i] += data['Q2'][i]
            avg_data['Q4'][i] += data['Q4'][i]
            avg_data['F'][i] += data['F'][i]
    
    # Divide by number of samples to get the mean
    for i in range(nbetas):
        avg_data['Q2'][i] /= n_samples
        avg_data['Q4'][i] /= n_samples
        avg_data['F'][i] /= n_samples
    
    # Prepare error‐arrays
    errors = {
        'Q2': [0.0] * nbetas,
        'Q4': [0.0] * nbetas,
        'binder': [0.0] * nbetas
    }
    
    # Compute sample variance for Q2 and Q4, then standard error
    for i in range(nbetas):
        var_q2 = sum((data['Q2'][i] - avg_data['Q2'][i])**2 for data in data_list) / (n_samples - 1)
        var_q4 = sum((data['Q4'][i] - avg_data['Q4'][i])**2 for data in data_list) / (n_samples - 1)
        
        sigma_q2 = math.sqrt(var_q2 / n_samples)
        sigma_q4 = math.sqrt(var_q4 / n_samples)
        
        errors['Q2'][i] = sigma_q2
        errors['Q4'][i] = sigma_q4
        
        # Now propagate error into g = Q4 / (Q2)^2
        q2 = avg_data['Q2'][i]
        q4 = avg_data['Q4'][i]
        # ∂g/∂q2 = -2 * q4 / q2^3
        dg_dq2 = -2.0 * q4 / (q2**3)
        # ∂g/∂q4 = 1 / q2^2
        dg_dq4 = 1.0 / (q2**2)
        
        errors['binder'][i] = math.sqrt((dg_dq2 * sigma_q2)**2 + (dg_dq4 * sigma_q4)**2)
    
    return avg_data, errors


def find_crossing_points(L_data):
    """
    L_data: dictionary where keys are integer L (e.g. 8,12,16…) and values are dicts with:
      'binder': [g_i for i=0..nbetas-1],
      'betas':  [β_i for i=0..nbetas-1].
    Returns a list of crossings; each crossing is a dict:
      { 'L1': L, 'L2': Lnext, 'beta': β_cross, 'g': g_cross }.
    """
    L_values = sorted(L_data.keys())
    crossings = []
    
    for idx in range(len(L_values) - 1):
        L1 = L_values[idx]
        L2 = L_values[idx + 1]
        
        binder1 = L_data[L1]['binder']
        binder2 = L_data[L2]['binder']
        betas = L_data[L1]['betas']
        
        # Walk through adjacent β‐bins to find sign changes in (g1 - g2)
        for j in range(len(betas) - 1):
            diff1 = binder1[j]   - binder2[j]
            diff2 = binder1[j+1] - binder2[j+1]
            if diff1 * diff2 < 0:
                # They cross between betas[j] and betas[j+1].  Do linear interpolation.
                x1, x2 = betas[j], betas[j+1]
                y1_1, y1_2 = binder1[j],   binder1[j+1]
                y2_1, y2_2 = binder2[j],   binder2[j+1]
                
                # Slope of g1 in [x1,x2]:
                m1 = (y1_2 - y1_1) / (x2 - x1)
                # Slope of g2 in [x1,x2]:
                m2 = (y2_2 - y2_1) / (x2 - x1)
                
                # Intercepts at x=0: b1 = y1_1 - m1*x1, b2 = y2_1 - m2*x1
                b1 = y1_1 - m1 * x1
                b2 = y2_1 - m2 * x1
                
                if m1 != m2:
                    beta_cross = (b2 - b1) / (m1 - m2)
                    g_cross = m1 * beta_cross + b1
                    crossings.append({
                        'L1': L1,
                        'L2': L2,
                        'beta': beta_cross,
                        'g': g_cross
                    })
    return crossings


def main():
    if len(sys.argv) < 2:
        print("Usage: python finite_size_scaling_plot.py /path/to/base/directory [Lz] [nbetas]")
        print("Example: python finite_size_scaling_plot.py /ada/.../MCQSG/output 2048 16")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    Lz = int(sys.argv[2]) if len(sys.argv) > 2 else 2048
    nbetas = int(sys.argv[3]) if len(sys.argv) > 3 else 16
    
    # 1) Discover which L‐folders exist under base_dir, chosen from a small list:
    L_values = []
    for L in [8, 12, 16, 20]:
        candidate = f"{base_dir}/{L}x{Lz}"
        if os.path.isdir(candidate):
            L_values.append(L)
    
    print(f"Found system sizes: {L_values}")
    print(f"Expected Lz: {Lz}, nbetas: {nbetas}\n")
    
    # Container to hold averaged data for each L
    L_data = {}
    
    # 2) For each L, find all overlap files under I000/BIT*/:
    for L in L_values:
        print(f"Processing L = {L} ...")
        pattern = f"{base_dir}/{L}x{Lz}/I000/BIT*/overlap.PBCZ.GPU"
        files = sorted(glob.glob(pattern))
        
        if not files:
            print(f"  No files found for L={L}")
            continue
        
        print(f"  Found {len(files)} disorder realizations for L={L}")
        
        all_data = []
        for fpath in files:
            try:
                data = read_overlap_binary(fpath, nbetas)
                all_data.append(data)
                bit_num = fpath.split('/BIT')[1].split('/')[0]
                print(f"    → Read BIT{bit_num}")
            except Exception as e:
                print(f"    → Error reading {fpath}: {e}")
        
        if not all_data:
            continue
        
        # 3) Average over disorder and compute errors
        avg_data, errors = average_disorder_realizations(all_data)
        binder = calculate_binder_cumulant(avg_data['Q2'], avg_data['Q4'])
        
        L_data[L] = {
            'data': avg_data,
            'errors': errors,
            'binder': binder,
            'betas': avg_data['betas'],
            'n_samples': len(all_data)
        }
    
    # 4) Write out finite_size_data.txt exactly as before
    with open('finite_size_data.txt', 'w') as fout:
        fout.write("# Finite-size scaling data for spin glass\n")
        fout.write(f"# System sizes: {sorted(L_data.keys())}\n")
        fout.write("# Format: L beta Q2 Q4 Binder error_Binder\n\n")
        
        for L in sorted(L_data.keys()):
            fout.write(f"# L = {L} (averaged over {L_data[L]['n_samples']} disorder realizations)\n")
            npts = L_data[L]['data']['nbetas']
            for i in range(npts):
                β_i = L_data[L]['betas'][i]
                Q2_i = L_data[L]['data']['Q2'][i]
                Q4_i = L_data[L]['data']['Q4'][i]
                g_i  = L_data[L]['binder'][i]
                σg_i = L_data[L]['errors']['binder'][i]
                fout.write(f"{L} {β_i:.8f} {Q2_i:.10f} {Q4_i:.10f} {g_i:.10f} {σg_i:.10f}\n")
            fout.write("\n")
    
    # 5) Find and print crossing points
    crossings = find_crossing_points(L_data)
    
    print("\n" + "="*60)
    print("FINITE-SIZE SCALING ANALYSIS RESULTS")
    print("="*60)
    for L in sorted(L_data.keys()):
        betas = L_data[L]['betas']
        binder = L_data[L]['binder']
        σbinder = L_data[L]['errors']['binder']
        
        print(f"\nL = {L}:")
        print(f"  Number of disorder realizations: {L_data[L]['n_samples']}")
        print(f"  Beta range: [{min(betas):.4f}, {max(betas):.4f}]")
        
        # Show Binder at β ≈ 0.29 for reference (closest point)
        target_beta = 0.29
        idx = min(range(len(betas)), key=lambda i: abs(betas[i] - target_beta))
        print(f"  Binder at β≈{target_beta:.2f}: {binder[idx]:.4f} "
              f"± {σbinder[idx]:.4f}")
    
    if crossings:
        print("\n" + "-"*40)
        print("CROSSING POINTS:")
        print("-"*40)
        for c in crossings:
            print(f"  L={c['L1']} × L={c['L2']}: β_c = {c['beta']:.6f}, g_c = {c['g']:.6f}")
        
        avg_beta = sum(c['beta'] for c in crossings) / len(crossings)
        avg_g = sum(c['g'] for c in crossings) / len(crossings)
        print(f"\n  Average critical point:")
        print(f"    β_c = {avg_beta:.6f}")
        print(f"    g_c = {avg_g:.6f}")
    
    # 6) Now: produce a Matplotlib figure and save as 'binder_scaling.png'
    # Prepare one plot with error bars for each L
    plt.figure(figsize=(8, 6))
    
    for L in sorted(L_data.keys()):
        betas = L_data[L]['betas']
        binder = L_data[L]['binder']
        σbinder = L_data[L]['errors']['binder']
        
        # Plot with marker 'o', solid line, vertical error bars
        plt.errorbar(
            betas, binder, yerr=σbinder,
            fmt='o-',  # circle markers joined by solid lines
            label=f"L = {L}",
            capsize=3,  # small horizontal cap on each error bar
            markersize=4
        )
    
    plt.xlabel("β (inverse temperature)", fontsize=12)
    plt.ylabel("Binder cumulant $g = \\langle q^4\\rangle / \\langle q^2\\rangle^2$", fontsize=12)
    plt.title("Finite-Size Scaling of Binder Cumulant", fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc='best', fontsize=10)
    plt.tight_layout()
    
    # Save the figure
    output_png = "binder_scaling.png"
    plt.savefig(output_png, dpi=300)
    plt.close()
    
    print(f"\nSaved Binder‐vs‐β plot as: {output_png}")


if __name__ == "__main__":
    main()
