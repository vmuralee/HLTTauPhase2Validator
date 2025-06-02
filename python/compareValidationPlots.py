import sys
import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta
import os

# Binomial efficiency with asymmetric errors
def binomial_confidence(n_pass, n_total, conf_level=0.6827):
    alpha = 1 - conf_level
    eff = np.divide(n_pass, n_total, out=np.zeros_like(n_pass, dtype=float), where=n_total != 0)

    lo = beta.ppf(alpha / 2, n_pass, n_total - n_pass + 1)
    hi = beta.ppf(1 - alpha / 2, n_pass + 1, n_total - n_pass)

    lo = np.nan_to_num(lo, nan=0.0)
    hi = np.nan_to_num(hi, nan=1.0)

    err_low = eff - lo
    err_up = hi - eff

    err_low[n_total == 0] = 0
    err_up[n_total == 0] = 0

    return eff, err_low, err_up

# File paths — set these to your actual files
file1_path = sys.argv[1]
file2_path = sys.argv[2]
file1_label = "Dataset A"
file2_label = "Dataset B"

channels = ["ditau", "mutau", "eletau"]
output_base = "compare_plots"
os.makedirs(output_base, exist_ok=True)

histogram_sets = [
    ("hTrig1PassPt", "hTrig1TotalPt", "Trig1", "pt"),
    ("hTrig2PassPt", "hTrig2TotalPt", "Trig2", "pt"),
    ("hTrig1PassEta", "hTrig1TotalEta", "Trig1", "eta"),
    ("hTrig2PassEta", "hTrig2TotalEta", "Trig2", "eta"),
]

# Initialize HTML
html_lines = [
    "<html><head><title>Trigger Efficiency Comparison</title></head><body>",
    "<h1>Trigger Efficiency: Comparison Between Files</h1>"
]

# Open both ROOT files
file1 = uproot.open(file1_path)
file2 = uproot.open(file2_path)

# Loop over each channel
for channel in channels:
    channel_dir = os.path.join(output_base, channel)
    os.makedirs(channel_dir, exist_ok=True)

    for pass_name, total_name, trig_label, axis_type in histogram_sets:
        # Read histograms from both files
        h1_pass = file1[f"{channel}/{pass_name}"].to_numpy()
        h1_total = file1[f"{channel}/{total_name}"].to_numpy()

        h2_pass = file2[f"{channel}/{pass_name}"].to_numpy()
        h2_total = file2[f"{channel}/{total_name}"].to_numpy()

        values1_pass, bin_edges = h1_pass[0], h1_pass[1]
        values1_total = h1_total[0]
        values2_pass = h2_pass[0]
        values2_total = h2_total[0]

        bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

        eff1, err1_low, err1_up = binomial_confidence(values1_pass, values1_total)
        eff2, err2_low, err2_up = binomial_confidence(values2_pass, values2_total)

        # Plot comparison
        plt.figure(figsize=(8, 5))

        plt.errorbar(
            bin_centers, eff1, yerr=[err1_low, err1_up],
            fmt='o', linestyle='-', color='blue', label=file1_label
        )
        plt.errorbar(
            bin_centers, eff2, yerr=[err2_low, err2_up],
            fmt='s', linestyle='--', color='red', label=file2_label
        )

        axis_label = r"$p_T^\tau$ (GeV)" if axis_type == "pt" else r"$\eta^\tau$"
        plt.xlabel(axis_label)
        plt.ylabel("Efficiency")
        plt.title(f"{trig_label} Efficiency Comparison – {channel} ({axis_type})")
        plt.ylim(0, 1.1)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()

        # Save plot
        plot_filename = f"{trig_label}_eff_{axis_type}.png"
        plot_path = os.path.join(channel_dir, plot_filename)
        plt.savefig(plot_path)
        plt.close()

        # Add to HTML
        relative_path = os.path.join(channel, plot_filename)
        html_lines.append(f"<h2>{trig_label} – {channel} ({axis_type})</h2>")
        html_lines.append(f'<img src="{relative_path}" alt="{trig_label} {channel} {axis_type}"><br>')

# Finalize HTML
html_lines.append("</body></html>")
html_path = os.path.join(output_base, "efficiencies.html")
with open(html_path, "w") as f:
    f.write("\n".join(html_lines))

print(f"Comparison plots saved. Open '{html_path}' to view.")
