import sys
import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta
import os

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

# Configuration
file_path = sys.argv[1]  # Replace with actual file
channels = ["ditau", "mutau", "eletau"]
output_base = "plots"
os.makedirs(output_base, exist_ok=True)

# Define histogram sets
histogram_sets = [
    ("hTrig1PassPt", "hTrig1TotalPt", {"ditau":"hltHpsPFTauTrack","mutau":"hltL3crIsoL1TkSingleMu22","eletau":"hltEle30WPTightGsfTrackIsoL1SeededFilter"}, "pt"),
    ("hTrig2PassPt", "hTrig2TotalPt", {"ditau":"hltHpsDoublePFTau35MediumDitauWPDeepTau","mutau":"hltHpsPFTau27LooseTauWPDeepTau","eletau":"hltHpsPFTau30LooseTauWPDeepTau"}, "pt"),
    ("hTrig1PassEta", "hTrig1TotalEta", {"ditau":"hltHpsPFTauTrack","mutau":"hltL3crIsoL1TkSingleMu22","eletau":"hltEle30WPTightGsfTrackIsoL1SeededFilter"}, "eta"),
    ("hTrig2PassEta", "hTrig2TotalEta", {"ditau":"hltHpsDoublePFTau35MediumDitauWPDeepTau","mutau":"hltHpsPFTau27LooseTauWPDeepTau","eletau":"hltHpsPFTau30LooseTauWPDeepTau"}, "eta"),
]

# Start HTML output
html_lines = [
    "<html><head><title>Trigger Efficiencies</title></head><body>",
    "<h1>Trigger Efficiency Plots</h1>"
]

# Loop through each channel
for channel in channels:
    channel_dir = os.path.join(output_base, channel)
    os.makedirs(channel_dir, exist_ok=True)

    with uproot.open(file_path) as file:
        for pass_name, total_name, trig_label, axis_type in histogram_sets:
            h_pass = file[f"{channel}/{pass_name}"].to_numpy()
            h_total = file[f"{channel}/{total_name}"].to_numpy()

            values_pass, bin_edges = h_pass[0], h_pass[1]
            values_total = h_total[0]
            bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

            eff, err_low, err_up = binomial_confidence(values_pass, values_total)

            # Plot
            plt.figure(figsize=(8, 5))
            plt.errorbar(
                bin_centers,
                eff,
                yerr=[err_low, err_up],
                fmt='s',
                #linestyle='-',
                color='black',
                label=f"{trig_label[channel]} Efficiency"
            )

            axis_label = r"$p_T^\tau$ (GeV)" if axis_type == "pt" else r"$\eta^\tau$"
            plt.xlabel(axis_label)
            plt.ylabel("Efficiency")
            plt.title(f"{trig_label[channel]} Efficiency – {channel} ({axis_type})")
            plt.ylim(0, 1.1)
            plt.grid(True)
            plt.legend()
            plt.tight_layout()

            # Save plot
            plot_filename = f"{trig_label[channel]}_eff_{axis_type}.jpg"
            plot_path = os.path.join(channel_dir, plot_filename)
            plt.savefig(plot_path)
            plt.close()

            # Add image to HTML
            relative_path = os.path.join(channel, plot_filename)
            html_lines.append(f"<h2>{trig_label[channel]} – {channel} ({axis_type})</h2>")
            html_lines.append(f'<img src="{relative_path}" alt="{trig_label[channel]} {axis_type} {channel}"><br>')

# Finalize HTML
html_lines.append("</body></html>")
html_path = os.path.join(output_base, "efficiencies.html")
with open(html_path, "w") as f:
    f.write("\n".join(html_lines))

print(f"All efficiency plots saved and HTML file generated at '{html_path}'")
