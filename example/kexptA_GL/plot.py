import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob, sys

# Function to parse one file and return dataframe
def parse_simulation_file(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('#') or line.startswith('[') or line.strip() == '':
                continue  # skip headers/comments
            tokens = line.strip().split()
            if len(tokens) > 17:
                continue  # ignore malformed lines
            try:
                entry = {
                    'runs': int(tokens[0]),
                    'nx': int(tokens[1]),
                    'kdim': int(tokens[2]),
                    'dt': float(tokens[3]),
                    'tol': float(tokens[4]),
                    'n_acc': float(tokens[5]),
                    'n_rej': float(tokens[6]),
                    'kp': float(tokens[7]),
                    'err_l2': float(tokens[8]),
                    'rk_time_mean': float(tokens[9]),
                    'rk_time_std': float(tokens[10]),
                    'krylov_time_mean': float(tokens[11]),
                    'krylov_time_std': float(tokens[12]),
                    'elapsed_time': float(tokens[13])
                }
                data.append(entry)
            except Exception as e:
                print(f"Error parsing line in {filename}: {line}")
                print(e)
    return pd.DataFrame(data)

# Load all matching files
all_files = glob.glob("../../kexpm*.txt")  # Adjust the pattern as needed
df_list = [parse_simulation_file(f) for f in all_files]
df = pd.concat(df_list, ignore_index=True)

# Add derived columns
df['tol_sci'] = df['tol'].apply(lambda x: f"{x:.0e}")
df['rk_total'] = df['rk_time_mean'] + df['rk_time_std']
df['n_total'] = df['n_acc'] + df['n_rej']

# Format tolerance in scientific notation for grouping
df['tol_sci'] = df['tol'].apply(lambda x: f"{x:.0e}")

# Normalize dt for alpha scaling (log scale)
dt_min = df['dt'].min()
dt_max = 1.0  # as per your spec: max dt for full alpha = 1

def compute_alpha(dt):
    # Map log(dt_min) -> 0.4, log(1) -> 1.0
    log_dt = np.log10(dt)
    log_min = np.log10(dt_min)
    log_max = 0  # log10(1) = 0
    alpha = 0.4 + 0.6 * (log_dt - log_min) / (log_max - log_min)
    return np.clip(alpha, 0.4, 1.0)

# Begin the combined plot
plt.figure(figsize=(10, 6))

# Loop over nx values
for nx_val in sorted(df['nx'].unique()):
    subset_nx = df[df['nx'] == nx_val]

    # Loop over dt values
    for dt_val in sorted(subset_nx['dt'].unique()):
        subset = subset_nx[subset_nx['dt'] == dt_val]
        grouped = subset.groupby('tol').mean(numeric_only=True).reset_index()
        alpha = compute_alpha(dt_val)

        plt.plot(
            grouped['kp'], grouped['tol'],
            marker='o',
            linestyle='-',
            label=f'nx={nx_val}, dt={dt_val:.0e}',
            alpha=alpha
        )

# Set log scale and labels
plt.yscale('log')
plt.xlabel('Average Krylov Steps (kp)')
plt.ylabel('Tolerance (tol)')
plt.title('kp vs tol for different nx and dt (alpha = dt scale)')
plt.grid(True)
plt.legend(fontsize=8, loc='best', ncol=2)
plt.tight_layout()
tols = df['tol_sci'].unique()
ncols = 3
nrows = int(np.ceil(len(tols) / ncols))

unique_tols = sorted(df['tol'].unique())
ncols = 3
nrows = int(np.ceil(len(unique_tols) / ncols))

fig, axs = plt.subplots(nrows, ncols, figsize=(16, 4 * nrows))

for idx, tol in enumerate(unique_tols[::-1]):
    ax = axs.flat[idx]
    sub_df = df[df['tol'] == tol]
    
    for i, nx_val in enumerate(sorted(sub_df['nx'].unique())):
        d = sub_df[sub_df['nx'] == nx_val].sort_values('dt')
        
        color = plt.cm.tab10(i % 10)  # Use color cycle for consistency

        # === Krylov plot ===
        ax.plot(
            d['dt'], d['rk_time_mean']/d['krylov_time_mean'], 
            label=f'Krylov nx={nx_val}', 
            color=color, linestyle='-'
        )

    ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.set_title(f'tol = {tol:.0e}')
    ax.set_xlabel('dt')
    ax.set_ylabel('Speedup (Krylov vs. RK)')
    ax.grid(True)
    ax.legend()

# Remove unused subplots
for j in range(len(unique_tols), nrows * ncols):
    fig.delaxes(axs.flat[j])

plt.tight_layout()

fig, axs = plt.subplots(nrows, ncols, figsize=(15, 4 * nrows), sharex=False, sharey=False)

for i, tol in enumerate(sorted(tols, key=lambda x: float(x))):
    ax = axs[i // ncols, i % ncols]
    sub_df = df[df['tol_sci'] == tol]
    sub_df = sub_df.sort_values(by='dt')
    for nx_val in sorted(sub_df['nx'].unique()):
        subset = sub_df[sub_df['nx'] == nx_val]
        n_total = subset['n_acc'] + subset['n_rej']
        ax.plot(subset['dt'], n_total, marker='o', label=f'nx={nx_val} (n_acc + n_rej)')
        ax.plot(subset['dt'], subset['kp'], marker='x', linestyle='--', label=f'nx={nx_val} (kp)')
    ax.set_xscale('log')
    ax.set_title(f"tol = {tol}")
    ax.set_xlabel("dt")
    ax.set_ylabel("n_acc + n_rej and kp")
    ax.grid(True)
    ax.legend()

# Remove empty subplots
for j in range(i + 1, nrows * ncols):
    fig.delaxes(axs[j // ncols, j % ncols])

plt.tight_layout()

plt.show()
