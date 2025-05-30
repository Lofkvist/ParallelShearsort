import pandas as pd
import matplotlib.pyplot as plt

def plot_strong_efficiency(filename, ax):
    df = pd.read_csv(filename)
    df = df.sort_values('proc')
    
    # Get unique problem sizes
    problem_sizes = sorted(df['N'].unique())
    
    # Plot each problem size separately
    for N in problem_sizes:
        subset = df[df['N'] == N].copy()
        subset = subset.sort_values('proc')
        
        # Calculate efficiency: speedup / number of processes
        t1 = subset.loc[subset['proc'] == 1, 'time'].values[0]
        subset['speedup'] = t1 / subset['time']
        subset['efficiency'] = subset['speedup'] / subset['proc']
        
        # Plot this problem size
        ax.plot(subset['proc'], subset['efficiency'], marker='^', 
                label=f'N={N}', linewidth=2)
    
    # Add ideal efficiency line (should be 1.0)
    ax.axhline(y=1.0, linestyle='--', color='gray', alpha=0.7, label='Ideal Efficiency')
    
    ax.set_xlabel('Number of Processes', color='black')
    ax.set_ylabel('Efficiency', color='black')
    ax.set_title('Strong Efficiency', color='black', fontsize=14)
    ax.set_ylim(0, 1.1)
    ax.tick_params(axis='x', labelcolor='black')
    ax.tick_params(axis='y', labelcolor='black')
    ax.set_xticks(sorted(df['proc'].unique()))
    ax.legend()
    ax.grid(True, alpha=0.3)

if __name__ == "__main__":
    strong_file = 'strong_scal.txt'
    
    # Create figure for strong efficiency plot only
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_strong_efficiency(strong_file, ax)
    
    plt.show()