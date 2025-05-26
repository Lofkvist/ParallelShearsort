import pandas as pd
import matplotlib.pyplot as plt

def plot_strong_scalability(filename, ax):
    df = pd.read_csv(filename)
    df = df.sort_values('procs')

    t1 = df.loc[df['procs'] == 1, 'time'].values[0]
    df['speedup'] = t1 / df['time']

    ax.plot(df['procs'], df['speedup'], marker='o', label='Measured Speedup', color='tab:blue')
    ax.plot(df['procs'], df['procs'], linestyle='--', label='Ideal Speedup', color='gray')  # Ideal line

    ax.set_xlabel('Number of Processes', color='black')
    ax.set_ylabel('Speedup', color='black')
    ax.set_title('Strong Scalability', color='black')
    ax.tick_params(axis='x', labelcolor='black')
    ax.tick_params(axis='y', labelcolor='black')
    ax.set_xticks(df['procs'])
    ax.legend()

    plt.tight_layout()

def plot_weak_scalability(filename,ax):
    df = pd.read_csv(filename)
    df = df.sort_values('procs')

    t1 = df.loc[df['procs'] == 1, 'time'].values[0]
    df['efficiency'] = t1 / df['time']

    ax.plot(df['procs'], df['efficiency'], marker='x', color='tab:red', label='Efficiency')

    ax.set_xlabel('Number of Processes', color='black')
    ax.set_ylabel('Efficiency', color='black')
    ax.set_title('Weak Scalability', color='black')
    ax.set_ylim(0, 1.05)
    ax.tick_params(axis='x', labelcolor='black')
    ax.tick_params(axis='y', labelcolor='black')
    ax.set_xticks(df['procs'])
    ax.legend()

    plt.tight_layout()

if __name__ == "__main__":
    strong_file = 'strong_scalability_results.txt'
    weak_file = 'weak_scalability_results.txt'
    
    fig1, ax1 = plt.subplots()
    plot_strong_scalability(strong_file, ax1)

    fig2, ax2 = plt.subplots()
    plot_weak_scalability(weak_file,ax2)
    plt.show()
