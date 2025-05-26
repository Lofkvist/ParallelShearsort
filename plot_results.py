import pandas as pd
import matplotlib.pyplot as plt

def plot_strong_scalability(filename):
    df = pd.read_csv(filename)
    # Baseline time = time at 1 core
    t1 = df.loc[df['Cores'] == 1, 'ExecutionTime(seconds)'].values[0]

    df['Speedup'] = t1 / df['ExecutionTime(seconds)']
    df['Efficiency'] = df['Speedup'] / df['Cores']

    fig, ax1 = plt.subplots()

    ax1.plot(df['Cores'], df['Speedup'], marker='o', label='Speedup')
    ax1.set_xlabel('Number of Cores')
    ax1.set_ylabel('Speedup', color='tab:blue')
    ax1.tick_params(axis='y', labelcolor='tab:blue')
    ax1.set_xticks(df['Cores'])

    ax2 = ax1.twinx()
    ax2.plot(df['Cores'], df['Efficiency'], marker='x', color='tab:red', label='Efficiency')
    ax2.set_ylabel('Efficiency', color='tab:red')
    ax2.tick_params(axis='y', labelcolor='tab:red')

    plt.title('Strong Scalability')
    fig.tight_layout()
    plt.show()

def plot_weak_scalability(filename):
    df = pd.read_csv(filename)
    # Baseline time at 1 core
    t1 = df.loc[df['Cores'] == 1, 'ExecutionTime(seconds)'].values[0]

    # Speedup = t1 / current time
    df['Speedup'] = t1 / df['ExecutionTime(seconds)']
    # Efficiency = Speedup / sqrt(number_of_cores) or just speedup / cores (common varies)
    # But weak scaling efficiency is often just (time at 1 core) / (time at n cores)
    df['Efficiency'] = df['Speedup'] / df['Cores']

    fig, ax1 = plt.subplots()

    ax1.plot(df['Cores'], df['Speedup'], marker='o', label='Speedup')
    ax1.set_xlabel('Number of Cores')
    ax1.set_ylabel('Speedup', color='tab:blue')
    ax1.tick_params(axis='y', labelcolor='tab:blue')
    ax1.set_xticks(df['Cores'])

    ax2 = ax1.twinx()
    ax2.plot(df['Cores'], df['Efficiency'], marker='x', color='tab:red', label='Efficiency')
    ax2.set_ylabel('Efficiency', color='tab:red')
    ax2.tick_params(axis='y', labelcolor='tab:red')

    plt.title('Weak Scalability')
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    strong_file = 'strong_scalability_results.txt'
    weak_file = 'weak_scalability_results.txt'

    print("Plotting strong scalability...")
    plot_strong_scalability(strong_file)

    print("Plotting weak scalability...")
    plot_weak_scalability(weak_file)
