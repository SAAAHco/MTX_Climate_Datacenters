# Run the complete pipeline in order. Set MADRID_DATA and TEXAS_DATA (or place the data files in data/) first.
import subprocess
import sys
import time

STEPS = ['loaddata', 'run_base', 'analysis', 'sens', 'gsa', 'extras', 'figs', 'fig2', 'propcheck']

if __name__ == '__main__':
    steps = sys.argv[1:] or STEPS
    for step in steps:
        t0 = time.time()
        print(f'--- {step}.py', flush=True)
        subprocess.run([sys.executable, f'{step}.py'], check=True)
        print(f'--- {step}.py finished in {time.time() - t0:.0f} s', flush=True)
