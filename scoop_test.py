import socket
import scoop.futures
from multiprocessing import cpu_count
import time
import numpy as np

M = np.random.rand(2000, 2000)
n = list(range(1000))

def print_procs(x):
    # print(socket.gethostname())
    # print('num_cores:' + str(cpu_count()))
    return np.sum(x * M * M)

def main():
    start_time = time.perf_counter()
    out = list(scoop.futures.map(print_procs, n))
    print(time.perf_counter() - start_time)

if __name__ == "__main__":
    main()
