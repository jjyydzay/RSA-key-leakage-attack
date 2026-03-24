from multiprocessing import Pool, cpu_count
import time

def square(n):
    time.sleep(1)  # Simulate a time-consuming task
    return n * n

def execute_task_parallel():
    with Pool(5) as p:
        return p.map(square, range(1, 6))
    
def execute_task_sequential():
    return [square(i) for i in range(1, 6)]

def measure_time(execution_function, execution_type):
    start_time = time.time()
    results = execution_function()
    end_time = time.time()
    
    print(f"{execution_type} execution:")
    for i, result in enumerate(results, start=1):
        print(f"Square of {i}: {result}")
    print(f"Total time: {end_time - start_time}")

def main():
    print("cpu number:", cpu_count())
    measure_time(execute_task_sequential, "Sequential")
    measure_time(execute_task_parallel, "Parallel")

if __name__ == "__main__":
    main()