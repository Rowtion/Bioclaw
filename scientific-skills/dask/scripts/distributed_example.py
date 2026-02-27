"""
Dask Distributed Example: Parallel Task Execution
Demonstrates using Dask Distributed for custom parallel workflows.
"""
from dask.distributed import Client, as_completed
import time


def setup_cluster():
    """Set up a local Dask cluster."""
    client = Client(n_workers=4, threads_per_worker=2)
    print(f"Dashboard: {client.dashboard_link}")
    return client


def process_data_chunk(data_chunk):
    """Simulate data processing."""
    time.sleep(1)  # Simulate work
    return sum(data_chunk)


def parallel_map_example(client):
    """Process data in parallel using client.map."""
    # Create data chunks
    data = [range(i, i + 100) for i in range(0, 1000, 100)]
    
    # Submit all tasks
    futures = client.map(process_data_chunk, data)
    
    # Gather results as they complete
    results = []
    for future in as_completed(futures):
        result = future.result()
        results.append(result)
        print(f"Completed: {result}")
    
    return results


def parameter_sweep(client):
    """Run parameter sweep in parallel."""
    def evaluate_model(alpha, beta):
        # Simulate model evaluation
        time.sleep(0.5)
        score = alpha * 2 + beta * 3
        return {'alpha': alpha, 'beta': beta, 'score': score}
    
    # Define parameter grid
    alphas = [0.1, 0.5, 1.0, 2.0]
    betas = [0.01, 0.1, 1.0]
    
    # Submit all combinations
    futures = []
    for alpha in alphas:
        for beta in betas:
            future = client.submit(evaluate_model, alpha, beta)
            futures.append(future)
    
    # Get all results
    results = client.gather(futures)
    best = max(results, key=lambda x: x['score'])
    
    print(f"Best parameters: {best}")
    return results


if __name__ == "__main__":
    client = setup_cluster()
    try:
        print("Running parallel map example...")
        results = parallel_map_example(client)
        print(f"Total sum: {sum(results)}")
        
        print("\nRunning parameter sweep...")
        sweep_results = parameter_sweep(client)
    finally:
        client.close()
