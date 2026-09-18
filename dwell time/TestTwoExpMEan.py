import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import expon
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

# Error handling for threadpoolctl issue
os.environ["OMP_NUM_THREADS"] = "1"


def neg_log_likelihood_fixed_lambdas(params, data, lambda1, lambda2):
    p1 = params[0]
    p2 = 1 - p1

    likelihoods = (p1 * expon.pdf(data, scale=1 / lambda1) +
                   p2 * expon.pdf(data, scale=1 / lambda2))

    neg_log_lik = -np.sum(np.log(likelihoods))
    return neg_log_lik


def search_best_lambdas(data, initial_lambda1, initial_lambda2, range_factor=2.0, num_steps=20):
    best_neg_log_lik = np.inf
    best_lambdas = (initial_lambda1, initial_lambda2)
    best_p1 = None

    lambda1_range = np.linspace(initial_lambda1 * 0.5, initial_lambda1 * range_factor, num_steps)
    lambda2_range = np.linspace(initial_lambda2 * 0.5, initial_lambda2 * range_factor, num_steps)

    for lambda1 in lambda1_range:
        for lambda2 in lambda2_range:
            result = minimize(neg_log_likelihood_fixed_lambdas, [0.5],
                              args=(data, lambda1, lambda2), bounds=[(0, 1)])
            neg_log_lik = result.fun
            p1 = result.x[0]

            if neg_log_lik < best_neg_log_lik:
                best_neg_log_lik = neg_log_lik
                best_lambdas = (lambda1, lambda2)
                best_p1 = p1

    best_lambda1, best_lambda2 = best_lambdas
    return best_lambda1, best_lambda2, best_p1, best_neg_log_lik


def run_fixed_lambda_mixture_model(data_path, initial_lambda1, initial_lambda2):
    data = pd.read_csv(data_path)['Bound Time'].values / 10  # Convert to seconds by dividing by 10

    # Search for best lambdas within the specified range
    best_lambda1, best_lambda2, p1, neg_log_lik = search_best_lambdas(data, initial_lambda1, initial_lambda2)

    p2 = 1 - p1

    # Corrected mean times excluding photobleaching effect
    mean_time1 = 1 / best_lambda1
    mean_time2 = 1 / best_lambda2

    n = len(data)
    aic = 2 * 1 + 2 * neg_log_lik  # Only optimizing 1 parameter (p1)

    # Calculate weight fractions
    weight_fraction1 = (mean_time1 * p1) / (mean_time1 * p1 + mean_time2 * p2)
    weight_fraction2 = (mean_time2 * p2) / (mean_time1 * p1 + mean_time2 * p2)

    results = {
        "Estimated Mean Time 1": mean_time1,
        "Estimated Decay Rate 1": best_lambda1,
        "Estimated Mean Time 2": mean_time2,
        "Estimated Decay Rate 2": best_lambda2,
        "Fraction of Exponential Component 1": p1,
        "Fraction of Exponential Component 2": p2,
        "Weight Fraction 1": weight_fraction1,
        "Weight Fraction 2": weight_fraction2,
        "AIC": aic
    }
    return results, data, best_lambda1, best_lambda2, p1, p2


def plot_results(data, lambda1, lambda2, p1, p2):
    bins = np.linspace(0, data.max(), 50)
    hist, bins = np.histogram(data, bins=bins, density=True)
    bin_centers = (0.5 * (bins[1:] + bins[:-1]))

    n_samples = len(data)  # Get the number of data points

    with PdfPages(os.path.join(os.path.dirname(data_path), 'exponential_fits_best_lambdas_d12.pdf')) as pdf:
        plt.figure(figsize=(12, 6))

        # Plot histogram
        plt.hist(data, bins=bins, density=True, alpha=0.6, color='grey', label=f'Data (n={n_samples})',
                 edgecolor='black')

        # Plot double exponential fit
        pdf1 = p1 * expon.pdf(bin_centers, scale=1 / lambda1)
        pdf2 = p2 * expon.pdf(bin_centers, scale=1 / lambda2)
        plt.fill_between(bin_centers, pdf2, pdf1 + pdf2, color='red', alpha=0.5,
                         label=f'Exp 1 (Mean={1 / lambda1:.3f} λ={lambda1:.3f})')
        plt.fill_between(bin_centers, 0, pdf2, color='blue', alpha=0.5,
                         label=f'Exp 2 (Mean={1 / lambda2:.3f} λ={lambda2:.3f})')
        plt.plot(bin_centers, pdf1 + pdf2, 'k-', linewidth=2, label='Total Fit (Double Exponential)')

        plt.xlabel('Bound Time (seconds)')
        plt.ylabel('Probability Density Function')
        plt.legend()
        plt.title('Two Exponential Fit with Best Lambdas')

        plt.tight_layout()
        pdf.savefig()  # Save the current figure
        plt.close()

        # Create text page with results
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.axis('off')
        text = ""
        for key, value in results.items():
            text += f"{key}: {value}\n"
        ax.text(0.1, 0.5, text, transform=ax.transAxes, fontsize=12, verticalalignment='center',
                horizontalalignment='left', wrap=True)
        plt.title('Results Summary')
        pdf.savefig()  # Save the current figure
        plt.close()


# Initial fixed lambdas
initial_lambda1 = 1/1  # Example fixed value
initial_lambda2 = 1/4 # Example fixed value

# Run the mixture model analysis with searching for the best lambdas
data_path = 'C:\\Users\\JpRas\\OneDrive\\Escritorio\\SearchAnalysis\\wt\\timelapse\\set1\\_ColBD_LIFE_FInalestM\\_ColBD_LIFE_rebind-strict-boundtime.csv'
results, data, lambda1, lambda2, p1, p2 = run_fixed_lambda_mixture_model(data_path, initial_lambda1, initial_lambda2)

print("Results:")
for key, value in results.items():
    print(f"{key}: {value}")

plot_results(data, lambda1, lambda2, p1, p2)
