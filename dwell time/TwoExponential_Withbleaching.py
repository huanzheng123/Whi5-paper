import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import expon
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from matplotlib.backendsu.backend_pdf import PdfPages
import os

# Error handling for threadpoolctl issue
os.environ["OMP_NUM_THREADS"] = "1"

def neg_log_likelihood_single(params, data):
    lambda1 = params[0]
    likelihoods = expon.pdf(data, scale=1 / lambda1)
    neg_log_lik = -np.sum(np.log(likelihoods))
    return neg_log_lik

def neg_log_likelihood_with_bleaching(params, data, lambda_bleach):
    lambda1, lambda2, p1 = params
    p2 = 1 - p1

    # Combined decay rates including photobleaching
    lambda1_total = lambda1
    lambda2_total = lambda2

    # Adjusted likelihoods for the mixture of two exponentials with photobleaching
    likelihoods = (p1 * expon.pdf(data, scale=1 / lambda1_total) +
                   p2 * expon.pdf(data, scale=1 / lambda2_total))

    neg_log_lik = -np.sum(np.log(likelihoods))
    return neg_log_lik

def get_initial_parameters(data):
    data_reshaped = data.reshape(-1, 1)
    kmeans = KMeans(n_clusters=2, n_init=10, random_state=0).fit(data_reshaped)
    labels = kmeans.labels_
    cluster_centers = kmeans.cluster_centers_.flatten()
    lambda1_initial = 1 / cluster_centers[0]
    lambda2_initial = 1 / cluster_centers[1]
    p1_initial = np.sum(labels == 0) / len(data)
    return lambda1_initial, lambda2_initial, p1_initial

def run_single_exponential_model(data):
    initial_params = [1 / np.mean(data)]
    result = minimize(neg_log_likelihood_single, initial_params, args=(data,), bounds=[(1e-6, None)])
    lambda1 = result.x[0]
    neg_log_lik = neg_log_likelihood_single([lambda1], data)
    n = len(data)
    bic = 2 * neg_log_lik + np.log(n) * 1  # 1 parameter
    mean_time_single = 1 / lambda1
    return lambda1, mean_time_single, bic

def run_mixture_model(data_path, lambda_bleach, max_iterations=100, tolerance=0.05):
    data = pd.read_csv(data_path)['Bound Time'].values/10

    for i in range(max_iterations):
        lambda1_initial, lambda2_initial, p1_initial = get_initial_parameters(data)
        initial_params = [lambda1_initial, lambda2_initial, p1_initial]

        result = minimize(neg_log_likelihood_with_bleaching, initial_params,
                          args=(data, lambda_bleach), bounds=[(1e-6, None), (1e-6, None), (0, 1)])

        lambda1, lambda2, p1 = result.x
        p2 = 1 - p1

        # Ensure lambda1 is the larger (lower rate, higher mean) and lambda2 is the smaller (higher rate, lower mean)
        if lambda1 < lambda2:
            lambda1, lambda2 = lambda2, lambda1
            p1, p2 = p2, p1

        # Corrected mean times excluding photobleaching effect
        mean_time1_corrected = (1 / lambda1) * (1/lambda_bleach) /((1/lambda_bleach)-(1/lambda1))
        mean_time2_corrected = (1 / lambda2) * (1/lambda_bleach) /((1/lambda_bleach)-(1/lambda2))

        # Non-corrected mean times including photobleaching effect
        mean_time1_non_corrected = 1 / lambda1
        mean_time2_non_corrected = 1 / lambda2

    neg_log_lik = neg_log_likelihood_with_bleaching([lambda1, lambda2, p1], data, lambda_bleach)
    n = len(data)
    bic = 2 * neg_log_lik + np.log(n) * 3  # 3 parameters

    # Calculate weight fractions

    weight_fraction1 = (mean_time1_non_corrected * p1) / (mean_time1_non_corrected * p1 + mean_time2_non_corrected * p2)
    weight_fraction2 = (mean_time2_non_corrected * p2) / (mean_time1_non_corrected * p1 + mean_time2_non_corrected * p2)

    results = {
        "Estimated Mean Time 1 (corrected)": mean_time1_corrected,
        "Estimated Decay Rate 1 (corrected)": 1/mean_time1_corrected,
        "Estimated Mean Time 2 (corrected)": mean_time2_corrected,
        "Estimated Decay Rate 2 (corrected)": 1/mean_time2_corrected,
        "Estimated Mean Time 1 (non-corrected)": mean_time1_non_corrected,
        "Estimated Decay Rate 1 (non-corrected)": lambda1,
        "Estimated Mean Time 2 (non-corrected)": mean_time2_non_corrected,
        "Estimated Decay Rate 2 (non-corrected)": lambda2,
        "Fraction of Exponential Component 1": p1,
        "Fraction of Exponential Component 2": p2,
        "Weight Fraction 1": weight_fraction1,
        "Weight Fraction 2": weight_fraction2,
        "BIC Double Exponential": bic
    }
    return results, data, lambda1, lambda2, p1, p2

def plot_results(data, lambda1, lambda2, p1, p2, lambda1_single=None, mean_time_single=None):
    bins = np.linspace(0, data.max(), 50)
    hist, bins = np.histogram(data, bins=bins, density=True)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])

    plt.figure(figsize=(10, 6))

    # Plot histogram
    plt.hist(data, bins=bins, density=True, alpha=0.6, color='grey', label='Data', edgecolor='black')

    if lambda1_single is not None and mean_time_single is not None:
        # Plot single exponential fit
        pdf_single = expon.pdf(bin_centers, scale=1 / lambda1_single)
        plt.plot(bin_centers, pdf_single, 'g--', linewidth=2, label=f'Single Exponential Fit (mean={1/lambda1_single:.3f} λ={lambda1_single:.3f})')
    else:
        # Plot double exponential fit
        pdf1 = p1 * expon.pdf(bin_centers, scale=1 / lambda1)
        pdf2 = p2 * expon.pdf(bin_centers, scale=1 / lambda2)
        plt.fill_between(bin_centers, pdf2, pdf1 + pdf2, color='red', alpha=0.5, label=f'Exp 1 (mean={1/lambda1:.3f}, λ={lambda1:.3f})')
        plt.fill_between(bin_centers, 0, pdf2, color='blue', alpha=0.5, label=f'Exp 2 (mean={1/lambda2:.3f}, λ={lambda2:.3f})')
        plt.plot(bin_centers, pdf1 + pdf2, 'k-', linewidth=2, label='Total Fit (Double Exponential)')

    plt.xlabel('Bound Time')
    plt.ylabel('Density')
    plt.legend()
    plt.title('Histogram with Exponential Fits')

    plt.tight_layout()

def plot_results_with_trendlines(data, lambda1, lambda2, p1, p2, lambda1_single=None, mean_time_single=None):
    bins = np.linspace(0, data.max(), 50)
    hist, bins = np.histogram(data, bins=bins, density=True)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])

    plt.figure(figsize=(10, 6))

    # Plot histogram
    plt.hist(data, bins=bins, density=True, alpha=0.6, color='grey', label='Data', edgecolor='black')

    # Plot double exponential fit as trendlines
    pdf1 = p1 * expon.pdf(bin_centers, scale=1 / lambda1)
    pdf2 = p2 * expon.pdf(bin_centers, scale=1 / lambda2)

    plt.plot(bin_centers, pdf1, 'r-', linewidth=2, label=f'Exp 1 Trendline (mean={1/lambda1:.3f}, λ={lambda1:.3f})')
    plt.plot(bin_centers, pdf2, 'b-', linewidth=2, label=f'Exp 2 Trendline (mean={1/lambda2:.3f}, λ={lambda2:.3f})')
    #mean_time1_correctedplt.plot(bin_centers, pdf1 + pdf2, 'k-', linewidth=2, label='Total Fit (Double Exponential)')
    '''
    # Optionally, plot single exponential fit
    if lambda1_single is not None and mean_time_single is not None:
        pdf_single = expon.pdf(bin_centers, scale=1 / lambda1_single)
        plt.plot(bin_centers, pdf_single, 'g--', linewidth=2, label=f'Single Exponential Fit (mean= {1/lambda1_single:.3f}, λ={lambda1_single:.3f})')
    '''
    plt.xlabel('Diffusion Time')
    plt.ylabel('Density')
    plt.legend()
    plt.title('Histogram with Double Exponential Trendlines')

    plt.tight_layout()

# Known photobleaching rate
lambda_bleach = 1 / 10

# Run the mixture model analysis
data_path = 'C:\\Users\\JpRas\\OneDrive\\Escritorio\\SearchAnalysis\\MutA7\\timelapse\\AnalysisRebindCBC_start0_Quality5\\_ColBD_LIFE_FInalestM_withset3\\individual\pooled\\_ColBD_LIFE_rebind-strict-boundtime.csv'
results, data, lambda1, lambda2, p1, p2 = run_mixture_model(data_path, lambda_bleach)
lambda1_single, mean_time_single, bic_single = run_single_exponential_model(data)
results["BIC Single Exponential"] = bic_single
results["Estimated Mean Time (Single Exponential)"] = mean_time_single
results["Estimated Decay Rate (Single Exponential)"] = lambda1_single

# Calculate weight fraction for single exponential
weight_fraction_single = mean_time_single / (mean_time_single * len(data))
results["Weight Fraction (Single Exponential)"] = weight_fraction_single

# Save results and plots to a single PDF
output_folder = os.path.dirname(data_path)
pdf_path = os.path.join(output_folder, 'Expontential_fit_diffusion_B.pdf')

with PdfPages(pdf_path) as pdf:
    # Save text results as a figure
    plt.figure(figsize=(8, 4))
    plt.text(0.1, 0.05, '\n'.join([f"{key}: {value}" for key, value in results.items()]), fontsize=12)
    plt.axis('off')
    pdf.savefig()
    plt.close()

    # Save single exponential plot
    if bic_single < results["BIC Double Exponential"]:
        plot_results(data, lambda1, lambda2, p1, p2, lambda1_single, mean_time_single)
        pdf.savefig()
        plt.close()

    # Save double exponential plot
    plot_results(data, lambda1, lambda2, p1, p2)
    pdf.savefig()
    plt.close()

    # Save trendlines plot
    plot_results_with_trendlines(data, lambda1, lambda2, p1, p2, lambda1_single, mean_time_single)
    pdf.savefig()
    plt.close()

print(f"Results and plots saved to {pdf_path}")
