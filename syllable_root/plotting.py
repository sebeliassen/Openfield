import matplotlib.pyplot as plt
import os
from data_processing import process_snippets
import numpy as np

def plot_syllable_signal(syllable, signal_type, syllable_snippets, save_figs=False, output_dir=''):
    """
    Plot average signals with SEM shading for a given syllable and signal type.

    Parameters
    ----------
    syllable : int
        The syllable index to plot.
    signal_type : str
        Type of photometry signal ('DS' or 'VS').
    syllable_snippets : defaultdict
        Nested dictionary containing signal snippets.
    save_figs : bool, optional
        Whether to save the figures as PNG files. Default is False.
    output_dir : str, optional
        Directory to save figures if save_figs is True. Default is '' (current directory).

    Returns
    -------
    None
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    fig.suptitle(f'Syllable {syllable} - {signal_type}_470 Signal')

    genotypes = ['WT', 'IFxDN HET']
    injection_phases = ['pre', 'post']

    # Define time axis based on snippet length
    time_axis = np.arange(-150, 151)  # Default values
    if syllable_snippets:
        any_genotype = next(iter(syllable_snippets))
        any_phase = next(iter(syllable_snippets[any_genotype]))
        any_syllable_snippets = syllable_snippets[any_genotype][any_phase].get(syllable, np.array([]))

        # Check if any_syllable_snippets is a 2D array and not empty
        if any_syllable_snippets.ndim == 2 and any_syllable_snippets.size > 0:
            snippet_length = any_syllable_snippets.shape[1]  # Use shape[1] since it's a 2D array
            m = (snippet_length - 1) // 2
            n = (snippet_length - 1) - m
            time_axis = np.arange(-m, n + 1)

    for i, genotype in enumerate(genotypes):
        for j, injection_phase in enumerate(injection_phases):
            snippets = syllable_snippets[genotype][injection_phase].get(syllable, np.array([]))

            # Ensure snippets is a 2D array and not empty
            if snippets.ndim == 2 and snippets.size > 0:
                mean, ci, n = process_snippets(snippets)
                ax = axes[i, j]
                if mean is not None:
                    color = 'blue' if signal_type == 'DS' else 'darkviolet'
                    ax.plot(time_axis, mean, label=f'{signal_type}_470', color=color)
                    ax.fill_between(time_axis, mean - ci, mean + ci, color=color, alpha=0.3)
                    ax.axvline(0, color='gray', linestyle='--')
                    ax.set_title(f'{genotype}, {injection_phase} injection (n={n})')
                    ax.set_ylabel('Signal (Normalized)')
                    if i == 1:
                        ax.set_xlabel('Frames Relative to Initiation')
                else:
                    ax.set_visible(False)
            else:
                axes[i, j].set_visible(False)  # Hide axes if there are no snippets

    # Adjust the layout
    plt.tight_layout()
    fig.subplots_adjust(top=0.88)

    if save_figs:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        # Construct a unique filename based on syllable and signal type
        filename = f'syllable_{syllable}_{signal_type}.png'
        filepath = os.path.join(output_dir, filename)
        plt.savefig(filepath)
        plt.close()  # Close the figure to free memory
    else:
        plt.show()


def mse_syllable_snippet(syllable_snippets, genotype, injection_phase, syllable):
    """
    Calculate the Mean Standard Error (MSE) for a specified genotype, injection phase, and syllable
    within a nested dictionary of syllable snippets.

    Parameters
    ----------
    syllable_snippets : defaultdict
        Nested dictionary containing signal snippets.
    genotype : str
        The genotype to access within the dictionary.
    injection_phase : str
        The injection phase to access within the dictionary.
    syllable : int
        The syllable index to access within the dictionary.

    Returns
    -------
    mse : float
        Mean Standard Error of the signal snippets, representing overall certainty.
    """
    # Access the array for the specific genotype, injection phase, and syllable
    snippets = syllable_snippets.get(genotype, {}).get(injection_phase, {}).get(syllable, np.array([]))

    # Ensure snippets is a 2D array and has data
    if snippets.ndim != 2 or snippets.size == 0:
        return None  # Return None if there are no snippets or they are not in the expected 2D format

    # Calculate standard deviation across snippets at each time point
    std_dev = np.std(snippets, axis=0)
    
    # Calculate standard error at each time point
    sem = std_dev / np.sqrt(snippets.shape[0])
    
    # Mean Standard Error (MSE) across all time points
    mse = np.mean(sem)

    return mse