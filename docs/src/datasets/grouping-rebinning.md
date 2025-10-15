# Grouping and Rebinning Spectral Data

When working with spectral data, it's often necessary to group or rebin channels to improve signal-to-noise ratios or reduce the number of data points. SpectralFitting.jl provides several utilities for this purpose.

## Overview

Grouping spectral channels is a common operation in X-ray spectral analysis. The grouping is controlled by a `grouping` vector where:
- `1` (or `NEW_GRP`) indicates the next channel (if it exists) will starts a new group
- `0` (or `CONTINUE_GRP`) indicates the channel continues the previous group

## Grouping by Minimum Counts

You can group channels to achieve a minimum number of counts per bin:

```julia
using SpectralFitting

# Load your data
data = OGIPDataset("spectrum.pha")

# Group to minimum 25 counts per bin
regroup!(data; min_counts = 25)
```

This is useful when you want to ensure sufficient statistics in each bin for χ² fitting.

## Grouping by Signal-to-Noise Ratio

New in this version, you can group channels to achieve a minimum signal-to-noise ratio per bin:

```julia
# Group to minimum S/N = 5 per bin
regroup!(data; min_snr = 5.0)
```

The signal-to-noise ratio is calculated following the [XMM SAS specgroup](https://xmm-tools.cosmos.esa.int/external/sas/current/doc/specgroup/node10.html) convention:

### Without Background

For spectra without background subtraction:

```math
\text{S/N} = \frac{S}{\sqrt{S}}
```

where $S$ is the total source counts in the accumulated bin.

### With Background

For background-subtracted spectra:

```math
\text{S/N} = \frac{S - B \times \text{areanorm}}{\sqrt{S + B \times \text{areanorm}^2}}
```

where:
- $S$ is the total source counts
- $B$ is the total background counts
- $\text{areanorm} = \frac{\text{BACKSCAL}_{\text{src}}}{\text{BACKSCAL}_{\text{bkg}}} \times \frac{\text{EXPOSURE}_{\text{src}}}{\text{EXPOSURE}_{\text{bkg}}}$

The background is automatically used if present in the dataset.

### Units

The S/N grouping automatically handles both count and count rate units:
- `u"counts"`: Uses counts directly
- `u"counts / s"`: Converts to counts using exposure time

## Low-Level Interface

For more control, you can use the low-level functions directly on `Spectrum` objects:

```julia
# Group by minimum counts
SpectralFitting.group_min_counts!(spectrum, 25)

# Group by minimum S/N without background
SpectralFitting.group_min_snr!(spectrum, 5.0)

# Group by minimum S/N with background
SpectralFitting.group_min_snr!(spectrum, 5.0; background = background_spectrum)
```

## Example: Comparing Grouping Methods

```julia
using SpectralFitting, Plots

# Load data
data = OGIPDataset("spectrum.pha")

# Create copies for different grouping strategies
data_counts = deepcopy(data)
data_snr = deepcopy(data)

# Apply different grouping
regroup!(data_counts; min_counts = 25)
regroup!(data_snr; min_snr = 3.0)

# Compare the number of bins
println("Original bins: ", length(data.spectrum.data))
println("After min_counts=25: ", count(==(1), data_counts.spectrum.grouping))
println("After min_snr=3: ", count(==(1), data_snr.spectrum.grouping))

# Visualize
p1 = plot(data_counts, title="Grouped by counts ≥ 25")
p2 = plot(data_snr, title="Grouped by S/N ≥ 3")
plot(p1, p2, layout=(2,1))
```

## Choosing Between Methods

- **Minimum counts**: Use when you need a specific minimum number of counts for χ² statistics (typically 15-25 counts per bin)
- **Minimum S/N**: Use when you want consistent statistical significance across your spectrum, especially useful when background varies significantly

## API Reference

```@docs
SpectralFitting.regroup!
SpectralFitting.group_min_counts!
SpectralFitting.group_min_snr!
```
