# DarkDivNiches - Dark Diversity with Niche Distance Weighting

This is a fork of the original [DarkDiv package](10.1111/geb.13203), published by Carbona and Pärtel (2021), that extends the functionality to allow niche weighting of co-occurrence patterns using a distance matrix (e.g. based on Ellenberg indicator values or other ecological traits).

The rationale for adding this is that species that have strong abiotic relationships, e.g. very wet loving species, may be detected in plots with generalist species, i.e. species that occur in a wide range of conditions. However, this co-occurrence means that the strength of the strong abiotic relationship may be diluted by its co-occurrence with a generalist species. By weighting the co-occurrence patterns by niche similarity, we can give more weight to co-occurrences between species with similar niches, and less weight to co-occurrences between species with dissimilar niches.

This hopefully retains the strengths of using co-occurrence patterns to add additional biotic information, while reducing the dilution of strong abiotic relationships by generalist species.

## New Feature: Niche Weighting Support

The enhanced DarkDivNiches package supports weighting co-occurrence patterns by niche similarity through a new `niche_weighting` parameter, allowing for more ecologically meaningful dark diversity estimates.

### Key Features

- **Niche similarity weighting**: Co-occurrences between ecologically similar species are given more importance
- **Flexible input**: Accepts any square distance/dissimilarity matrix with proper species names (e.g., Ellenberg indicator values, functional traits)
- **Method compatibility**: Works with all DarkDivNiches methods (Hypergeometric, RawBeals, Favorability, ThresholdBeals)
- **Backward compatibility**: All existing functionality preserved when `niche_weighting = NULL` (default)

### How It Works

1. **Distance → Similarity**: Converts distance matrix to similarity weights using exponential decay: `exp(-distance/scale)`
2. **Weighted Co-occurrences**: Multiplies raw co-occurrence counts by niche similarity weights
3. **Enhanced Indication**: Creates more ecologically meaningful indication values between species pairs
4. **Improved Predictions**: Dark diversity estimates reflect both co-occurrence patterns and niche relationships

## Usage

### Basic Usage (Original Functionality)

```r
library(DarkDivNiches)
library(vegan)
data(dune)

# Standard dark diversity estimation
result_basic <- DarkDivNiches(x = dune, method = "Hypergeometric")
```

### With Niche Weighting (New Feature)

```r
# Create or load your niche distance matrix
# (e.g., based on Ellenberg indicator values)
species_names <- colnames(dune)
niche_distances <- your_ellenberg_distance_matrix  # Your actual distance matrix

# Dark diversity with niche weighting
result_weighted <- DarkDivNiches(x = dune, 
                          method = "Hypergeometric", 
                          niche_weighting = niche_distances)

# Compare results
summary(result_basic$Dark)
summary(result_weighted$Dark)
```

### Supported Methods

All original DarkDivNiches methods work with niche weighting:

```r
# Hypergeometric (default)
dd_hyper <- DarkDivNiches(x = dune, method = "Hypergeometric", niche_weighting = distances)

# Raw Beals
dd_beals <- DarkDivNiches(x = dune, method = "RawBeals", niche_weighting = distances)

# Favorability
dd_favor <- DarkDivNiches(x = dune, method = "Favorability", niche_weighting = distances)

# Threshold Beals
dd_thresh <- DarkDivNiches(x = dune, method = "ThresholdBeals", niche_weighting = distances)
```

## Implementation Details

### Requirements for Niche Distance Matrix

- **Square matrix**: Same number of rows and columns
- **Species names**: Row and column names must match species names in your data
- **Distance interpretation**: Smaller values = more similar niches
- **Symmetric**: Matrix should be symmetric (distances from A to B = B to A)
- **Zero diagonal**: Distance from each species to itself should be 0

### Validation and Error Handling

The implementation includes comprehensive validation:

- Checks matrix dimensions and symmetry
- Validates species name matching
- Provides informative error messages
- Handles edge cases gracefully

## Benefits

- **Ecologically meaningful**: Co-occurrences between similar species weighted more heavily
- **Flexible**: Accepts any square distance/similarity matrix with proper species names
- **Robust**: Comprehensive validation and error handling
- **Consistent**: Same approach works across all DarkDivNiches methods
- **Tested**: Extensive test suite ensures reliability and backward compatibility

## Testing and Validation

The implementation has been thoroughly tested with:

- ✅ Multiple test cases using real Ellenberg ecological data
- ✅ Compatibility with all DarkDivNiches methods
- ✅ Proper input validation and error handling  
- ✅ Preservation of all basic DarkDivNiches properties
- ✅ Backward compatibility (all existing functionality unchanged)

## References

Original DarkDiv package: Carbona, M., & Pärtel, M. (2021). Dark diversity in drylands: A new tool for understanding biodiversity. Global Ecology and Biogeography, 30(8), 1581-1594. [DOI: 10.1111/geb.13203](https://doi.org/10.1111/geb.13203)

## Installation

Install the development version from GitHub with:

```r
devtools::install_github("davidyshen/DarkDiv-with-Ellenberg")
```
