# Brobustedger 1.1.0

## Major changes
- Replaced the `missForest`-based missing and outlier value imputation with
  the `missRanger` package to improve computational efficiency and scalability.

## Improvements
- Faster imputation for large datasets due to ranger-based implementation.
- Reduced memory usage compared to `missForest`.

## Internal changes
- Updated imputation functions to support `missRanger::missRanger()`.
