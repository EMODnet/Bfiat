# Package index

## Density

Functions to estimate density under fishing

- [`density_perturb()`](https://EMODnet.github.io/Bfiat/reference/density.md)
  [`density_metier()`](https://EMODnet.github.io/Bfiat/reference/density.md)
  [`density_logistic()`](https://EMODnet.github.io/Bfiat/reference/density.md)
  :

  Estimates density at specific times for the logistic and perturbation
  model. `density_metier` allows to specify different metiers.

## Steady-state

Functions to estimate steady-state under fishing

- [`steady_perturb()`](https://EMODnet.github.io/Bfiat/reference/steadydensity.md)
  [`steady_metier()`](https://EMODnet.github.io/Bfiat/reference/steadydensity.md)
  [`steady_logistic()`](https://EMODnet.github.io/Bfiat/reference/steadydensity.md)
  : Estimate steady-state density for logistic and perturbation model
  for one and more metiers.

## Trawl

Functions to estimate density at trawling events

- [`trawl_perturb()`](https://EMODnet.github.io/Bfiat/reference/trawlingdensity.md)
  [`trawl_metier()`](https://EMODnet.github.io/Bfiat/reference/trawlingdensity.md)
  : Estimates density before and after a trawling event, or the average
  density between trawls.

## Models as in deSolve

Functions to estimate density using deSolve input

- [`run_perturb()`](https://EMODnet.github.io/Bfiat/reference/trawlingfun.md)
  [`run_logistic()`](https://EMODnet.github.io/Bfiat/reference/trawlingfun.md)
  : Functions to estimate the impact of bottom disturbances on benthic
  taxa. Allows for non-regularly spaced events.

## Model parameters

Functions to estimate the model parameters

- [`par_r()`](https://EMODnet.github.io/Bfiat/reference/trawlingpars.md)
  [`par_d()`](https://EMODnet.github.io/Bfiat/reference/trawlingpars.md)
  [`par_K()`](https://EMODnet.github.io/Bfiat/reference/trawlingpars.md)
  [`par_m()`](https://EMODnet.github.io/Bfiat/reference/trawlingpars.md)
  : Estimate model parameters for use in the fishing impact model.

## Ecosystem functions

Functions to estimate ecosystem functions or traits

- [`get_Db_model()`](https://EMODnet.github.io/Bfiat/reference/getDbModel.md)
  [`get_irr_model()`](https://EMODnet.github.io/Bfiat/reference/getDbModel.md)
  [`get_sptrait_model()`](https://EMODnet.github.io/Bfiat/reference/getDbModel.md)
  : From model output (time x taxon density) to time x
  species_bioturbation or time x species_bioirrigation potentials.
  get_sptrait_model extracts the contribution of a species to an
  arbitrary trait
- [`get_trait_model()`](https://EMODnet.github.io/Bfiat/reference/getTraitModel.md)
  : Function to obtain (time x trait) data by combining model output
  (time x taxon) with trait (taxon x trait) data.

## Sensitivity functions

Sensitivity of an area or taxon to fishing parameters

- [`sensitivity_taxon()`](https://EMODnet.github.io/Bfiat/reference/sensitivity.md)
  [`sensitivity_area()`](https://EMODnet.github.io/Bfiat/reference/sensitivity.md)
  : Sensitivity of an area or taxon to fishing parameters

## Critical values

Critical values for fishing or taxon parameters

- [`critical_sar()`](https://EMODnet.github.io/Bfiat/reference/critical.md)
  [`critical_gpd()`](https://EMODnet.github.io/Bfiat/reference/critical.md)
  [`critical_r()`](https://EMODnet.github.io/Bfiat/reference/critical.md)
  [`critical_d()`](https://EMODnet.github.io/Bfiat/reference/critical.md)
  : Critical values for fishing or taxon parameters.

## Data sets

- [`SAR`](https://EMODnet.github.io/Bfiat/reference/SARdata.md) :
  Swept-area ratios in the Northsea, OSPAR data.
- [`NIOZdata`](https://EMODnet.github.io/Bfiat/reference/BFIATdata.md) :
  Biological data from a transect across the North sea for running the
  BFIAT tool.

- [`Bfiat-package`](https://EMODnet.github.io/Bfiat/reference/BFIAT_package.md)
  [`Bfiat`](https://EMODnet.github.io/Bfiat/reference/BFIAT_package.md)
  : Bottom Fishing Impact Assessment Tools
