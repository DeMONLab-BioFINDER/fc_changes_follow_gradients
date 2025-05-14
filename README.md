Code repository for reproducing "Age and Alzheimer’s disease affect functional connectivity along separate axes of functional brain organization"

To install all dependencies, for _the whole_ project run 
`renv::restore()`

However, you only need the python dependencies if you want to extract timeseries from images.

Unfortunately, we cannot host the real data. So, to reproduce the paper using synthetic data, set `real_data` to `FALSE` in main.R

There are three other arguments that you may need to consider:

```
from_start <- FALSE
create_brain_permutations <- FALSE
extract_timeseries <- FALSE
```
If you run main.R for the first time you should likely set `from_start` to `TRUE`. `create_brain_permutations` and `extract_timeseriesa` can
both be `FALSE` as the synthetic data are timeseries already, and the brain permutations have already been created and are commited in the repo. 

The code is _not_ thoroughly tested and more sparsely documented than I would have preferred at this time. If you run into any problems, please create an issue and I will help to the best of my ability :)
