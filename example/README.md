This folder contains some example datasets to test CoaTran:
* `small_transmissions.tsv` is a small transmission network with a single seed
  * Can be used with `small_times.tsv`, `small_times_multisample.tsv`, or `small_times_subsample.tsv`
* `small_transmissions_multiseed.tsv` is a small transmission network with multiple seeds
  * Can be used with `small_times.tsv`, `small_times_multisample.tsv`, or `small_times_subsample.tsv`
* `small_times.tsv` is a small sample times file with a single sample time per individual
  * Can be used with `small_transmissions.tsv` or `small_transmissions_multiseed.tsv`
* `small_times_multisample.tsv` is a small sample times file with multiple sample times per individual
  * Can be used with `small_transmissions.tsv` or `small_transmissions_multiseed.tsv`
* `small_times_subsample.tsv` is a small sample times file in which a seed's entire transmission chain was not sampled
  * Should only be used with `small_transmissions_multiseed.tsv` (it's not interesting with `small_transmissions.tsv`)
* `big_transmissions.tsv` is a large transmission network with a single seed
  * Can only be used with `big_times.tsv`
* `big_transmissions_multiseed.tsv` is a large transmission network with multiple seeds
  * Can only be used with `big_times_multiseed.tsv`
