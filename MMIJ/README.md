File descriptions:

* `csv2nc.py`: script to transform the MMIJ `.csv` file to a netCDF4 file (most recommended)
* `csv2nc-compound-sem.py`: script to transform the MMIJ `.csv` file to a netCDF4 file (not recommended, uses compound data structures, for which software support is still lacking)
* `csv2nc-compound-nosem.py`: script, variant of `csv2nc-compound-sem.py`, to transform the MMIJ `.csv` file to a netCDF4 file (not recommended, uses compound data structures, for which software support is still lacking)
* `metadata-*.yaml`: metadata files used by the transformation scripts
* `structure.yaml`: structure file used by the transformation scripts
* `automated-check-snippets.py`: script for (semi-)automated quality checks, such as outlier detection
* `discrete_values_check-example.py`: script showing an example of a live session to check discrete values
* `get_uncertainty+bias.py`: script calculating average relative statistic uncertainties and average relative bias in the sample standard deviation
