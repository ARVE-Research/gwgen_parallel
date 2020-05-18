This version of the globally parameterized, globally applicable, monthly-to-daily weather generator GWGEN is described in the following publication (freely available for download [here](https://www.geosci-model-dev.net/10/3771/2017/)):

Sommer, P. S., & Kaplan, J. O. (2017). A globally calibrated scheme for generating daily meteorology from monthly statistics: Global-WGEN (GWGEN) v1.0. *Geoscientific Model Development*, 10(10), 3771-3791. doi:10.5194/gmd-10-3771-2017

GWGEN-grid is designed to work with netCDF formatted gridded input and output, and depends on an installation of [NetCDF-Fortran](https://www.unidata.ucar.edu/downloads/netcdf/index.jsp) (version 4.1-beta2 or greater).

For new clones, run `./autogen.sh` to generate the `configure` script. Then `./configure --help` to check the configuration options, and finally: `./configure [options] ; make` to compile.
