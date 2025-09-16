## RBS 

Package for [BEAST 2](http://beast2.org).

This package contains the reversible jump based substitution model (but bModelTest is recommended instead) and the AutoPartition model described in

NB: for the Free rate model, use the [FreeRateModel](https://github.com/BEAST2-Dev/FreeRateModel) package.

Evolutionary rates and HBV: issues of rate estimation with Bayesian molecular methods
R Bouckaert, MV Alvarado-Mora, JRR Pinho
Antiviral therapy 18 (3_part_2), 497-504

It also contains the variable selection based substitution  model (VS model) as described in the book:

Alexei J. Drummond and Remco R. Bouckaert. 
Bayesian evolutionary analysis with BEAST 2. 
CUP. 2014


## Installation

The latest version is for BEAST v2.7, which can be obtained from [here](http:/beast2.org).

Make sure `package-extra-2.7.xml` is in the list of package repositories:

* In BEAUti, select `File/Manage packages`
* Click `Package repositories` at the bottom of the screen
* A window pops up with the list of package repositories -- if `package-extra-2.7.xml` is not there click `Add URL`
* A window pops up where you can enter "https://raw.githubusercontent.com/CompEvol/CBAN/refs/heads/master/packages-extra-2.7.xml" (without quotes). Click OK
* Click the `Close button` to return to the package list.

The RBS package should now be listed in the packagemanager window, so select the row with RBA and click "Install/Upgrade" to install.
