# Data from: Partitioning the effects of plant diversity on ecosystem functions at different trophic levels

[https://doi.org/10.5061/dryad.zs7h44jbp](https://doi.org/10.5061/dryad.zs7h44jbp)

## Description of the data and file structure

This data is from the PaNDiv experiment. Methodological details can be found in the associated manuscript.

### Files and variables

#### File: Plot\_information.csv

**Description:** Treatment information for the plots

##### Variables

* Bock: unique Block number(1-4)
* Plot: unique plot number (1-336)
* Nitrogen: fertilization information (0 = no fertilizer, 1 = nitrogen added)
* Fungicide: pathogen exclusion treatment (0 = control, 1 = fungicide)
* composition: plant species composition index (string indicating species richness (first number), unique species combination (second number), and whether the comination of species was assembled from the fast (f), slow (s) or mixed (m) species pool). Monocultures have the species abbreviation to indicate composition.
* Species_diversity: number of species sown in a plot. (1, 4, 8 or 20)
* Sown sla: Mean SLA of sown species
* Sown mpd sla: Mean pairwise difference in SLA of sown species.
* Notes: Excluded monocultures for the analysis are flagged.

#### File: Species\_information.csv

**Description:** Information on plant species in PaNDiv

##### Variables

* Species: Species abbreviation
* Binomial\_ name:	Binomial species name
* Pool: Species pool to which the species belongs	(F: fast growing species, S: slow growing species)
* Notes: Species that were excluded for the analysis in the associated manuscript are flagged

#### File: Abundance.csv

**Description:** relative abundance (cover) of species per plot and timepoint

##### Variables

* Plot: unique plot number (1-336)
* Harvest: time point when the data was collected (J17 = June 2017, A17 = August 2017, J18 = June 2018, A18 = August 2018)
* Species: Species abbreviation (full name can be found in the species information data sheet)
* Abundance: relative abundance of the species (proportion of cover. visually estimated)

#### File: Biomass.csv

**Description:** plant biomass collected 5cm above ground in each plot and time point

##### Variables

* Plot: unique plot number (1-336)
* Harvest: time point when the data was collected (J17 = June 2017, A17 = August 2017, J18 = June 2018, A18 = August 2018)
* Species: Species abbreviation (full name can be found in the species information data sheet)

- Biomass: Species biomass calculated based on plot biomass and relative abundance [g/m^2]
- Notes: notes

#### File: Pathogen\_infection.csv

**Description:** Pathogen infection measured as proportion of individuals with disease symptoms per species, plot and timepoitn

##### Variables

* Plot: unique plot number (1-336)
* Harvest: time point when the data was collected (J17 = June 2017, A17 = August 2017, J18 = June 2018, A18 = August 2018)
* Species: Species abbreviation (full name can be found in the species information data sheet)

- Infection: Proportion of individuals with diseases symptoms
- Notes: notes

#### File: Insect\_herbivory.csv

**Description:**  Herbivory measured as proportion of individuals with herbivory symptoms per species, plot and timepoitn

##### Variables

* Plot: unique plot number (1-336)
* Harvest: time point when the data was collected (J17 = June 2017, A17 = August 2017, J18 = June 2018, A18 = August 2018)
* Species: Species abbreviation (full name can be found in the species information data sheet)

- Herbivory: Proportion of individuals with herbivory symptoms
- Notes: notes

## Access information

Other publicly accessible locations of the data:

* NA

Data was derived from the following sources:

* own data

