# DNA_metabarcoding
Clean(ish) repo of the DNA metabarcoding datasets

### Data
The data folder has two subfolders, one for data about the collection of the samples and one for the sequenced DNA results.

#### data/collection_data

files to start with:

* _krat_fecal_data.csv_: data about the collection of fecal samples from only the kangaroo rats. It is the same information that is in the fecal_sample_collection.csv but I've filtered based on sample_type (trap only), notes (nothing from the trap_bait_test), and species (DO and DM only)

additional files:

* _fecal_sample_collection_.csv: data for every single fecal sample
* _plant_voucher_collection_: data for every plant voucher collected
* _vial_id.csv_: dataset with the all of the vial ids (same as vial_barcode and SampleID in various dataframes) and whether they are plant samples or fecal samples

#### data/sequenced_data

files to start with:

* _krat_fecal_trnl_no_millet.csv_: cleaned DNA data from the trnL primer. This only includes  DNA from fecal samples from kangaroo rats. Any OTUs that have been identified as potentially being millet have been removed.
  + trnL is the primer that is not good at identifying grasses (Poaceae)
  + you can check the script called find_millet_OTUs.r to see how I decided which OTUs to remove. Basically, anything that was an unidentified Poaceae (family), unidenified Panicoideae (subfamily), or unidentified Panicum (genus) was removed
* _krat_fecal_its_no_millet.csv_: leaned DNA data from the ITS2 primer. This only includes  DNA from fecal samples from kangaroo rats. Any OTUs that have been identified as potentially being millet have been removed.
  + ITS2 is the primer that is much better at identifying grasses (Poaceae)
  + you can check the script called find_millet_OTUs.r to see how I decided which OTUs to remove. Basically, any Panicum (genus) that was either unidentified at the species level or was identified as a Panicum species not found at the site (i.e., Panicum hirticule) was removed

additional files:

* _krat_fecal_trnl.csv_: cleaned DNA data from the trnL primer. This only includes  DNA from fecal samples from kangaroo rats but potential millet OTUs have not been removed
* _krat_fecal_its.csv_: cleaned DNA data from the ITS2 primer. This only includes  DNA from fecal samples from kangaroo rats but potential millet OTUs have not been removed
* _trnL_reads_WeeTU.csv_: trnL reads data for all samples (fecal from all species and plant vouchers)
* _ITS2_reads_WeeTU.csv_: ITS2 reads data for all samples (fecal from all species and plant vouchers)

##### What's a WeeTU?
WeeTUs (or WTU) are the Weecology re-calculation of the OTUs (organismal taxonomic unit).
DNA sequences that are identified as the same species (or other taxonomic level) can have different OTUs based on minor changes in the nucleotides. 
At any given taxonomic level, all of the matching identities are given the same WTU number.

For example:

* Any OTU that is idenified down to the species _Muhlenbergia porteri_ is given the same WTU.species number: 613
* Any OTU identified to a different species of _Mulhenbergia_ is given a unique WTU.species. _Mulhenbergia richardsonis_ is WTU.species 614 while _Mulhenbergia arenacea_ is WTU.species 609.
* OTUs which are identified to the genus _Mulhenbergia_ but not identified at a species level (species = NA) are all given the same WTU.species number (618)
* Although all of the OTUs above will have different WTU.species numbers, they will all have the same WTU.genus number (and all other WTUs) because they are all identified to the same genus: WTU.genus = 395; WTU.family = 54.

If you decide to use WTUs, I would suggest using WTU.genus for most analyses. Most of the species which are identified do not match up to species we have at our site. 

### Scripts
I've uploaded a couple scripts, though they probably aren't too helpful at this point. I'm hoping to upload a script with some example of how to manipulate the data and use the functions I've written that might be useful.

* _find_millet_OTUs_: script for determining which OTUs/WTUs could be millet
* _functions_: script with some functions that might eventually be helpful
* _test_script_: a quick script to read in the original datasets, select kangaroo rat fecal samples and remove millet. I used this script to create the krat_* files above for you to use.