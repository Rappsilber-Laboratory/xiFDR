xiFDR
=====

You can download the latest release of xiFDR **[here](https://www.rappsilberlab.org/software/xifdr/)**. The application can then be run by clicking on "startWindows", "startUnix" or "startMacOS". xiFDR is implemented as a java application and requires java 8 or above to run.

xiFDR is an application for estimating false discovery rates (FDRs) in crosslinking mass spectrometry. It filters crosslinked peptide spectra matches (CSMs) to a list of identifications and derives associated confidence values.

It performs a generic FDR calculations for CSMs and resulting peptide pairs, crosslinks and protein pairs. It complies with the standards for data reporting set by the HUPO proteomics standards initiative and can output results in .mzIdentML 1.3.0 format for deposition in databases. It is search engine-agnostic and can therefore perform FDR filtering on results from [xiSEARCH](https://github.com/Rappsilber-Laboratory/xisearch) but also other crosslinking MS search engines. The output can then be directly uploaded to [xiView.org](https://xiview.org/xiNET_website/index.php) for spectral analysis, network visualization and mapping to structures.


For questions regarding usage of xiFDR, please open a [discussion](https://github.com/Rappsilber-Laboratory/XiSearch/discussions).

### Background
Correct estimation of false discovery rates in crosslinked peptide identifications presents several quirks that must be taken into account. The 2 main issues are 1) correct estimation of FDR from the target-decoy approach at the level of CSMs 2) Correct handling of the random space for self and heteromeric crosslinks and 3) propagation of error from CSMs to peptide pairs, crosslinked residue pairs and protein-protein interactions. xiFDR handles both of these issues allowing for accurate FDR estimation from a search result file.

1. Estimation of FDRs in crosslinking MS via the target-decoy approach

The theoretical basis of adapting the proteomic target-decoy approach for crosslinking MS are covered in [Waltzoheni et al. 2012](https://www.nature.com/articles/nmeth.2103) and [Fisher et al. 2017](https://pubs.acs.org/doi/pdf/10.1021/acs.analchem.6b03745) and updated for heterobifunctional crosslinkers (crosslinkers with different reactivities at either end) in [Fisher et al. 2018](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0196672). The approach outlined in those papers is implemented in xiFDR. Briefly, decoy protein sequences are added to the database and spectra can then be matched to target linear peptides (T), decoy linear peptides (D) or to crosslinked peptide pairs made up of 2 target sequences (TT), a target-decoy (TD) sequence pair or a decoy-decoy pair (DD).

The formula to estimate FDR (provided the chance of matching decoys accurately models the chance of random matching) for CSMs is then

    FDR=(TD-DD)/TT

xiFDR performs this calculation by finding the score cutoff that corresponds to a desired FDR number at the level of analysis (see next section). 

2. Self and heteromeric crosslinks

The chance of random matching a crosslinked peptide pair spectrum within a protein sequence (self crosslink) is different from that of matching a crosslink between 2 different proteins (heteromeric crosslink) as shown in [Lenz et al. 2021](https://www.nature.com/articles/s41467-021-23666-z) Fig. 1. Moreover, spectra of heteromeric crosslinks tend to have lower signal-to-noise than those of self crosslinks. For these reasons, FDR calculation on self and heteromeric crosslinks should be performed separately. xiFDR automatically splits FDR calculations between self and heteromeric crosslinks reporting them separately in the "results" page. 

3. Error propagation

Very often, the final product of a crosslinking MS experiment is not a list of spectral matches, but rather a list of which residues were found crosslinked to which other residues ("crosslinks" or "links" or "crosslinked residue pairs"). To obtain this list, CSMs must be aggregated into which peptides are paired with each other (multiple spectra can come from the same peptide pair), and those need to be aggregated into crosslinked residues (multiple peptide pairs can come from the same crosslinked residues, because of modifications and miscleavages, for example). Finally, if we are interested in producing a protein-protein interaction (PPI) network, error has to be propagated from residue pairs to PPIs. Correct error propagation from lower to higher levels of result aggregation has big implications for the error rate of the final 
result, as
covered in [Lenz et al. 2021](https://www.nature.com/articles/s41467-021-23666-z) and [Yugandhar et al. 2020](https://www.nature.com/articles/s41592-020-0959-9). Thus, in reporting crosslinking MS data on which residues are crosslinked to each other, FDR filtering should be set at the link level and not at the CSM level.

![ErrorPropagation](https://user-images.githubusercontent.com/44289027/203276815-046f76bf-2c1d-49a0-ac7a-479c73b32a3a.PNG)

(target-target matches in grey, matches involving decoys in red. Error control should account for aggregation of matches from lower levels to higher levels.)

xiFDR allows the user to filter for the desired FDR at the level or levels of interpretation of the results. For example, data may be filtered at 5% FDR at the link level, and 10% at the PPI level. Error is propagated by aggregating target and decoy matches from lower levels with a sum of squares approach. The "boosting" feature sets thresholds at lower error levels to reduce the number of decoys at higher error levels, resulting in more matches passing FDR at residue pair or protein pair levels.



#### Terminology
| Term                                     | Description                                                                                                                                     |
|------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------|
| CSM                                      | crosslink spectra match                                                                                                                         |
| PSM                                      | peptide (or crosslink) spectra match                                                                                                            |
| peptide pair                             | 2 peptides crosslinked to each other, no matter the crosslink site. Multiple CSMs can support the same peptide pair                             |
| Link/crosslink/residue pair              | 2 residues crosslinked to each other. Multiple peptide pairs can support the same crosslink                                                     |
| Protein pair/PPI                         | 2 proteins found crosslinked to each other, no matter the number/identity of the residue pairs. Multiple residue pairs can support the same PPI |
| linear peptide                           | a non-crosslinked peptide                                                                                                                       |
| Self link                                | a link within a protein sequence (includes links between multiple copies of the same protein)                                                   |
| Heteromeric link                         | a link between 2 different protein sequences                                                                                                    |
| Ambiguity                                | when peptide pairs or residue pairs can be mapped to multiple proteins because these share the same sequence                                    |
| Protein Group                            | the group of proteins sharing the sequence making up the ambiguous group                                                                        |
| Prefilter                                | a score filter applied prior to FDR calculation to both target and decoy matches                                                                |
| Local FDR/poterior error probability/PEP | the estimation of FDR of a particular match based on a windowed approach or a fitting of the score distributions.                               |
| Boosting                                 | a grid search approach optimising settings to maximise the number of matches passing a given FDR threshold                                      |
| DeltaScore                               | The score of the best explanation of a CSM/residue pair etc. divided by the second best explanation                                             |
| Conservative                             | An explanation of a spectral feature where non-lossy matches are weighted more heavily than lossy ones.                                         |
| Coverage                                 | number of fragments matched / max number of theoretical fragments                                                                               |

## Calculating FDR for crosslinking MS data with xiFDR

xiFDR reads in search results from xiSEARCH in .csv format or from other search engines. It then performs FDR filtering and writes out results in .csv and .mzIdentML formats.

### The interface

The interface should be loaded via the executable files "startWindows.bat"/"startUnix.sh"/"startMacOS.command" depending on whether the program is run in Windows, Linux or Mac. It is not advisable to run the .jar file directly as this may not detect the correct encoding, which will impact writing the .mzIdentML file.

The interface provides several tabs. 

The first tab, "input", is used to read in results from search engines, define column names if necessary, and apply prefilters prior to FDR estimation. 

The second, "FDR settings", performs the actual FDR calculation at the desired error level.

The "Results" tab shows a summary of the FDR calculation. It also allows the user to save FDR-filtered crosslinking MS results in .csv or [HUPO-compliant mzIdentML 1.3.0 format](https://www.psidev.info/mzidentml#mzid13) for deposition in databases.

There is then a tab displaying a log of the FDR calculation and an "about" tab with software information.

### Loading crosslinking MS search results

The top half of the input tab relates to reading in results, while the bottom relates to mapping column names, which is needed only in case of loading results that are not coming from xiSEARCH.

#### Loading search results from xiSEARCH

To load results from xiSEARCH, only use the top half of the "input" tab containing csv-file, xi config and fasta.

xiSEARCH outputs results in a .csv file that contains all crosslinked peptide spectrum matches (CSMs) and scores. It includes target-target matches, target-decoy and decoy-decoy matches.

Select the xiSEARCH .csv results file. The .fasta database and the .config file used in the search should also be uploaded in the respective boxes if available. Fasta and config file are required  to write out results in .mzIdentML format. 
 **These should be set before pressing "read".** See the [xiSEARCH documentation](https://github.com/Rappsilber-Laboratory/XiSearch?tab=readme-ov-file#the-files-tab) for details on what .fasta files generate meaningful accessions, protein names and descriptions in mzIdentML.

Any prefilters on spectral quality (see below) should be set by ticking the "filter" box.

Press "read" to initiate reading of search result. Monitor  the bottom left of the window for the message "finished reading file" and the bottom right for memory usage. 

If xiFDR is very slow or crashing, restart the program by increasing the allocated memory by editing the java -Xmx option in the startup .bat/.sh/.command file. The default is -Xmx3G, providing 3 Gb of RAM. Often, this will not be enough for searches with dozens of runs and thousands of matches. Large searches may require tens of Gb to read in. In this case, prefilters (see below) are a useful way to reduce memory requirements if needed.


#### prefilters
xiSEARCH provides many features of CSMs that may be used to prefilter the results prior to FDR estimation. Doing so equally on targets and decoys prior to FDR estimation retains the accuracy of FDR, while doing score filtering or other filters post FDR estimation generates results with unknown error rates. We recommend doing this only after having a look at the FDR-filtered results without any prefilters. 
If running without pre-filters already generates warnings about not enough TTs at the level of interest (links or protein pair), pre-filters are often not making the FDR more reliable.

The prefilters may be toggled in the "input" tab by clicking the "filter" option on the top left side of the interface. These are generally used after a first look at results without prefilters if spectra of low quality are still passing the FDR. A sufficient number of target-target CSMs are still needed to estimate an accurate FDR. xiFDR will therefore warn you once you attempt to calculate FDR if there are too many decoy or not enough target hits at each particular FDR level ("not enough TT"). Some of the commonly used prefilters (for FDR calculation performed on xiSEARCH results, especially on large scale searches) are

| Filter name                                         | Description                                                                                                                                                    | Commonly set to               |
|-----------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------|
| peptide 1/2 unique matched                          | minimum observed fragments for peptide 1 or 2                                                                                                                  | >3 (4-5 to be more stringent) |
| fragment unique crosslinked matched conservative    | minimum fragments across both peptides containing the crosslink site                                                                                           | >0                            |
| peptide 1/2 unique crosslinked matched conservative | impose to have only spectra where both peptides have a fragment of the other                                                                                   | >0                            |
| delta                                               | only spectra where the best match is x times better than the second best by score                                                                              | >1.2                          |
| peptide1/2 CCPepFragmentDoubletCount                | For cleavable crosslinkers, only consider spectra where at least X doublets are found on peptide 1/2                                                           | >0                            |
| fragment CCPepDoubletCount                          | For cleavable crosslinkers, only consider spectra where at least X doublets are found on either peptide. usually set to greater than 0 in large scale searches | >0                            |


Notice that these are not meant to be used blindly all at once! In fact, at most one or two of these should be used depending on the nature of the dataset.

MS-cleavable crosslinkers present several advantages. Their signature crosslinker stubs and peptide doublets help to provide extra confidence, that the peptide masses are correct. This increases the chance that the peptides themselves are correctly identified. xiFDR can make the most out of these features by prefiltering spectra on a minimum of crosslinker stubs observed, and then boosting on stubs and doublets.


#### Loading search results from other crosslinking MS search engines
xiFDR can directly read search result files in te standard .mzIdentML format in the mzIdentML tab. 

If the results are not in .mzIdentML format, search results should be read in via the csv tab. To do so, a minimal csv has to be generated containing the following information:



Following are lists of columns that have to be or can be provided. The name of the columns do not need to match exactly - but can be manually assigned. I.e. if instead of `scan ` the column in the CSV file is called `spectrum number` then these columns can be associated in the interface. Internally  xiFDR has a list of known alternative names for a some of the columns and will try to automatically match the right columns. This can be expanded by selecting "intelligent column matching" before selecting a file - but in both one should check that the correct columns are used.

The minimal set of columns that need to be present are:

| column name | Description |
|-------------| ----------- |
| run         | raw file name that the spectrum derived from |
| scan        | scan number within that run |
| peptide1 | sequence of the first peptide |
| peptide2 | sequence of the second peptide |
| peptide link 1 | which residue of the first peptide does the cross-linker attach to |
| peptide link 2 | which residue of the second peptide does the cross-linker attach to |
| is decoy 1 | is the first peptide from the decoy database |
| is decoy 2 | is the second peptide from the decoy database |
| precursor charge | charge state of the precursor ion |
| accession1 | protein accession number(s) for the source of the first peptide |
| accession2 | protein accession number(s) for the source of the second peptide |
| peptide position 1 | position(s) of the first peptide in the protein(s) |
| peptide position 2 | position(s) of the second peptide in the protein(s) |
| score | score of the spectrum match |

The decoys have to be reported in the search results and the accession of decoy proteins need to pairable with target proteins. Meaning a decoy protein protein needs to have the same accesion number as a target protein plus an optional prefix. xiFDR recognises following prefixes `REV_`, `RAN_` and `DECOY:` to the accession-numbers.

For ambiguous peptides xiFDR assumes that all proteins for the given peptide are either target- or decoy-proteins.

There are several more columns mappable.

Some column will get their value assigned/guessed from the top columns:

| column name | Description |
|-------------| ----------- |
| psmid       | a ue ID for the PSM - if not given will be defined based on run and scan |
| peptide length 1 | length (in amino acids) of the first peptide - if not given, then it will be guessed from the sequence|
| peptide length 2 | length (in amino acids) of the second peptide - if not given, then it will be guessed from the sequence |

When psmid is given run and scan number become optional.


Additional columns may be provided. These can be used by xiFDR's boost or filter function to improve the results of the FDR calculation:


| column name | Description |
|-------------| ----------- |
| delta score | the delta score of the match |
| peptide coverage1 | how well is peptide explained |
| peptide coverage2 | how well is peptide 2 explained |
| minimum peptide coverage | how well is the worse explained peptide explained |
| peptide1 fragments | How many ue fragments where assigned to peptide 1 |
| peptide2 fragments | How many ue fragments where assigned to peptide 2 |
| peptides with stubs | For how many peptides was an cross-linker stub identified |
| peptides with doublets | For how many peptides was an doublet of cross-linker stubs identified |
| negative grouping | if some matches have an inherently higher chance to be false positive then they can be flagged here |
| positive grouping | if some matches have an inherently lower chance to be false positive then they can be flagged here |

The last two columns are used to group matches for which the FDR is calculated separately. This requires that enough data are present and should be considered experimental.

Additionally there are some columns that are mainly just passed through to provide a "nice" output:

| column name | Description |
|-------------| ----------- |
| crosslinker | name of the cross-linker involved in this PSM |
| crossLinkerModMass | mass difference between the sum of the non-cross-linked peptides and the cross-linked peptides |
| experimental mz | experimental precursor M/Z |
| calculated mass | calculated mass of the precursor |
| Protein name1 | a name for the first protein |
| Protein name2 | a name for the second protein |
| description1 | description of the first protein |
| description2 | description of the second protein |
| info | arbitrary info field to be forwarded to the results table |

For xiSEARCH results there some more columns that are important for the mzIdentML export that get automatically mapped. 


#### other settings in input tab

| setting      | Description                                                                                                             | When to use                                                                                                                                                                                                                   |
|--------------|-------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Flag modifications | splits FDR calculation for modified and unmodified peptides                                                             | EXPERIMENTAL: try if spectral quality of modified and unmodified peptides differs and there are sufficient numbers of matches in both modified and unmodified to perform FDR calculation                                                        |
| Normalize    | Normalize score distributions according to FDR or the decoy distribution. Suggested settings "FDR based" or "AutoScore" | generate results to be merged with other searches with different settings. For example, normalize scores to later combine FDR results from searches with HCD and ETD fragmentation, that yield different score distributions. |
| Forward Columns    | Text box for a regular expression to take columns from the input .csv file and include it in the xiFDR output files | include for example retention time in the result file to perform additional rescoring based on machine learning, for example with [xiRT](https://www.nature.com/articles/s41467-021-23441-0)|



### Performing the FDR calculation

In the "FDR settings", one can perform the actual FDR filtering. 

By default, the view is set to "reduced FDR", which shows just the basic settings. The cutoff is set at 5% at the residue pair level, meaning the  error will be propagated so that 5% of the residue pairs correspond to a wrong/random match. The "boosting" features is enabled (see below for more details). These are perfectly acceptable FDR filtering settings for experiments aimed at characterising the crosslinks in a purified protein complex and should give a good idea of the number of crosslinks detectable with reasonable certainty in the sample. In analyses of cellular fractions or searches with hundreds of proteins in the database, it is advisable to also include an FDR cutoff at the "Protein Pairs" level. Similarly, in analyses devoted to method development on the quality of spectra, a filter at the "PSM" is advised.
Whenever multiple cutoffs are selected (e.g. both residue pairs and PPI level at 5%), xiSEARCH will report only the entries at lower levels that also pass the upper level. 
For example, only the residue pairs that will pass both the residue pair and PPI level thresholds will be considered as passing the FDR. Notice that in these cases the results files at lower levels will have less decoy lines than the stated FDR threshold at that particular level (as only the targets and decoys passing also the higher level are reported). Thus, these files cannot be used to recompute an estimate of the false positive rate, as the decoy-based FDR would drastically underestimate the FPR. Thus: with both residue pair and PPI level FDR, the number of decoys in the "Links" csv cannot be used as a basis of analyses that recompute apparent decoy-based FDR. For these analyses, the use of a threshold at a single level (the one of interest) is recommended.


If these settings are satisfactory, press "calculate". 

Otherwise, more advanced options are available if one ticks the "complete FDR" setting. Further ticking the box "more" will present the full set of options for FDR filtering in xiFDR. Hovering the mouse over checkboxes or up/down arrows provides a brief explanation of the setting.

All FDR settings:

| Setting                                                     | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           | when to use                                                                                                           |
|-------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------|
| max FDR                                                     | the FDR setting at that particular level (PSM, peptide pair etc.)                                                                                                                                                                                                                                                                                                                                                                                                                                                     |                                                                                                                       |
| Local FDR                                                   | whether to compute a per-PSM/per-peptide pair etc. score based on a windowed FDR approach. Black box=compute local FDR, ticked box=calculate and filter by local FDR, white box=do not calculate local FDR                                                                                                                                                                                                                                                                                                            | to analyse distributions of scores and FDR or to generate a per-crosslink error estimate e.g. for structure modeling. |
| min Pep. length                                             | minimum length of peptide 1 and peptide 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                             | increase to 6 if noise matching of short sequences                                                                    |
| min supporting peptide pairs Protein Group/Residue Pair/PPI | minimum number of peptide pairs necessary to aggregate to a given level                                                                                                                                                                                                                                                                                                                                                                                                                                               | define to have more robust network especially at the PPI level.                                                       |
| min TD chance                                               | require enough target matches so that these number of decoy matches would still pass the target FDR. In other words, a minimum number for reliable FDR calculation                                                                                                                                                                                                                                                                                                                                                    | Use for more accurate FDR calculation. By default, violations of this only cause warning messages, because "ignore validity checks" is ticked  |
| Ignore validity checks                                      | Allow FDR calculation even if number of decoys is not sufficient for sufficiently accurate estimation of FDR                                                                                                                                                                                                                                                                                                                                                                                                          | Untick and use in combination with min TD chance to only get results with accurate FDR calculation                    |
| Ambiguity                                                   | increase number to only allow matches with at most N ambiguity. E.g. if "ambiguity" is set to 1 in the residue pairs, residue pairs that could come from more than 1 place in the sequence or more than 2 ue proteins are discarded                                                                                                                                                                                                                                                                               | Set to a number if not interested in considering ambiguous crosslinks at all |
| Unique PSMs                                                 | Aggregate spectra matches so that only the top scoring spectrum match per peptide pair, charge state, and linkage sites combination is considered                                                                                                                                                                                                                                                                                                                                                                     | on by default. If unticked the CSM/PSM FDR can become unreliable                                                                       |
| Group by PSMs                                               | split FDR calculations into groups of protein pairs having different numbers of matches. In theory splittng FDR for highly abundant proteins and low-abundant                                                                                                                                                                                                                                                                                                                                                         | off by default. Experimental                                                                                                                 |
| No consecutive                                              | remove matches from consecutive peptides. When the two peptides of a crosslinked peptide match are consecutive in sequence, then the same precurosor mass could also come from a single missedcleaved peptide, spanning the joined sequence,  with a crosslinker modification (hydrolised crosslinker). These two cases not only share the same precusor mass but also a lot of possible MS2 fragments, making these sometimes hard to distinguish.  These crosslinks are typically structurally not very informative | Consider switching on in searches were protein where individual protein structure is the focus and short sequence distance crosslinks have little value |

Settings in the "more" panel:

These are minimum filters that essentially act as prefilters prior to FDR calculation. The difference between setting them here versus in the "input" tab is that these parameters can then be part of the boosting grid search (see below).

| Setting                | Description                                                                | when to use                                                                                         |
|------------------------|----------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------|
| min. peptide fragments | minimum observed fragments per peptide                                     | 0 by default, increase number in SDA searches                                                       |
| min. coverage          | number between 0 and 1. Minimum fraction of theoretical fragments required | included in boosting                                                                                |
| min. coverage          | number between 0 and 1. Minimum fraction of theoretical fragments required | included in boosting                                                                                |
| min. peptide stubs     | minimum number of fragment stubs observed                                  | important for MS-cleavable crosslinkers, included in boosting                                       |
| min. peptide doublets  | minimum number of  peptide doublets observed                               | important for MS-cleavable crosslinker searches, included in boosting                               |
| min. score             | minimum spectral match score                                               | for searches with a small database where FDR may not be properly computed, not included in boosting |
| boost separately|boosting in two steps - first on the lower level FDRs and in a second step change optimize these pre-filter| on by default                                                     |

The settings in "define groups" are currently very beta and unsupported. They allow for splitting the dataset further down into custom groups for FDR calculation.

#### Boosting
xiFDR's boosting feature performs a grid search to optimise FDR settings at lower levels to reach the maximum number of matches passing validation at the desired FDR level. For example, enabling boosting for a residue pair-level FDR of 5% will tweak PSM, peptide pair and protein group level FDR cutoffs to maximise the number of residue pairs passing the 5% FDR threshold.

Boosting is performed with a grid search of parameters as described in [Fisher et al. 2017](https://pubs.acs.org/doi/pdf/10.1021/acs.analchem.6b03745). 

The user may control which parameters are part of boosting by changing the selection in the "boosting includes" button.

The "steps" controls how many steps of the grid search per parameter are tested each round of optimization. The "between" box ensures that boosting is performed to maximise the number of heteromeric residue pairs/PPIs etc. passing FDR rather than the overall number. This is recommended for searches where the goal is to produce a protein-protein interaction network and where large numbers of heteromeric crosslinks are available.

We recommend leaving boosting on and selecting "between" if desired. For experiments with MS-cleavable crosslinkers, we suggest also boosting on minimum peptide doublets by toggling those on in the "boost includes" menu.

FDR calculations with boosting enabled can take some minutes to conclude.

#### EC-filter
Ticking the EC-filter box in the complete FDR view enables the EC-filter, which filters protein herteromeric residue pairs (between) based on whether they occur between proteins that are also identified with crosslinker modified linear peptides or with self links. See [Fischer&Rappsilber 2024](https://www.embopress.org/doi/full/10.1038/s44320-024-00057-2) for more details.
The EC-filter is mostly useful for searches involving a large number of proteins in the database (at least hundreds) of which a good portion is not expected to yield residue pairs.
This typically means searches involving in situ/in cell/in lysate crosslinking.

### FDR suggestions & common mistakes

The extreme flexibility of xiFDR requires some careful use. Here are some of our suggestions and ways to spot that the FDR calculation may not be giving sensible results. When in doubt, keep it simple! 

#### Tips & tricks
- For searches with cleavable crosslinkers such as DSSO, set prefilters (see below) on the number of peptide doublets and boost on that parameter
- For searches with unspecific crosslinkers like SDA, include prefilters (or boost) on minimum number of fragments. Include also a prefilter on "fragment unique crosslinked matched conservative" and exclude consecutives. This ensures that only spectral matches with at least one fragment proving a crosslinked peptide pair are considered.
- Sometimes, filtering for identifications supported by fragments with no neutral losses greatly cleans up the dataset (prefilter on "peptide1 uniqued matched non lossy"/"peptide2 uniqued matched non lossy")

#### Watch out for:
- The program warns if there aren't sufficient targets to estimate FDR accurately. If this warning refers to the level of analysis you are interested in, there is likely almost nothing in the data. Unless prefilters were set way too stringent.
- The FDR calculation rests on the assumption that target-decoy pairs have twice the chance of random matching than decoy-decoy pairs (TD+DT vs DD). If in your results you have more decoy-decoy (DD) than target-decoy (TD), your FDR evaluates to a negative number and becomes meaningless. This may be a sign that there are very few crosslinks in the datasets, or that prefilters are not working as intended.


### Results summary
After the calculation is complete (check the log tab, the lower bottom left of the window or the "calculate" button becoming clickable again), the results tab presents a summary of all matches passing the FDR threshold with the number of decoys and target-decoys in each category. 
Remember the ratio of targets and decoys passing determines the accuracy of the FDR estimation. 
This is where one must **carefully** check if the number of target-target (TT) identifications as a fraction of total identifications makes sense with the stated FDR.

As an example, seeing the following in the residue pair between box:

    123 (115 TT)

makes sense if one requests an FDR of 5% at the residue pair level. On the other hand, seeing

    50 (20 TT)

makes no sense and should only happen as a result of too many decoy-decoy passing the filter (which triggers a warning) and breaking down the basic FDR assumptions. 
Because of the FDR formula, too many decoy-decoy may actually artificially reduce the apparent FDR (TD-DD is in the numerator). 
Thus, if the number of target-target hits does not square with the FDR (and the program should have printed a warning about it), the results should be considered invalid for the level of analysis requested.

### Writing out search results

Results may be written out in xiFDR .csv format (csv tab) or in mzIdentML1.3.0 format. 
The csv format writes human-readable tables of identifications and the mzIdentML writes a standards-compliant results for later deposition or upload to [xiView.org](https://www.xiview.org).

The mzIdentML can only be generated if xiFDR is provided with .fasta file and search config file (from xiSEARCH) prior to reading in results.

#### csv output 
will generate several files:

- Summary file: a file containing a summary of the FDR calculation, results and all settings
- CSM file: a file containing all crosslink spectra matches passing all thresholds. This file can be uploaded to xiview.org for visualization
- peptide pairs file: a file containing all peptide pairs passing all thresholds
- protein groups file: all protein groups passing all thresholds (i.e. including ambiguity when a peptide cannot be assigned to an individual protein)
- Links file: all the crosslinked residue pairs passing all thresholds
- PPI file: all the protein-protein interactions passing all thresholds.

Each file contains information about the FDR, posterior error probability (local FDR), target/decoy nature of each match, as well as the peakfile of origin.

Beware that if multiple FDR thresholds are set, *only* the matches passing all levels of aggregation are included in the results. In other words, the csv files for a search with both a 5% residue pair FDR and 5% protein pair FDR will only include residue pairs with a 5% or better FDR that also lead to a protein pair with a 5% or better FDR.

This can result in a lower number of decoys at a given FDR level than the reported FDR number. Meaning, when filtering for both 5% residue pair FDR and 5% protein pair FDR the resulting residue pair file will contain only decoys that reach the protein pair level. This should not be misconstrued as having only as many false positives left at the residue pair level as you have decoys residue pairs: the protein pair filter basically invalidates any follow up decoy based error estimation at the lower residue pair level. Thus, any FDR at a higher level of aggregation (eg. protein pair vs CSM) will prevent a meaningful analysis of false positives via decoys at lower levels.

#### mzIdentML output
will generate a single file .mzIdentML compliant with standards. The file can be deposited in ProteomeXchange repositories (PRIDE) or uploaded to xiview.org for visualization. It contains information about the search results, peaks and validation. 
Reading in the xiSEARCH "config" file in the "input" tab is required for mzIdentML output. In a PRIDE deposition, selecting it as a crosslinking MS dataset will deposit the results in [PRIDE crosslinking](https://www.ebi.ac.uk/pride/markdownpage/crosslinking)
The mzIdentML output will also generate a summary csv file with all FDR settings.

On some new versions of java (such as 17.0.6 2023-01-17), xiFDR produces the "ExceptionInitializerError" when writing mzIdentML output. In this case, edit the java command in the launcher script (on windows, startWindows.bat) to:

    java -Dfile.encoding=UTF-8 -Xmx3g --add-opens java.base/java.lang=ALL-UNNAMED -jar %mypath%xiFDR-2.1.5.2.jar

Which should fix the error (remember to edit the Xmx option for memory allocation as appropriate).


### Running xiFDR from the command line

xiFDR may be run in the command line for automated processing or cluster jobs.

The full range of options and their description is available with

    java -jar xiFDR.jar --help

Allocate memory with the java -Xmx flag. E.g. 

    java -Xmx5g -jar xiFDR.jar

To allocate 5 Gb of RAM. 

Prefilters are not supported in the command line version of xiFDR. However, prefiltering the .csv input file as desired may be done in python/pandas, R, or any other tool prior to loading the input file into xiFDR.
