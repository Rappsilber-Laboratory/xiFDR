
2.4.1

* Enforce tsv/txt/csv extension for write out of results
* GUI now defaults to self+between boosting
* BugFix: csv.getDouble() was initialised with int value

2.4

* Self and heteromeric matches can now boosted sequentially and (mostly) independently
* Protein homomultimeric self links are now treated separately
* BUGFIX config modifications missing leading to unknoqwn modification warning during mzIdentML export

2.3.10

* new csv argument --filter
* GUI filter can be switched between simple  text field and GUI input
* mzIdentML protein name will be set to accession if no name is given

2.3.9

* Write out mzIdentML crosslinker stub information to donor and acceptor
* no need to restart anymore - should produce valid UTF-8 mzIdentML without it
* db reading fix for large scan_ids

2.3.8

* Write out crosslinker stub information to mzIdentML
* updated to XLMOD java library 1.1 to reduce "unknown modifications"
* BugFix: mzIdentML 1.3 empty protein descriptions no longer written out
* BugFix: mzIdentML 1.3 export still had some mzIdentML 1.2 namespace attributes
   * currently using a modified version of jmzIdentML (pull request pending)
   * allow matching to crosslinker without mass definition (e.g. DSSO)
* BugFix: mzIdentML export for non-covalent peptide pairs

2.3.7

* BugFix: xiview csv file can miss data

2.3.6

* mzIdentML 1.3 export
* mzIdentML import/export possibility to only consider matches flagged as passing
* dedicated xview csv file as alternative to mzIdentML
* --forward command line argument also for DB interface
* CSV input without individual petide scores:changed how CSM score gets distributed to proteins
* reset button for FDR settings in the complete interface
* BugFix: if one of delimiter and quote where given both where tried to be guessed again
* Bugfix: summary-file protein group FDR setting was wrongly reported
* Bugfix: odd unicode error with NaN and infinite
* Bugfix: mzIdentML export with newer java-versions
* Bugfix: iFDR at last element turns negative

2.3.5

* dedicated xview csv file as alternative to mzIdentML
* --forward command line argument also for DB interface
* Bugfix: commandline mzIdentML writing without owner information
* Bugfix: reading multiple files with different forwarded columns
* Bugfix: database reading

2.3.4

* Validity check command-line arguments
* Bugfix: crash if no xi version was defined/detected
* Bugfix: xi Version input field to small

2.3.3

* Write out the highest CSM score for each raw-file also into the PPI result file
* BugFix: if several csv-files where read in only the highest score for raw files
   found in one of the source CSV files got reported in the peptide pair and residue pair result files

2.3.2

* Bugfix: shared topranking matches, validation status was only checked for first one

2.3.1

* Support for writing xiSEARCH version to mzIdentML
  * versions can be defined by --xiversion= argument or in the GUI - but if not given it will be guessed from the xiSEARCH result file name
  * limitation: currently only the xiVersion associated with the first input file is forwarded to mzIdentML
* CSV input permits different xiSEARCH configs to be associated with input different files (also mainly relevant for mzIdentML export)  

2.3

* xiSEARCH2 mzIdentML support

2.2.5

* select correct xiSEARCH score column by default

2.2.4

* support for matching ensemble accession numbers
* Descriptions1 column in ppi file renamed to Description

2.2.3

* Speed up
* CSV parser related changes
* BugFix: warning about different number of protein accession positions and descriptions are only show at most 5 times

2.2.2

* BugFix: command line version used the filter to unique PSMs additionally also as the boost between flag.
* BugFix: Defaults for minimum peptide length was different for command line and GUI
* BugFix: Defaults for length group was different for command line and GUI 

2.2.1

* CSM-level boosting from command line
* Decoy protein names in mzIdentML are changed to "decoy"
* automatically restart xiFDR with needed arguments to run on newer JDK versions
* some changes that hopefully addresses crosslinker related issues with mzIdentML export

2.2

* explicitly definable prefix used in front of decoy accessions
* mzIdentML owner info arguments for command-line version
* Support for expected crosslink filter (ec-filter)
* changes to protein FDR (validity checks are more stringent here by default and can't be disabled)
* Boosting between should not result in zero self links do to validity checks
* commandline arguments for mzIdentML document owner information
* xiSEARCH2 DB interface
* prevent database disconnect during running query
* speed up database reading (actually database result processing)
* BugFix when starting a query for a large search from the database the connection could get closed before the query was finished
* BugFix mzIdentML export protein n-terminal modifications where not flagged as such
* BugFix re-parsing peaklist info from peaklist-files failed (needed potentially for mzIdentML export)
* BugFix optional boosting on cleavable crosslinker stubs failed
* BugFix reading peptide fragment subscores from xiSEARCH csv files

2.1.5.5

BugFix for reading from instable db connections

2.1.5.4

* Don't show a popup every time we have a DB-disconnect during reading in

2.1.5.3

* Reading from DB with unstable connection improved

2.1.5.2

* CSVParser updated

2.1.5.1

* When reading from CSV-automatically shift scores if negative scores are observed
* BugFix: "Ignore groups" setting was ignored during boosting

2.1.5

* BugFix mzIdentML export
* new commandline option for xisearch input "--lastowner" to not ask for mzIdentML owner information

2.1.4

* BugFix CSV: for empty description fields with ambiguous matches
* CSV: More tolerant against unquoted ";" in description
* Filter for peptide stubs and doublets are exposed in the GUI
* Changed how currently hidden filter settings are treated when starting a new FDR calculation
* BugFix Boosting: Hopefully fixed issue of boosting settings sometimes incompletely transferred to final result
* BugFix DB: Filter for cleavable peptide filter had an error
* BugFix DB: crosslinker-mass not always assigned correctly to match
* BugFix DB: Reading search from DB without stored sub-scores fails
* BugFix DB: Overwriting Validation states produced wrong results
* BugFix in 2.1.3 was incomplete.
* DB: LinkSiteScore is read out from DB
* DB: Overwriting validation states in DB now first deletes all validation values of the desired state

2.1.3

* BugFix: for writing to many matches out for validation.

2.1.2

* BugFix for error introduced in 1.4.3 that can result in drastically fewer identifications

2.1.1

* easy arbitrary subscore filter when reading from DB

2.1

* Non-covalent matches are written out into separate files
* If no linear matches where found, then the output files for these ar not generated
* BugFix: Boosting of PSMs failed

2.0

* Version to go along with PPI FDR Paper
* command-line version exposes the filter against consecutive peptides
* protein accessions and descriptions can now contain ";" but need to be quoted correctly
* extended the number of recognised column-names
* some preliminary support for cleavable crosslinker specific boosting
* if retention time is found in the input it will be forwarded to the output(column-names checked: retention time, elution time, elution time start, rt [case-insensitive and spaces can be left out or replaced by '_'])
* Between and self FDR shown for residue pairs and PPIs

1.4.6

* BugFix: command-line version did always boost

1.4.5

* BugFix: DB reading of linear peptides did not update the protein information
* BugFix: command-line bugfix

1.4.4

* renamed internal to self
* Type of Normalisation can be selected
* BugFix: some mzIdentML related fixes
* BugFix: for ambiguous matches read in from CSV
* BugFix: Normalisation had some trouble

1.4.3.1

* BugFix: command-line version was broken

1.4.3

* Changed Database read-in to enable additional sets of filters
* BugFix: If the desired FDR in a subgroup was not reached elements up to the last TD got ignored
* BugFix: Calculate button fails to function
* BugFix: Linear Ambiguous matches lead to duplication

1.4.2.1

* BugFix: Calculate Ranges not having all settings forwarded
* BugFix: if xi-config is provided also "odd" modifications are handled correctly

1.4.2

* BugFix: Boosting of MinPeptideFragmentsFilter not correctly reported back
* BugFix: Dropdown box in column assignment got disabled
* expanded the FASTA matching for reading CSV files
* PEP estimation slightly changed
* Arbitrary Columns can be forwarded from CSV-input to output
* Columns starting with "peak_" are forwarded as double values
* FDR calculation can be reported for a range of FDRs
* Boosting happens in two steps if additional columns for boosting are used

1.4.1

* Some cleavable cross-linker related columns are read in
  * specific grouping for these can be selected
* BugFix for spaces in ambiguous accession list
* Validity check disabled by default  but still reports warning
* Validity check warnings more verbose
* When writing csv-files from the command line mzIdentML files get written out as well
  * a window asking for document-owner will be displayed
* decoy:psms better marked
* Improved synchronisation between settings-panel
* Some adaptation to forward settings from the command-line to the GUI

1.4.0

* Medium complexity level of FDR interface replaces (previous simple one)
* Some cleavable cross-linker related columns are read in
  * specific grouping for these can be selected
* Validity check disabled by default  but still reports warning
* Validity check warnings more verbose
* When writing csv-files from the command line mzIdentML files get writen out as well
  * a window asking for document-owner will be displayed
* decoy:psms better marked
* Improved synchronisation between settings-panel
* Some adaptation to forward settings from the command-line to the GUI
* BugFix for spaces in ambiguous accession list

1.3.37

* BUGFIX local FDR not working correctly
* Forward arguments from command-line to GUI
* min number of fragments readable from csv
* validity checks can be switched off
* gives feedback on what sub-groups got dropped

1.3.36

* changing the logging detail should no longer hang hang up the gui
* Noncovalents are no longer forwarded to Links and PPIs.
* Noncovalents are by default filtered out during reading from the database
* Added optional Filter and boosting by min number of fragments per peptide
* BugFix for zero decoy lists not being reported

1.3.35

* BugFix for ambiguous peptides

1.3.34

* support for local FDR (PEP)
* will read in and try to boost on delta score and peptide coverage
* searches can be normalized by a version of an interpolated PSM-FDR
  * mainly useful for combining searches that have different score-ranges
* using fastutil to reduce some memory requirements
* support for directional crosslinker has been (partially) removed
* bugfix for automatic recognition of column-names from command-line version 
* new output column iFDR for interpolated global FDR

1.2.33

* bugfix - some more settings where ignored for the commandline-version

1.2.32

* bugfix - report factor is now ignored when running from command-line

1.2.31

* BugFix commandline argument outputlocal not recognised

1.2.30

* Option to filter out consecutive peptides from PSMs before FDR
* Renamed BoostIgnores into BoostIncludes as that is actually what a tick means
* Boosting can now be stopped from either setting panel
* Changed parsing of supplied mgf-files a bit now it also supports apl-files
* made 2 the default min TD count
* Grouping between by internal for peptide pairs links and PPIs

1.2.29

* when exporting to mzIdentML - some potentially missing information can be recovered from the originally searched MGF-files
* xiFDR can now work with non-english e.g. german number formats if told that a input/output file is in the according format
* peptide-files contain the peptide positions

1.1.26

* two new columns for CSV-import
  * filescanindex
  * matchrank

1.1.25

* BugFix csv-import

1.1.24

* group names defined as string 
* input can explicitly set positive and negative groups for matches
* minimum potential TD count per group to be reported

1.0.23

* peptide and protein ids are nopw long values

1.0.22

* linear peptides where not checked for preceding or suucceding amino acids

1.0.20

* read multiple csv together
* boosting ignores decoy counts
* summary file contains more infos



