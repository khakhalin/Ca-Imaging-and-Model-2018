Misc comments on the project
======================== 

File: ca protocol calculations.xls for a comparison of concentrations that different people from different labs use in their experiments.

## Early processing, comments

### Intensity measurement

### Data export and transformation

Beware: by default NES exports to Excel not full data, but only a brief report on data, which looks like a nice square table, but is limited to 74 first ROIs. To get all data, switch export model to "full data", in which case you won't get a square table, but rather information for each ROI given one under another. An unfortunate surprise you should be prepared for is that in "Report" mode time axis is given in seconds, while in "Full data" mode it is given in ms. I spent about 3 hours troubleshooting my program, trying to understand why it works on old format data, but does not work on new, before I realized that the time axis is all wrong.