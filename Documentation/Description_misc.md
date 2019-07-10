Misc internal comments for the project
======================== 

# Technical details and early processing

## Solutions and concentrations

See file "ca protocol calculations.xls" for a comparison of concentrations that different groups used in their experiments.

### 131213 experiment: Latency calibration

Old-school (ephys paper –compatible) stimuli, recorded by the camera, and compared. Stored in the main folder, not in respective Data folder.
 
Based on this data (average lightness of the fiber, negative thresholdoing at a 99.9% value), average latencies compared to a flash (in ms) are: [0 225 172 468] for flash, crash, grid and back respectively. The reason for latency adjustment for "grid" and "crash" is that for these stimuli a one-pixel-size object (in the very beginning of the stimulus) would not actually be shown on the screen, even though the program would try to produce it (due to computer->VGA->NTSC->LCD conversion). Therefore actual stimulus is shown a tiny bit later than Matlab starts generating it. For "back" stimulus the situation is even more dire, as the stimulus starts off-limits of the fiber, and only some time later enters the fiber image.

Note that this approach differs from what was used in "Tectal processing" paper, as there we assumed "flash", "grid", "ramp" and "crash" having equal latencies. Check out the file "crampf - Eye recordings.xls", tab "Stimuli" for comparison.

I think it would make sense to introduce discontinuity with the previous paper (ephys), and re-adjust stimuli, so that they start and end simultaneously (otherwise dynamics comparisons will be really problematic). Some pilot experiments however (131209 and 131213) still use old stimuli, and thus require a heavy latency adjustments.

### Data export and transformation

Beware: by default NES exports to Excel not full data, but only a brief report on data, which looks like a nice square table, but is limited to 74 first ROIs. To get all data, switch export model to "full data", in which case you won't get a square table, but rather information for each ROI given one under another. An unfortunate surprise you should be prepared for is that in "Report" mode time axis is given in seconds, while in "Full data" mode it is given in ms. I spent about 3 hours troubleshooting my program, trying to understand why it works on old format data, but does not work on new, before I realized that the time axis is all wrong.
