Comments on how data is stored
============

A structure of the original dataset:

1.	Excel file with fluorescence signals
2.	TIFF containing ROIs – Not necessary, as ROI (x,y) coordinates are also stored in same Excel data file.
3.	A print-screen of ROIs (for further processing in imageJ)
4.	.mat file with a Matlab structure (S) containing all data from this Excel file, and also spikes reconstructed

### 131213 Latency calibration

Old-school (ephys paper –compatible) stimuli, recorded by the camera, and compared. Stored in the main folder, not in respective Data folder.
 
Based on this data (average lightness of the fiber, negative thresholdoing at a 99.9% value), average latencies compared to a flash (in ms) are: [0 225 172 468] for flash, crash, grid and back respectively. The reason for latency adjustment for "grid" and "crash" is that for these stimuli a one-pixel-size object (in the very beginning of the stimulus) would not actually be shown on the screen, even though the program would try to produce it (due to computer->VGA->NTSC->LCD conversion). Therefore actual stimulus is shown a tiny bit later than Matlab starts generating it. For "back" stimulus the situation is even more dire, as the stimulus starts off-limits of the fiber, and only some time later enters the fiber image.

Note that this approach differs from what was used in "Tectal processing" paper, as there we assumed "flash", "grid", "ramp" and "crash" having equal latencies. Check out the file "crampf - Eye recordings.xls", tab "Stimuli" for comparison.

I think it would make sense to introduce discontinuity with the previous paper (ephys), and re-adjust stimuli, so that they start and end simultaneously (otherwise dynamics comparisons will be really problematic). Some pilot experiments however (131209 and 131213) still use old stimuli, and thus require a heavy latency adjustments.