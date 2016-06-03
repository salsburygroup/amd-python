# Visual Statistics
The two python scripts here generate an image of a biopolymer with uncertainty in position shown as shadow.

The script vsigma.py takes a structure file, a trajectory and clustering data where a row contains the frames in a cluster, with the first entry in the row being the median structure in the cluster. It ouputs pictures with the median frame as solid and all frames within 1 standard deviation by RMSD as shadow.

The script vdistribution.py takes two input PDBs -- the median structure as one PDB and the part of the distribution to be displayed as shadow as the second PDB.

You will need to have vmd either in your shell search path or defined as an alias. That is, if you can enter the command "vmd" into your terminal and have vmd launch, you're good to go. The python scripts will spawn a new instance of your login shell in order to run the needed VMD commands. You may see VMD open and close several times. This is normal.

A lot of output will be dumped to your terminal, so you may want to consider sending output to a file using something like "command stuff >&log.txt"

The current process is memory (RAM) efficient but not time efficient. Expect at least 2 seconds per frame in your trajectory. It will also generate temporary files that may take up hard drive space of up to 10mb per trajectory frame in the output directory. This is a worst-case scenatrio, but you may want to set the output directory to somewhere on a high capacity drive. These files will be removed before the script exits, but may not be deleted in the event of an error or crash.

Edge case: If you move these scripts into the same directory as your psf(s), pdb(s) and dcd(s), the way I search for the helper files will fail. 

Example vsigma.py call: 
python /Users/melvrl13/Documents/AMD/AMD-PYTHON/Rendering/VisualStats/MacOS/vsigma.py -s /Volumes/RyanMdata/FUMP10/Folding/weightedSims3200/f10.psf  -t /Volumes/RyanMdata/FUMP10/Folding/weightedSims3200/weighted3200.dcd -c /Volumes/RyanMdata/FUMP10/Folding/weightedSims3200/clustering/cluster5.dat -l 3 -r NewCartoon

(Can take an optional -o for output directory. Defaults to cwd. See -h for details.)

Example vdistribution.py call
python /Users/melvrl13/Documents/AMD/AMD-PYTHON/Rendering/VisualStats/MacOS/vdistribution.py -s sigma.pdb -m mu.pdb -r NewCartoon -o .

On Mac and Windows, you may run into memory issues if your trajectory is larger than 3gb, as the current VMD distributions are 32bit on Windows and Mac. However, I have tried to be as memory efficient as I know how to be. 64bit Linux distributions with 64bit version of VMD should have no problem whatsoever.