# nnh processor

This folder is dedicated to the Marlin processor that will transform the mini-DST lcio files into TTrees in ROOT files.

Usage of the processor itself is very simple, just modify the [``./script/NNH_steer.xml``](./script/NNH_steer.xml) file by replacing ``input.slcio`` and ``output.root`` into, respectively, the mini-DST input file name and the desired ROOT output file name, and then do
```
$ Marlin NNH_steer.xml
```

The [``launchNNHProcessor.py``](./script/launchNNHProcessor.py) is just there to automatize the processing of multiple files. 


> **IMPORTANT: Similarly to the [mini-DST maker](../miniDSTMaker), if you need to access mini-DST files on lyogrid06, do not execute the usual init_ilcsoft.sh script, but instead do:**
>
>     $ . ./script/setEnv.sh
> and do not forget to create a proxy if not done already
>
>     $ voms-proxy-init --voms calice


In order for the script to run properly, the input files must be ordered like this:

```
inputFiles
└───402001
|   | 402001_file_0_mini-DST.slcio
|   | 402001_file_1_mini-DST.slcio
|   | ...
└───402002
|   | 402002_file_0_mini-DST.slcio
|   | 402002_file_1_mini-DST.slcio
|   | ...
└───402003
|   | 402003_file_0_mini-DST.slcio
|   | 402003_file_1_mini-DST.slcio
|   | ...
...
```

### Run the script
```
$ python3 ./script/launchNNHProcessor.py -h
> usage: launchNNHProcessor.py [-h] [-n NCORES] [-p PROCESSES [PROCESSES ...]] -i INPUTDIRECTORY
> [-r] -o OUTPUTDIRECTORY
> 
> optional arguments:
>   -h, --help            show this help message and exit
>   -n NCORES, --ncores NCORES
>                         Number of threads
>   -p PROCESSES [PROCESSES ...], --processes PROCESSES [PROCESSES ...]
>                         ProcessIDs to analyse
>   -i INPUTDIRECTORY, --inputDirectory INPUTDIRECTORY
>                         Path of input files
>   -r, --remote          indicate that files need to be downloaded
>   -o OUTPUTDIRECTORY, --outputDirectory OUTPUTDIRECTORY
>                         output directory
```
For example, it you want to run all the files that are present on lyogrid06 (from a lyoui):
```
$ python3 script/launchNNHProcessor.py -n 10 -f root://lyogrid06.in2p3.fr/dpm/in2p3.fr/home/calice/garillot/ILD/AHCAL -r -o /path/to/wherever/youwant
```
Do not forget to put the ``-r`` argument to tell the script that the files need to be downloaded (from lyogrid06).
It is unnecessary to ask for more than 10 threads on lyoui because the limiting factor is the transfer time between lyogrid06. The output path cannot be remote like the input path.

If you want to run only files from processesID 402007 and 402008, and the input files are present locally : 
```
$ python3 script/launchNNHProcessor.py -n 10 -p 402007 402008 -f /path/to/inputFiles -o /path/to/wherever/youwant
```

The ``launchNNHProcessor.py`` script will create one ROOT file per processID (all the results from each mini-DST file are merged). It will also create a ``logs/`` folder in the output directory to check if something wrong happened (the script will output an error message if it has encountered a problem with one file).

## Output variables

This is the list of the variables that are written in the output ROOT file:
- ``processID`` : processID of the event
- ``event`` : event number (for debugging)
- ``sqrtS`` : collision energy
- ``isValid_bb`` : is the event suitable for &nu;&nu;h (h &rarr; b bbar), set to false when less than 2 jets reconstructed
- ``isValid_ww`` : is the event suitable for &nu;&nu;h (h &rarr; WW* &rarr; qqqq), set to false when less than 4 jets reconstructed
- ``visible_e`` : total visible energy
- ``nParticles`` : number of reconstructed particles
- ``nIsoLep`` : number of isolated leptons
- ``eIsoLep`` : total energy of all isolated leptons
- ``higgs_e, higgs_pt, higgs_m`` : energy, transverse impulsion and mass of the reconstructed higgs
- ``higgs_cosTheta`` : cosine of angle of the higgs (with respect to z axis)
- ``higgs_recMass`` : recoil mass against the reconstructed higgs (is set to 0 if the value is impossible (square root of a negative number))
- ``higgs_bTag1`` : b-tagging of the most energetic jet (for 2-jet reconstruction)
- ``higgs_bTag2`` : b-tagging of the least energetic jet (for 2-jet reconstruction)
- ``b1_e, b1_pt, b1_m`` : energy, transverse impulsion and mass of the most energetic b quark (for 2-jet reconstruction)
- ``b2_e, b2_pt, b2_m`` : energy, transverse impulsion and mass of the least energetic b quark (for 2-jet reconstruction)
- ``w1_e, w1_pt, w1_m`` : energy, transverse impulsion and mass of the on-shell W boson (for 4-jet reconstruction)
- ``w2_e, w2_pt, w2_m`` : energy, transverse impulsion and mass of the off-shell W* boson (for 4-jet reconstruction)
- ``w1_cosBetw`` : cosine of the angle between the two jets idenfified as the on-shell W boson (for 4-jet reconstruction)
- ``w2_cosBetw`` : cosine of the angle between the two jets idenfified as the off-shell W* boson (for 4-jet reconstruction)
- ``higgs_bb_cosBetw`` : cosine of the angle between the two reconstructed b quarks (for 2-jet reconstruction)
- ``higgs_ww_cosBetw`` : cosine of the angle between the two reconstructed W bosons (for 4-jet reconstruction)
- ``y_12, y_23, y_34, y_45, y_56, y_67`` : jet parameters (actually -log10(y_ij))
- ``zz_z1_m`` : mass of the most energetic reconstructed Z boson (4-jet reconstruction : used to supress the ZZ background)
- ``zz_z2_m`` : mass of the least energetic reconstructed Z boson (4-jet reconstruction : used to supress the ZZ background)
- ``sl_w_m`` : mass of the reconstructed W boson (3-jet reconstruction : used to supress the WW semi-leptonic background)
- ``sl_rec_m`` : recoil mass against the reconstructed W boson (3-jet reconstruction : used to supress the WW semi-leptonic background)
- ``oblateness, sphericity, cosThrust, principleThrust, majorThrust, minorThrust`` : event shapes variables

These following variables are debug variables and should not be used for cuts of feed the BDT (This is true MC-information) :
- ``mc_ISR_e`` : total energy of ISR photons
- ``mc_ISR_pt`` : total transverse impulsion of ISR photons
- ``mc_nu_flavor`` : flavor of the two neutrinos
- ``mc_nu_e, mc_nu_pt, mc_nu_m`` : energy, transverse impulsion and mass of the two neutrinos system
- ``mc_nu_cosBetw`` : cosine of the angle between the two neutrinos
- ``mc_higgs_e, mc_higgs_pt, mc_higgs_m`` : true energy, transverse impulsion and mass of the higgs
- ``mc_higgs_recMass`` : true recoil mass against the higgs (is set to 0 if the value is impossible (square root of a negative number))
- ``mc_higgs_decay, mc_higgs_subDecay`` : see section [Higgs decay](#higgs-decay)
- ``mc_higgs_decay1_e, mc_higgs_decay1_pt, mc_higgs_decay1_m`` : true energy, transverse impulsion and mass of the most energetic decay from the higgs
- ``mc_higgs_decay2_e, mc_higgs_decay2_pt, mc_higgs_decay2_m`` : true energy, transverse impulsion and mass of the least energetic decay from the higgs
- ``mc_higgs_decay_cosBetw`` : cosine of the true angle between the two decays of the higgs

### Higgs Decay

The decay of the higgs is identified by 2 variables : ``mc_higgs_decay`` and ``mc_higgs_subDecay``.

The ``mc_higgs_decay`` variable is set to the [PDG code](https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf) of the decay, for example for h &rarr; b bbar : ``mc_higgs_decay == 5``, for h &rarr; &tau;&tau; : ``mc_higgs_decay == 15``, for h &rarr; ZZ* : ``mc_higgs_decay == 23``, etc... **Exception : for h &rarr; Z&gamma; it is set to 25.**

The ``mc_higgs_subDecay`` variable is there to represent the subdecay in h &rarr; WW* and h &rarr; ZZ* cases. 

- for WW* &rarr; qqqq and ZZ* &rarr; qqqq, ``mc_higgs_subDecay == 1``
- for ZZ* &rarr; qqll, ``mc_higgs_subDecay == 21``, ``22`` or ``23``, the second digit depending on the lepton flavor (1 for electron, 2 for muon and 3 for tau)
- for WW* &rarr; qql&nu;, ``mc_higgs_subDecay == 31``, ``32`` or ``33`` (same lepton flavor logic)
- for ZZ* &rarr; qq&nu;&nu;, ``mc_higgs_subDecay == 4`` (I did not care about the neutrino flavor here)
- for ZZ* &rarr; llll, ``mc_higgs_subDecay == 511``, ``512``, ``513``, ``521``, etc... (second and third digit corresponding to the two present lepton flavors)
- for WW* &rarr; l&nu;l&nu;, ``mc_higgs_subDecay == 611``, ``612``, ``613``, ``621``, etc... (same lepton flavor logic)
- for ZZ* &rarr; &nu;&nu;&nu;&nu;, ``mc_higgs_subDecay == 7``

For decays other than h &rarr; WW* and h &rarr; ZZ*, ``mc_higgs_subDecay == 0``.

These two variables are useful to distinguish between signal and background:
- for the &nu;&nu;h (h &rarr; b bbar) study, only the events with ``mc_higgs_decay == 5`` are signal events
- for the  &nu;&nu;h (h &rarr; WW* &rarr; qqqq) study, only the events with ``mc_higgs_decay == 24 && mc_higgs_subDecay == 1`` are signal events



