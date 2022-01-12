# nnhAnalysis

This repository contains code for the analysis of the channels: 

- e+e- &rarr; &nu;&nu;h (h &rarr; WW &rarr; qqqq)
- e+e- &rarr; &nu;&nu;h (h &rarr; b bbar)

The following backgrounds are considered:

- 2 fermions (leptonic and hadronic)
- 4 fermions (leptonic, semi-leptonic and hadronic)
- ZH (higgsstrahlung) events with the Z decaying into (charged leptons, neutrinos and quarks)

Each channel is identified with an unique number. The corresponding numbers for each channel can be found [here](https://ild.ngt.ndu.ac.jp/elog/genmeta/).

### Signal Channels

One of the big problems we have is that there is no way to distinguish between the &nu;&nu;h events that comes from WW-fusion (what we are interested in) and the events that comes from Higgsstrahlung (+ interference). But we at least know that the WW-fusion events occurs only for eL.pR polarisation.

The signal events can be found in those three channels:
- 402007 : &nu;&nu;h events for all higgs decays (inclusive decay).
- 402173 : &nu;&nu;h events with exclusive b bbar decay
- 402176 : &nu;&nu;h events with exclusive WW decay

We are not interested in all the events that are present in those channels: in the 402007 we do not want all the higgs decays. And, for example, if you are studing the (h &rarr; b bbar) case, the channel 402176 becomes a background. The procedure to select only the relevent events is discussed [here](#higgs-decays)

### Background Channels

Higgs background: 
- 402008 : &nu;&nu;h events with inclusive h decay (Higgsstrahlung only)
- 402185 : &nu;&nu;h with exclusive WW decay (Higgsstrahlung only)
- 402009, 402010 : &nu;&nu;h events with inclusive h decay (Higgsstrahlung only)
- 402001 to 402006 : llh events
- 402011, 402012 : qqh events
- 402013, 402014 : eeh events (ZZ-fusion)

2-4 fermions background:
- 500006, 500008 : 2 fermions leptonic 
- 500010, 500012 : 2 fermions hadronic
- 500062, 500064, 500066, 500068, 500070, 500072 : 4 fermions hadronic
- 500086, 500088, 500090, 500092, 500094, 500096, 500098, 500100, 500113, 500114, 500115, 500116, 500117, 500118, 500119, 500120, 500122, 500124, 500125, 500126, 500127, 500128 : 4 fermions leptonic
- 500074, 500076, 500078, 500080, 500082, 500084, 500101, 500102, 500103, 500104, 500105, 500106, 500107, 500108, 500110, 500112 : 4 fermions semi-leptonic

## Analysis workflow

The files we are interested in are the DST files.
The first step of the analysis is to launch some processors to compute some event shape variables (thrust, sphericity...), identify isolated leptons and photons, apply the jet algorithm and launch the b-tagging on those reconstructed jets. All those tasks are done in one loop using the mini-DST format ([mini-DST](https://github.com/shkawada/mini-DST)). The ["official" steering code](https://github.com/shkawada/mini-DST/blob/master/mini-DST.xml) to produce mini-DST files will launch the jet algorithm to produce from 2 to 6 jets and apply the b-tagging on all of those jets. The b-tagging takes a lot of time, so to speed up I choose to only launch the jets algorithm for 2,3 and 4 jets, because I did not really need the 5 and 6 jets reconstruction. The steering file that I used is located in the folder (PUT FOLDER)
Be careful that some paths may need to be modified.

All the mini-DST files I produced are stored here :
root://lyogrid06.in2p3.fr/dpm/in2p3.fr/home/calice/garillot/ILD/AHCAL

Once you have the mini-DST files ready we can proceed with the analysis.
The analysis code is a bit rushed and not very friendly. Because I wanted to experiment with scikit-learn instead of the usual TMVA, this introduces an additional layer of complexity. Basically the workflow of the analysis is :
- launch a Marlin processor for each mini-DST file to compute all the event variables (higgs energy, mass, pt, higgs decay mass, number of particles, jet parameters,...). This will produces one root file per initial mini-DST file.
- combine all these small root files into a single gigantic root file with a TTree containing all the events (signal and background)
- create training and testing splits for BDT
- train a BDT using scikit-learn on training set, then apply the BDT and store the results in a TTree in another root file.
- launch some macros to plot the results.













