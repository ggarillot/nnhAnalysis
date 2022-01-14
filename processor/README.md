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

In order for the script to run properly, the input files must be ordered like this:

```
inputPath
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

