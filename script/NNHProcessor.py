#!/usr/bin/env python

import os
import subprocess


class Params:
    def __init__(self):
        self.sqrtS = 250
        self.outputFileName = 'test.root'
        self.maxRecordNumber = 0
        self.skip = 0


def launch(params, files, logFileName=None):

    os.environ['MARLIN_DLL'] = f"{os.environ['NNH_HOME']}/lib/libnnhAnalysis.so"

    fileList = ''
    for name in files:
        fileList = f'{fileList}{name} '

    marlinCmd = f'''Marlin $NNH_HOME/script/NNH_steer.xml \\
                    --global.LCIOInputFiles={fileList} \\
                    --global.MaxRecordNumber={params.maxRecordNumber} \\
                    --global.SkipNEvents={params.skip} \\
                    --NNHProcessor.sqrtZ={params.sqrtS} \\
                    --NNHProcessor.RootFileName={params.outputFileName}'''

    if logFileName:
        with open(logFileName, 'a+') as logFile:
            logFile.write(marlinCmd)
    else:
        print(marlinCmd)

    marlin = subprocess.Popen(marlinCmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = marlin.communicate()

    if logFileName:
        with open(logFileName, 'a+') as logFile:
            logFile.write(stdout.decode('utf-8'))
    else:
        print(stdout.decode('utf-8'))

    if stderr:
        print(f'ERROR for {files} :')
        print(stderr.decode("utf-8"))

        if logFileName:
            with open(logFileName, 'a+') as logFile:
                logFile.write(stderr.decode('utf-8'))
        return 1

    return 0
