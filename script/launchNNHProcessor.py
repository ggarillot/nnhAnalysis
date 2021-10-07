#!/usr/bin/env python

import os
import abc
import time
import threading
import argparse

import NNHProcessor


class BaseThread(threading.Thread):

    __metaclass__ = abc.ABCMeta

    nCores = 1
    threads = []
    lock = threading.Lock()

    def __init__(self):
        threading.Thread.__init__(self)

    @abc.abstractmethod
    def launch(self):
        pass

    def run(self):
        while True:
            with BaseThread.lock:
                nCurrentThreads = len(BaseThread.threads)

                if nCurrentThreads < BaseThread.nCores:
                    t = threading.Thread(target=self.launch)
                    t.start()
                    BaseThread.threads.append(self)
                    break

            time.sleep(1)

        t.join()
        with BaseThread.lock:
            BaseThread.threads.remove(self)


class ProcessorThread(BaseThread):

    def __init__(self, _files, _outputFileName, remotePath=None):
        BaseThread.__init__(self)
        self.files = _files
        self.outputFileName = _outputFileName
        self.remotePath = remotePath

    def launch(self):
        print(f'Launch processor Files : {self.files}, output : {self.outputFileName}')

        if self.remotePath:
            for file in self.files:
                os.system(f'gfal-copy {self.remotePath}/{file} .')

        params = NNHProcessor.Params()
        params.outputFileName = self.outputFileName
        params.maxRecordNumber = 0
        params.skip = 0

        NNHProcessor.launch(params, self.files)

        if self.remotePath:
            for file in self.files:
                os.system(f'rm -f {file}')


class HaddThread(BaseThread):

    def __init__(self, _files, _outputFileName):
        BaseThread.__init__(self)
        self.files = _files
        self.outputFileName = _outputFileName

    def launch(self):
        os.system(f'rm -f {self.outputFileName}')

        haddCmd = f'hadd {self.outputFileName}'

        for file in self.files:
            haddCmd = f'{haddCmd} {file}'

        print(haddCmd)
        os.system(haddCmd)

        for file in self.files:
            os.system(f'rm -f {file}')


def launchAnalysis(processID, filesDirectory, remote=False):
    print(f"Launch analysis of process {processID}...")

    filesDirectory = f'{filesDirectory}/{processID}'
    fileList = []
    if not remote:
        fileList = [f'{filesDirectory}/{file}' for file in os.listdir(filesDirectory) if f'{processID}' in file]
    else:
        os.system(f'gfal-ls {filesDirectory} >> {processID}.txt')

        with open(f'{processID}.txt') as file:
            lines = file.readlines()
            fileList = [line.strip('\n') for line in lines]

        os.system(f'rm {processID}.txt')

    fileList.sort()

    fileList = [fileList[0]]

    print(f'{processID} : {len(fileList)} files to process')

    filesToMerge = []

    threads = []
    for i, file in enumerate(fileList):
        outputFileName = f'{processID}-{i}.root'
        if not remote:
            p = ProcessorThread([file], outputFileName)
        else:
            p = ProcessorThread([file], outputFileName, filesDirectory)
        filesToMerge.append(outputFileName)
        p.start()
        threads.append(p)

    for p in threads:
        p.join()

    outputFileName = f'{processID}.root'

    t = HaddThread(filesToMerge, outputFileName)
    t.start()
    t.join()


if __name__ == "__main__":

    # os.environ["MARLIN_DLL"] = '/home/garillot/nnhAnalysis/lib/libnnhAnalysis.so'

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--ncores', help='Number of threads', required=False, default=8)
    parser.add_argument('-p', '--processes', help='ProcessIDs to analyse', required=False, nargs='+')
    parser.add_argument('-f', '--filesDirectory', help='Path of remote files', required=True)
    parser.add_argument('-r', '--remote', help='indicate that files need to be downloaded', action='store_true', default=False)
    args = vars(parser.parse_args())

    BaseThread.nCores = int(args['ncores'])

    processesID = [402007, 402008, 402176, 402185, 402009, 402010, 402011, 402012, 402001, 402002, 402013, 402014, 402003, 402004,
                   402005, 402006, 500006, 500008, 500010, 500012, 500062, 500064, 500066, 500068, 500070, 500072, 500074, 500076,
                   500078, 500080, 500082, 500084, 500101, 500102, 500103, 500104, 500105, 500106, 500107, 500108, 500110, 500112,
                   500086, 500088, 500090, 500092, 500094, 500096, 500098, 500100, 500113, 500114, 500115, 500116, 500117, 500118,
                   500119, 500120, 500122, 500124, 500125, 500126, 500127, 500128]

    if args['processes']:
        processesID = []
        for p in args['processes']:
            processesID.append(p)
    processesID.sort()

    filesDirectory = args['filesDirectory']
    remote = False
    if args['remote']:
        remote = True

    threads = []
    for processID in processesID:
        p = threading.Thread(target=launchAnalysis, args=(processID, filesDirectory, remote))
        p.start()
        threads.append(p)

    for p in threads:
        p.join()
