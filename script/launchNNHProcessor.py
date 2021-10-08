#!/usr/bin/env python

import os
import abc
import time
import threading
import argparse
import subprocess

import NNHProcessor


class BaseThread(threading.Thread):

    __metaclass__ = abc.ABCMeta

    nCores = 1
    threads = []
    waitingThreads = []
    lock = threading.Lock()

    def __init__(self):
        threading.Thread.__init__(self)

    @abc.abstractmethod
    def launch(self):
        pass

    def run(self):
        with BaseThread.lock:
            BaseThread.waitingThreads.append(self)

        while True:
            with BaseThread.lock:
                nCurrentThreads = len(BaseThread.threads)

                if nCurrentThreads < BaseThread.nCores:
                    t = threading.Thread(target=self.launch)
                    t.start()
                    BaseThread.threads.append(self)
                    BaseThread.waitingThreads.remove(self)
                    break

            time.sleep(1)

        t.join()
        with BaseThread.lock:
            BaseThread.threads.remove(self)
            print(f'{len(BaseThread.waitingThreads)} threads remaining...')


class ProcessorThread(BaseThread):

    def __init__(self, _files, _outputFileName, _logFileName=None, remotePath=None):
        BaseThread.__init__(self)
        self.files = _files
        self.outputFileName = _outputFileName
        self.logFileName = _logFileName
        self.remotePath = remotePath

    def downloadFiles(self):

        okFiles = []
        for file in self.files:

            if self.logFileName:
                with open(self.logFileName, 'w') as logFile:
                    logFile.write(f'Download {self.remotePath}/{file}...\n')
            else:
                print(f'Download {self.remotePath}/{file}...')

            download = subprocess.Popen(f'gfal-copy {self.remotePath}/{file} .', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            stdout, stderr = download.communicate()

            if stderr:
                print(f'ERROR for {file} :')
                print(stderr.decode("utf-8"))

                if self.logFileName:
                    with open(self.logFileName, 'a+') as logFile:
                        logFile.write(stderr.decode('utf-8'))
            else:
                okFiles.append(file)
                if self.logFileName:
                    with open(self.logFileName, 'w') as logFile:
                        logFile.write(f'Download of {self.remotePath}/{file} complete\n')
                else:
                    print(f'Download of {self.remotePath}/{file} complete')

        return okFiles

    def launch(self):
        # print(f'Launch processor for files : {self.files}, to output : {self.outputFileName}')

        okFiles = self.files
        if self.remotePath:
            okFiles = self.downloadFiles()

        params = NNHProcessor.Params()
        params.outputFileName = self.outputFileName
        params.maxRecordNumber = 0
        params.skip = 0

        success = NNHProcessor.launch(params, okFiles, self.logFileName)

        if self.remotePath:
            for file in self.files:
                os.system(f'rm -f {file}')

        # print(f'Processor for files : {self.files}, to output : {self.outputFileName} completed')


class MergeThread(BaseThread):

    def __init__(self, _files, _logFiles, _outputFileNameBase, _outputDirectory):
        BaseThread.__init__(self)
        self.files = _files
        self.logFiles = _logFiles
        self.outputFileNameBase = _outputFileNameBase
        self.outputDirectory = _outputDirectory

    def launch(self):
        os.system(f'rm -f {self.outputFileNameBase}.root')
        os.system(f'rm -f {self.outputFileNameBase}.txt')

        currentFiles = os.listdir()

        processedFiles = []

        for file in self.files:
            if file in currentFiles:
                processedFiles.append(file)
            else:
                self.outputFileNameBase = f'incomplete_{self.outputFileNameBase}.root'

        if self.logFiles:
            with open(f'{self.outputFileNameBase}.txt', 'w') as mergedlogFile:

                self.logFiles.sort()
                for singlelogFileName in self.logFiles:

                    fileNameBase = singlelogFileName.split('.txt')[0]
                    mergedlogFile.write(f'{fileNameBase} : ==================================================================== \n\n')

                    with open(singlelogFileName, 'r') as singleLogFile:
                        for line in singleLogFile.readlines():
                            mergedlogFile.write(line)
                        mergedlogFile.write('\n\n\n')

        haddCmd = f'hadd {self.outputFileNameBase}.root'

        for file in processedFiles:
            haddCmd = f'{haddCmd} {file}'

        if self.logFiles:
            with open(f'{self.outputFileNameBase}.txt', 'a+') as mergedlogFile:
                mergedlogFile.write('===========================================================\n\n')
                mergedlogFile.write(f'{haddCmd}\n')
            haddCmd = f'{haddCmd} >> {self.outputFileNameBase}.txt'
        else:
            print(haddCmd)

        hadd = subprocess.Popen(haddCmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stdout, stderr = hadd.communicate()

        if self.logFiles:
            with open(f'{self.outputFileNameBase}.txt', 'a+') as mergedlogFile:
                mergedlogFile.write(stdout.decode('utf-8'))
        else:
            print(stdout.decode('utf-8'))

        if stderr:
            print(stderr.decode("utf-8"))
            if self.logFiles:
                with open(f'{self.outputFileNameBase}.txt', 'a+') as mergedlogFile:
                    mergedlogFile.write(stderr.decode('utf-8'))

        for file in processedFiles:
            os.system(f'rm -f {file}')

        for file in self.logFiles:
            os.system(f'rm -f {file}')

        os.system(f'mkdir -p {self.outputDirectory}/log')

        os.system(f'mv {self.outputFileNameBase}.root {self.outputDirectory}')

        if self.logFiles:
            os.system(f'mv {self.outputFileNameBase}.txt {self.outputDirectory}/log')


def launchAnalysis(processID, filesDirectory, outputDirectory, remote=False, log=None):
    # print(f"Launch analysis of process {processID}...")

    filesDirectory = f'{filesDirectory}/{processID}'
    fileList = []
    if not remote:
        fileList = [f'{filesDirectory}/{file}' for file in os.listdir(filesDirectory) if f'{processID}' in file]
    else:
        listFiles = subprocess.Popen(f'gfal-ls {filesDirectory} >> {processID}.txt', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        _, stderr = listFiles.communicate()

        if stderr:
            print(f'ERROR for {processID} :')
            print(stderr.decode("utf-8"))
            return

        with open(f'{processID}.txt') as file:
            lines = file.readlines()
            fileList = [line.strip('\n') for line in lines]

        os.system(f'rm {processID}.txt')

    fileList.sort()

    fileList = fileList[0:2]

    print(f'{processID} : {len(fileList)} files to process')

    filesToMerge = []
    logFilesToMerge = []

    threads = []
    for i, file in enumerate(fileList):
        outputFileName = f'{processID}-{i}.root'

        logFile = None
        if log:
            logFile = f'{processID}-{i}.txt'
            logFilesToMerge.append(logFile)

        if not remote:
            p = ProcessorThread([file], outputFileName, logFile)
        else:
            p = ProcessorThread([file], outputFileName, logFile, filesDirectory)

        filesToMerge.append(outputFileName)
        p.start()
        threads.append(p)

    for p in threads:
        p.join()

    outputFileNameBase = f'{processID}'

    t = MergeThread(filesToMerge, logFilesToMerge, outputFileNameBase, outputDirectory)
    t.start()
    t.join()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--ncores', help='Number of threads', required=False, default=8)
    parser.add_argument('-p', '--processes', help='ProcessIDs to analyse', required=False, nargs='+')
    parser.add_argument('-f', '--filesDirectory', help='Path of remote files', required=True)
    parser.add_argument('-r', '--remote', help='indicate that files need to be downloaded', action='store_true', default=False)
    parser.add_argument('-o', '--output', help='output directory', required=True)
    parser.add_argument('-l', '--log', help='output in text files instead of terminal', action='store_true', default=False)
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

    log = False
    if args['log']:
        log = True

    output = args['output']

    threads = []
    for processID in processesID:
        p = threading.Thread(target=launchAnalysis, args=(processID, filesDirectory, output, remote, log))
        p.start()
        threads.append(p)

    for p in threads:
        p.join()
