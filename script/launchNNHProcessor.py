#!/usr/bin/env python
import argparse
import os
import subprocess
import threading

from Transfer import TransferThread
import NNHProcessor
from NNHProcessor import NNHProcessorThread
from Merge import MergeThread

from Observer import Observer


class AnalysisFlow(Observer):

    mergedFiles = {}
    individualFiles = {}
    filesToDownload = {}
    failedFiles = {}

    lock = threading.Lock()

    def __init__(self):
        self.event = threading.Event()
        self.threadsToHandle = []

    def initAnalysisFlow(nThreadsTransfer, nThreadsProcess, nThreadsMerge):
        analysis = AnalysisFlow()

        for _ in range(nThreadsTransfer):
            thread = TransferThread()
            analysis.threadsToHandle.append(thread)
        for _ in range(nThreadsProcess):
            thread = NNHProcessorThread()
            analysis.threadsToHandle.append(thread)
        for _ in range(nThreadsMerge):
            thread = MergeThread()
            analysis.threadsToHandle.append(thread)

        for thread in analysis.threadsToHandle:
            thread.attachObserver(analysis)
            thread.start()

        return analysis

    def close():
        TransferThread.closeThreads()
        NNHProcessorThread.closeThreads()
        MergeThread.closeThreads()

    def addFile(self, targetMergedFileName, fileNamesToBeMerged, outputDirectory, remoteDirectory=None):

        mergedFile = {'notProcessed': [], 'processed': [], 'failed': [], 'outputDir': outputDirectory}
        with AnalysisFlow.lock:
            for i, file in enumerate(fileNamesToBeMerged):
                outputIndivName = targetMergedFileName.replace('.root', f'-{i}.root')
                AnalysisFlow.filesToDownload[file] = {'rootFileName': outputIndivName}
                AnalysisFlow.individualFiles[outputIndivName] = {'mergeInto': targetMergedFileName, 'slcioFile': file}
                mergedFile['notProcessed'].append(outputIndivName)

                if remoteDirectory:
                    self.launchTransfer(file, remoteDirectory, '.')
                else:
                    self.launchNNHProcessor(file, outputIndivName)

            AnalysisFlow.mergedFiles[targetMergedFileName] = mergedFile

    def launchTransfer(self, fileName, sourcePath, destinationPath):
        TransferThread.addTransfer(fileName, sourcePath, destinationPath)

    def launchNNHProcessor(self, inputFileName, outputFileName):
        params = NNHProcessor.Params()
        params.inputFileNames = [inputFileName]
        params.outputFileName = outputFileName
        NNHProcessorThread.addParams(params)

    def launchMerge(self, outputFileName, fileNamesToMerge):
        MergeThread.addMergeTask(outputFileName, fileNamesToMerge)

    def update(self, fileName, msg):

        if msg.startswith('TRANSFER'):
            with AnalysisFlow.lock:
                file = AnalysisFlow.filesToDownload.pop(fileName, None)
                remainingFilesToDownload = len(AnalysisFlow.filesToDownload)

            if file:
                rootFileName = file['rootFileName']

                if msg.endswith('SUCCESS'):
                    self.launchNNHProcessor(fileName, rootFileName)

                elif msg.endswith('FAIL'):
                    print(f'ERROR : transfer of {fileName} failed')
                    os.system(f'rm -f {fileName}')
                    with AnalysisFlow.lock:
                        mergedFile = AnalysisFlow.individualFiles[rootFileName]['mergeInto']

                        AnalysisFlow.mergedFiles[mergedFile]['notProcessed'].remove(rootFileName)
                        AnalysisFlow.mergedFiles[mergedFile]['failed'].append(rootFileName)
                else:
                    print(f'ERROR : message unknown : {msg}')

                print(f'{remainingFilesToDownload} remaining files to download')

            else:
                print(f'ERROR : unknown file :{fileName}')

        elif msg.startswith('NNH'):

            with AnalysisFlow.lock:
                file = AnalysisFlow.individualFiles.pop(fileName, None)
                remainingFilesToAnalyse = len(AnalysisFlow.individualFiles)

            if file:
                os.system(f'rm - r {file["slcioFile"]}')

                mergeFile = file['mergeInto']
                with AnalysisFlow.lock:
                    AnalysisFlow.mergedFiles[mergeFile]['notProcessed'].remove(fileName)

                    if msg.endswith('SUCCESS'):
                        AnalysisFlow.mergedFiles[mergeFile]['processed'].append(fileName)

                    elif msg.endswith('FAIL'):
                        print(f'ERROR : analyse of {fileName} failed')
                        AnalysisFlow.mergedFiles[mergeFile]['failed'].append(fileName)

                    else:
                        print(f'ERROR : message unknown : {msg}')

                    if not AnalysisFlow.mergedFiles[mergeFile]['notProcessed']:
                        self.launchMerge(mergeFile, AnalysisFlow.mergedFiles[mergeFile]['processed'])

                print(f'{remainingFilesToAnalyse} remaining files to analyse')

            else:
                print(f'ERROR : unknown file :{fileName}')

        elif msg.startswith('MERGE'):
            with AnalysisFlow.lock:
                file = AnalysisFlow.mergedFiles.pop(fileName, None)
                remainingFilesToMerge = len(AnalysisFlow.mergedFiles)

            if file:

                for ifile in file['processed']:
                    os.system(f'rm -f {ifile}')
                for ifile in file['failed']:
                    os.system(f'rm -f {ifile}')

                if msg.endswith('SUCCESS'):
                    os.system(f'mv {fileName} {file["outputDir"]}')
                elif msg.endswith('FAIL'):
                    os.system(f'rm -f {fileName}')
                else:
                    print(f'ERROR : message unknown : {msg}')

                print(f'{remainingFilesToMerge} remaining files to merge')

            else:
                print(f'ERROR : unknown file :{fileName}')

        else:
            print(f'ERROR : unknown message : {msg}')

        with AnalysisFlow.lock:
            if not AnalysisFlow.mergedFiles:
                self.event.set()

    def launchAnalysis(self, processID, filesDirectory, outputDirectory, remote=False):

        filesDirectory = f'{filesDirectory}/{processID}'
        fileNameList = []

        if not remote:
            fileNameList = [f'{filesDirectory}/{file}' for file in os.listdir(filesDirectory) if f'{processID}' in file]
        else:
            listFiles = subprocess.Popen(f'gfal-ls {filesDirectory}', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            stdout, stderr = listFiles.communicate()

            if stderr:
                print(f'ERROR for {processID} :')
                print(stderr.decode("utf-8"))
                return

            stdout = stdout.decode("utf-8")
            stdout = stdout.split('\n')

            for line in stdout:
                line = line.rstrip()
                if f'{processID}' in line:
                    fileNameList.append(line)

        fileNameList.sort()

        # fileNameList = fileNameList[0:2]

        print(f'{processID} : {len(fileNameList)} files to process')

        if not remote:
            self.addFile(f'{processID}.root', fileNameList, outputDirectory)
        else:
            self.addFile(f'{processID}.root', fileNameList, outputDirectory, filesDirectory)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--ncores', help='Number of threads', required=False, default=8)
    parser.add_argument('-p', '--processes', help='ProcessIDs to analyse', required=False, nargs='+')
    parser.add_argument('-f', '--filesDirectory', help='Path of remote files', required=True)
    parser.add_argument('-r', '--remote', help='indicate that files need to be downloaded', action='store_true', default=False)
    parser.add_argument('-o', '--outputDirectory', help='output directory', required=True)
    parser.add_argument('-l', '--log', help='output in text files instead of terminal', action='store_true', default=False)
    args = vars(parser.parse_args())

    nCores = int(args['ncores'])

    analysis = AnalysisFlow.initAnalysisFlow(1, nCores, 1)

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

    outputDirectory = args['outputDirectory']

    for processID in processesID:
        analysis.launchAnalysis(processID, filesDirectory, outputDirectory, remote)

    analysis.event.wait()
    AnalysisFlow.close()
