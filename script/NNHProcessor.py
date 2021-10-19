import threading
import queue
import subprocess

import os

from Observer import Observer


class Params:
    def __init__(self):
        self.inputFileNames = []
        self.sqrtS = 250
        self.outputFileName = ''
        self.maxRecordNumber = 0
        self.skip = 0


class NNHProcessorThread(threading.Thread):

    params = queue.Queue()
    lock = threading.Lock()
    threads = []

    def __init__(self):
        super().__init__()

        self.observers = []
        with NNHProcessorThread.lock:
            NNHProcessorThread.threads.append(self)

    def closeThreads():
        for i in range(len(NNHProcessorThread.threads)):
            NNHProcessorThread.params.put(None)

    def addParams(params: Params):
        NNHProcessorThread.params.put(params)

    def attachObserver(self, observer: Observer):
        if observer not in self.observers:
            self.observers.append(observer)

    def detachObserver(self, observer: Observer):
        self.observers.remove(observer)

    def notify(self, fileName, msg):
        for observer in self.observers:
            observer.update(fileName, msg)

    def launch(self, params: Params):

        os.environ['MARLIN_DLL'] = f"{os.environ['NNH_HOME']}/lib/libnnhAnalysis.so"

        fileList = ''
        for fileName in params.inputFileNames:
            fileList = f'{fileList}{fileName} '

        marlinCmd = f'''Marlin $NNH_HOME/script/NNH_steer.xml \\
                        --global.LCIOInputFiles={fileList} \\
                        --global.MaxRecordNumber={params.maxRecordNumber} \\
                        --global.SkipNEvents={params.skip} \\
                        --NNHProcessor.sqrtZ={params.sqrtS} \\
                        --NNHProcessor.RootFileName={params.outputFileName}'''

        marlin = subprocess.Popen(marlinCmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stdout, stderr = marlin.communicate()

        resultMsg = 'NNH_SUCCESS'
        if stderr:
            print(f'ERROR for {fileList} :')
            print(stderr.decode("utf-8"))

            resultMsg = 'NNH_FAIL'

        self.notify(params.outputFileName, resultMsg)

    def run(self):
        while True:
            params = NNHProcessorThread.params.get()
            if params is None:
                break

            self.launch(params)
            NNHProcessorThread.params.task_done()
