#!/usr/bin/env python

import os
import random
import string
import subprocess


class Params:
    def __init__(self):
        self.sqrtS = 250
        self.outputFileName = 'test.root'
        self.maxRecordNumber = 0
        self.skip = 0


def launch(params, files, logFileName=None):

    fileList = ''
    for name in files:
        fileList = f'{fileList}{name} '

    steeringFileName = ''.join(random.choice(
        string.ascii_letters) for i in range(10))

    xmlFileName = f'{steeringFileName}.xml'

    xml = f'''
<marlin>

	<execute>
		<processor name="NNHProcessor"/>
	</execute>

	<global>
		<parameter name="LCIOInputFiles">{fileList}</parameter>
		<parameter name="MaxRecordNumber" value="{params.maxRecordNumber}"/>
		<parameter name="SkipNEvents" value="{params.skip}" />
		<parameter name="SupressCheck" value="false" />
		<parameter name="Verbosity" options="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT"> MESSAGE </parameter>
 	</global>

 	<processor name="NNHProcessor" type="NNHProcessor">
        <parameter name="sqrtZ" type="float">{params.sqrtS}</parameter>
		<parameter name="RootFileName" type="string" >{params.outputFileName}</parameter>
	</processor>

</marlin>'''

    xmlFile = open(xmlFileName, 'w')
    xmlFile.write(xml)
    xmlFile.close()

    marlin = subprocess.Popen(f'Marlin {xmlFileName}', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = marlin.communicate()

    os.system(f'rm {xmlFileName}')

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
