#!/usr/bin/env python

import os


class Params:
    def __init__(self):
        self.outputFileName = 'test.root'
        self.maxRecordNumber = 0
        self.skip = 0


def launch(a, files):

    fileList = ''
    for name in files:
        fileList += name + ' '

    pid = os.getpid()

    xmlFileName = str(pid) + '.xml'

    xml = '''
<marlin>

	<execute>
		<processor name="MyProcessor_BG"/>
	</execute>

	<global>
		<parameter name="LCIOInputFiles">''' + fileList + '''</parameter>
		<parameter name="MaxRecordNumber" value="''' + str(
        a.maxRecordNumber) + '''"/>
		<parameter name="SkipNEvents" value="''' + str(
        a.skip) + '''" />
		<parameter name="SupressCheck" value="false" />
		<parameter name="Verbosity" options="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT"> MESSAGE </parameter> 
 	</global>

 	<processor name="MyProcessor_BG" type="MyProcessor_BG">
		<parameter name="RootFileName" type="string" >''' + a.outputFileName + '''</parameter>
	</processor>

</marlin>'''

    xmlFile = open(xmlFileName, 'w')
    xmlFile.write(xml)
    xmlFile.close()

    os.system('Marlin ' + xmlFileName)
    os.system('rm ' + xmlFileName)
