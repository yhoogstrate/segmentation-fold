#!/usr/bin/env python

"""
@file scripts/energy-estimation-utility/energy_estimation_utility/FoldController.py

@date 2015-07-20

@author Youri Hoogstrate

@section LICENSE
<PRE>
segmentation-fold can predict RNA 2D structures including K-turns.
Copyright (C) 2012-2015 Youri Hoogstrate

This file is part of segmentation-fold and originally taken from
yh-kt-fold.

segmentation-fold is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

segmentation-fold is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
</PRE>
"""



import shlex,subprocess,re,warnings,os



class FoldController:
    version_re = re.compile("\[Version\][\n ]+segmentation-fold[ \t]+([^\-\n]+)([^\n]*)")
    
    def __init__(self,binary,xml_file,sequence,associated_segments):
        self.binary = binary
        self.xml_file = xml_file
        
        self.set_sequence(sequence)
        self.associated_segments = associated_segments
    
    def set_sequence(self,sequence):
        self.sequence = sequence.upper().replace(" ","").replace("'","").replace("T","U").replace("\"","")
    
    def write_segments(self):
        fh = open(self.xml_file,"w")
        
        for segment_name,segment in self.associated_segments.items():
            self.write_segment(segment_name,segment,fh)
        
        fh.close()
    
    def write_segment(self,segment_name,segment,fh):
        self.we = round(segment['e'],4)
        fh.write(('<?xml version="1.0" encoding="UTF-8"?>\n'
                  '<root>\n'
                  '	<segments>\n'
                  '		<segment>\n'
                  '			<id>'+segment_name+'</id>\n'
                  '			<sequence_5prime>'+segment['5']+'</sequence_5prime>\n'
                  '			<bonds          >'+segment['x']+'</bonds>\n'
                  '			<sequence_3prime>'+segment['3']+'</sequence_3prime>\n'
                  '			<energy>'+str(round(segment['e'],4))+'</energy>\n'
                  '			<directions>\n'
                  '				<five_prime>true</five_prime>\n'
                  '				<three_prime>true</three_prime>\n'
                  '			</directions>\n'
                  '		</segment>\n'
                  '	</segments>\n'
                  '</root>'))
    
    def fold(self):
        self.write_segments()
        
        f = self.run_folding()
        os.remove(self.xml_file)
        
        return f
    
    def run_folding(self):
        argv = [self.binary,'-s',self.sequence,'-x',self.xml_file,'-t','2']
        output,error = subprocess.Popen(argv,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
        
        error = error.strip()
        if len(error) > 0:
            print "Error: "+error
            exit(1)
        
        output = output.split("\n")
        de = self.findDE(output[0])
        
        return {'free_energy':de,'dot_bracket':output[2].strip()}
    
    def findDE(self,string):
        e = string.split('dE:',1)[1].lstrip().split(' ',1)[0].strip()
        return round(float(e),4)
    
    def get_version(self):
        argv = [self.binary,'--version']
        
        try:
            output,error = subprocess.Popen(argv,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
        except:
            raise EnvironmentError('Binary \''+self.binary+'\' not executable, are you sure segmentation-fold is installed?')
        
        version_match = re.search(self.version_re,output)
        if version_match:
            version = version_match.group(1)
            suffix = version_match.group(2)
            if version[0:4] != "1.3.":
                raise EnvironmentError('Installed version of segmentation-fold ('+str(version)+') is out of date.')
            
            if suffix.find("(debug)") != -1:
                warnings.warn("Installed version of segmentation-fold was compiled in debug mode which is significantly slower than the release mode. We advice to recompile in relesae mode.", RuntimeWarning)
            
            return version
        else:
            raise EnvironmentError('Running \''+self.binary+' --version\' does not give the version in the right format. Aare you sure segmentation-fold is installed and up to date?')
