import sys
import os
import glob
import json
from string import split


#=================================================================================================================================================#
#=================================================================================================================================================#
#=================================================================================================================================================#

class Highlighter:
    def __init__(self):
        self._msgTypes={'black':'\033[0;30m', 'gray':'\033[1;30m', 'blue':'\033[0;34m', 'bold blue':'\033[1;34m', 'green':'\033[0;32m', 'bold green':'\033[1;32m', 'cyan':'\033[0;36m', 'bold cyan':'\033[1;36m', 'red':'\033[0;31m', 'bold red':'\033[1;31m', 'purple':'\033[0;35m', 'bold purple':'\033[1;35m', 'yellow':'\033[0;33m', 'bold yellow':'\033[1;33m', 'white':'\033[0;37m', 'bold white':'\033[1;37m'} ; self._reset='\033[0m' ; self._default='white'
    def ColorMsg(self,msg,msgLevel='white'):
        try:s=self._msgTypes[msgLevel]+msg+self._reset
        except: s=s=self._msgTypes[self._default]+msg+self._reset
        return s

def ColorOutput(msg,msgLevel='white'):
    return Highlighter().ColorMsg(msg,msgLevel) 

def printcolor(ToWrite,color,bool):
    if bool==True: sys.stdout.write(ColorOutput(ToWrite,color))
    else: sys.stdout.write(ToWrite)
    sys.stdout.flush()

def printerror(lstError):
    if len(lstError)!=0:
        for error in lstError: sys.stderr.write("\n[ERROR] "+error+"\n")
        exit("\n")

#=================================================================================================================================================#
#=================================================================================================================================================#
#=================================================================================================================================================#

def string_to_num(num):
    try: a = int(num) ; return "int"
    except ValueError:
        try: a = float(num) ; return "float"
        except ValueError: return "str"

def load_json(path):
    # Read json file
    json_data=open(path)
    data = json.load(json_data)
    json_data.close()
    return byteify(data)

def byteify(input):
    if isinstance(input, dict):
        return {byteify(key): byteify(value)
                for key, value in input.iteritems()}
    elif isinstance(input, list):
        return [byteify(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input


def explain_sam_flags(iFlags):
    lstFlags = [("read paired", 0x1),("read mapped in proper pair", 0x2),("read unmapped", 0x4),("mate unmapped", 0x8),("read reverse strand", 0x10),("mate reverse strand", 0x20),("first in pair", 0x40),("second in pair", 0x80),("not primary alignment", 0x100),("read fails platform/vendor quality checks", 0x200),("read is PCR or optical duplicate", 0x400),("supplementary alignment", 0x800)]
    decodeflag = ""
    for strFlagName, iMask in lstFlags:
        if iFlags & iMask:
            decodeflag = decodeflag+":"+strFlagName
    return decodeflag


def cigar_list_to_tuple(cigar_list,align_start):
    positionsLstTuple = []
    posRead = 0
    posmatch = 0
    posRef = align_start
    for cigar in cigar_list: #[(length,type), (length,type), (length,type)]
        for i in range(cigar[0]):
            if cigar[1]=="M":
                positionsLstTuple.append([posRead,posRef])
                posRead+=1
                posRef+=1
                posmatch+=1
            elif cigar[1]=="D":
                positionsLstTuple.append([None,posRef])
                posRef+=1
            else: # "I" & "S"
                positionsLstTuple.append([posRead,None])
                posRead+=1
    return positionsLstTuple,posRef,posmatch

#=================================================================================================================================================#
#=================================================================================================================================================#
#=================================================================================================================================================#

def mean_coverage(dicoInit,lstError):
    try:
        pathMeanCov = os.path.join(dicoInit['pathTmpDir'],"mean.cov")
        dicoMean = {}
        nb_sample = 0
        #***** BROWSE JSON files *****#
        for pathJson in glob.glob(os.path.join(dicoInit['pathTmpDir'],"*_SC.json")):
            dicoBam = load_json(pathJson)
            nb_sample+=1
            for key in dicoBam:
                try: dicoMean[key]+=dicoBam[key]['nb_reads_F']+dicoBam[key]['nb_reads_R']
                except: dicoMean[key] = dicoBam[key]['nb_reads_F']+dicoBam[key]['nb_reads_R']
        # Write Mean coverage file
        MEAN = open(pathMeanCov,'w')
        for pos in range(dicoInit["dicoGbk"]['refLength']-1):
            if str(pos) in dicoMean: ToWrite = "hsM\t"+str(pos)+"\t"+str(pos+1)+"\t"+str(float(dicoMean[str(pos)])/float(nb_sample))+"\n"
            else: ToWrite = "hsM\t"+str(pos)+"\t"+str(pos+1)+"\t0.0\n"
            MEAN.write(ToWrite)
        MEAN.close()
    except: 
        exc_type, exc_value, exc_traceback = sys.exc_info()
        lstError.append("mean_coverage \": "+str(exc_value)+" (line "+str(exc_traceback.tb_lineno)+")")
        
#=================================================================================================================================================#
#=================================================================================================================================================#
#=================================================================================================================================================#

def Alignment_downsampling(titleBam,dicoInit,lstError):
    try:
        dicoInit['dicoBam'][titleBam]['path_downsampling'] = os.path.join(dicoInit['pathTmpDir'],titleBam+"_reduce.bam")
        fraction = round((float(dicoInit['downCov'])/float(dicoInit['dicoBam'][titleBam]['nbReads'])),2)
        cmd = dicoInit['pathSamtools']+" view -b -q "+str(dicoInit["minQ"])+" -s "+str(fraction)+" "+dicoInit['dicoBam'][titleBam]['path']+" > "+dicoInit['dicoBam'][titleBam]['path_downsampling']
        os.system(cmd)
    except: 
        exc_type, exc_value, exc_traceback = sys.exc_info()
        lstError.append("Downsampling for BAM title \""+titleBam+"\": "+str(exc_value)+" (line "+str(exc_traceback.tb_lineno)+")")
