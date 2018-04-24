import json
import os
import sys
import pybam
from string import split
from eKLIPse_fct import *


#=================================================================================================================================================#
#=================================================================================================================================================#
#=================================================================================================================================================#

def Read_alignment(titleBam,dicoInit,lstError):
    try:
        pathSCjson = os.path.join(dicoInit['pathTmpDir'],titleBam+"_SC.json")
        pathSCfasta = os.path.join(dicoInit['pathTmpDir'],titleBam+".fasta")
        FASTA = open(pathSCfasta,'w')
        # Init dicoBam
        dicoBam = {}
        for pos in range(1,dicoInit["dicoGbk"]['refLength']+1,1): 
            dicoBam[pos] = { 'nb_reads_F':0, 'nb_reads_R':0,\
                             'nb_sc_reads_F':0, 'nb_sc_reads_R':0,\
                             'nb_sc_fasta_F':0, 'nb_sc_fasta_R':0 }
        #***** BROWSE READS & SEARCH SCR *****#
        # Switch to downsampled BAM if exist
        if dicoInit['dicoBam'][titleBam]['path_downsampling']!="" : dicoInit['dicoBam'][titleBam]['path'] = dicoInit['dicoBam'][titleBam]['path_downsampling']
        for alignment in pybam.read(dicoInit['dicoBam'][titleBam]['path']):
            if alignment.file_chromosomes[alignment.sam_refID]==dicoInit["dicoBam"][titleBam]['refName'] and alignment.sam_mapq>=dicoInit['minQ'] and not alignment.sam_cigar_string.__contains__("H"):# and not explain_sam_flags(alignment.sam_flag).__contains__("second in pair") and not explain_sam_flags(alignment.sam_flag).__contains__("supplementary"):
                #***** RETRIEVE positions tuple & lastMapped infos *****#
                positionsLstTuple,lastMappedPos,lastMappedPosRead = cigar_list_to_tuple(alignment.sam_cigar_list,alignment.sam_pos0)
                #***** FORWARD reads *****#
                if explain_sam_flags(alignment.sam_flag)=="" or explain_sam_flags(alignment.sam_flag).__contains__("mate reverse strand"):
                    # Count reads
                    for posTuple in positionsLstTuple:
                        try : dicoBam[posTuple[1]+1]['nb_reads_F']+=1
                        except: pass # None case
                    # Count softclipped (right soft-clipping)
                    length, operation = alignment.sam_cigar_list[len(alignment.sam_cigar_list)-1]
                    if operation=="S":
                        dicoBam[lastMappedPos]['nb_sc_reads_F']+=1
                        # Write to Fasta (apply filter) / not consider 'N'
                        if length-alignment.sam_seq[len(alignment.sam_seq)-length:].count("N")>=dicoInit['SCsize'] and lastMappedPos+length<=dicoInit["dicoGbk"]['refLength']:
                            nb_mapped_part = min(dicoInit["MappedPart"],lastMappedPosRead)
                            FASTA.write(">"+str(lastMappedPos)+"_"+str(nb_mapped_part)+"_scrF_"+alignment.sam_qname+"\n"+alignment.sam_seq[len(alignment.sam_seq)-length-nb_mapped_part:]+"\n")
                            dicoBam[lastMappedPos]['nb_sc_fasta_F']+=1
                #***** REVERSE reads *****#
                else:
                    # Count reads
                    for posTuple in positionsLstTuple:
                        try : dicoBam[posTuple[1]+1]['nb_reads_R']+=1
                        except: pass
                    # Count softclipped (left soft-clipping)
                    length, operation = alignment.sam_cigar_list[0]
                    if operation=="S":
                        dicoBam[alignment.sam_pos0+1]['nb_sc_reads_R']+=1
                        # Write to Fasta (apply filter)
                        if length-alignment.sam_seq[0:length].count("N")>=dicoInit['SCsize'] and alignment.sam_pos0+1-length>=0:
                            nb_mapped_part = min(dicoInit["MappedPart"],length)
                            FASTA.write(">"+str(alignment.sam_pos0+1)+"_"+str(nb_mapped_part)+"_scrR_"+alignment.sam_qname+"\n"+alignment.sam_seq[0:length+nb_mapped_part]+"\n")
                            dicoBam[alignment.sam_pos0+1]['nb_sc_fasta_R']+=1
        # CLOSE files
        FASTA.close()
        # WRITE .json results
        JSON = open(pathSCjson,'wb')
        JSON.write(json.dumps(dicoBam))
        JSON.close()
    except:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        lstError.append("ReadThread \""+titleBam+"\": "+str(exc_value)+" (line "+str(exc_traceback.tb_lineno)+")")

#=================================================================================================================================================#
#=================================================================================================================================================#
#=================================================================================================================================================#

def SC_blast(titleBam,dicoInit,lstError):
    try:
        pathBLASTjson = os.path.join(dicoInit['pathTmpDir'],titleBam+"_BLAST.json")
        pathSCfasta = os.path.join(dicoInit['pathTmpDir'],titleBam+".fasta")
        pathSCDELBLAST = os.path.join(dicoInit['pathTmpDir'],titleBam+"_blast.out")
        pathSCjson = os.path.join(dicoInit['pathTmpDir'],titleBam+"_SC.json")
        dicoDel = {}
        # Require sc json files to modify positionsRef
        dicoBam = load_json(pathSCjson)

        #***** LAUNCH BlastN *****# (if fasta not empty)
        if os.path.getsize(pathSCfasta)!=0:
            cmd = dicoInit['pathBlastN']+" -task blastn-short -query "+pathSCfasta+" -db "+os.path.join(dicoInit['pathTmpDir'],"ref.fa")+" -outfmt \"10 qseqid qstart qend sstart send qlen length\" -out "+pathSCDELBLAST+" -perc_identity "+str(dicoInit['blastIdThreshold'])+" -max_hsps 1 -gapopen "+str(dicoInit['blastGapOpen'])+" -gapextend "+str(dicoInit['blastGapExt'])
            os.system(cmd)
            #***** Read Blast output file *****#
            IN = open(pathSCDELBLAST,'r')
            lst_lines = split(IN.read(),"\n")
            IN.close()
            for line in lst_lines:
                if line!="":
                    # Blast line parser
                    splitLine = split(line,",")
                    splitQread = split(splitLine[0],"_")
                    SCposRead = int(splitQread[0])
                    SCmappedpart = int(splitQread[1])
                    readOrient = splitQread[2]
                    qstart = int(splitLine[1])
                    qend = int(splitLine[2])
                    sstart = int(splitLine[3])
                    send = int(splitLine[4])
                    qlen = int(splitLine[5])
                    alignLength = int(splitLine[6])
                    covPercent = min((float(alignLength+SCmappedpart)*100.0)/float(qlen),100.0)
                    # Find correspunding positions
                    if covPercent>=dicoInit['blastCovThreshold'] and sstart<send:
                        if readOrient=="scrF":
                            delStart = SCposRead-SCmappedpart+qstart-1
                            delEnd = sstart
                            delName = str(delStart)+"-"+str(delEnd)
                            # Filter deletions
                            if (delStart<delEnd and delEnd-delStart>=dicoInit['SCsize']) or (delStart>delEnd and delStart-delEnd>=dicoInit["mitosize"]):
                                if delName in dicoDel:
                                    dicoDel[delName]['scrF']['nbBlast']+=1
                                    dicoDel[delName]['scrF']['initial_SCposRead'].append(SCposRead)
                                    if send>dicoDel[delName]['scrF']['limit']: dicoDel[delName]['scrF']['limit'] = send
                                else: dicoDel[delName] = { 'scrF':{'nbBlast':1,'limit':send,'initial_SCposRead':[SCposRead]}, 'scrR':{'nbBlast':0,'limit':99999,'initial_SCposRead':[]}, 'freqF':0.0, 'freqR':0.0 }
                        else:
                            delStart = send-(SCmappedpart-qlen+alignLength)+1
                            delEnd = SCposRead
                            delName = str(delStart)+"-"+str(delEnd)
                            # Filter deletions
                            if (delStart<delEnd and delEnd-delStart>=dicoInit['SCsize']) or (delStart>delEnd and delStart-delEnd>=dicoInit["mitosize"]):
                                if delName in dicoDel:
                                    dicoDel[delName]['scrR']['nbBlast']+=1
                                    dicoDel[delName]['scrR']['initial_SCposRead'].append(SCposRead)
                                    if sstart<dicoDel[delName]['scrR']['limit']: dicoDel[delName]['scrR']['limit'] = sstart
                                else: dicoDel[delName] = { 'scrF':{'nbBlast':0,'limit':0,'initial_SCposRead':[]}, 'scrR':{'nbBlast':1,'limit':sstart,'initial_SCposRead':[SCposRead]}, 'freqF':0.0, 'freqR':0.0 }
   
        #***** WRITE .json results file *****#
        JSON = open(pathBLASTjson,'wb')
        JSON.write(json.JSONEncoder().encode(dicoDel))
        JSON.close()
    except: 
        exc_type, exc_value, exc_traceback = sys.exc_info()
        lstError.append("SCBlastThread \""+titleBam+"\": "+str(exc_value)+" (line "+str(exc_traceback.tb_lineno)+")")

#=================================================================================================================================================#
#=================================================================================================================================================#
#=================================================================================================================================================#

def deletionPrediction(titleBam,dicoInit,lstError):
    try:
        pathSCjson = os.path.join(dicoInit['pathTmpDir'],titleBam+"_SC.json")
        pathBLASTjson = os.path.join(dicoInit['pathTmpDir'],titleBam+"_BLAST.json")
        pathCumuljson = os.path.join(dicoInit['pathTmpDir'],titleBam+"_cumul.json")
       
        #***** LOAD JSON files *****#
        dicoBam = load_json(pathSCjson)
        dicoDel = load_json(pathBLASTjson)
       
        #**** SHIFT Deletions *****#
        # Block one position and apply shift
        dico_shift_del = {}
        for delName in dicoDel.keys():
            start = int(split(delName,"-")[0])
            end = int(split(delName,"-")[1])
            for delName2 in dicoDel.keys():
                if delName!=delName2:
                    start2 = int(split(delName2,"-")[0])
                    end2 = int(split(delName2,"-")[1])
                    # Same deletion start => shift end
                    if (start==start2 and (abs(end-end2)<=dicoInit["delShift"])) or (end==end2 and (abs(start-start2)<=dicoInit["delShift"])):
                        dico_shift_del[len(dico_shift_del)] = set([delName,delName2])
        # Merge set
        for key in dico_shift_del.keys():
            for key2 in dico_shift_del.keys():
                if key!=key2:
                    if len(dico_shift_del[key].intersection(dico_shift_del[key2]))>0:
                        dico_shift_del[key] = dico_shift_del[key].union(dico_shift_del[key2])
                        dico_shift_del[key2] = set()
        # Merge to one maximal blast position
        for key in dico_shift_del.keys():
            if len(dico_shift_del[key])>0:
                nbBlast_F = 0 ; nbBlast_R = 0
                lst_start_sc_fasta = [] ; lst_end_sc_fasta = []
                sum_start_sc_fasta = 0 ; sum_end_sc_fasta = 0
                maxBlast_F = 0 ; maxBlast_R = 0
                maxBlast_start = "" ; maxBlast_end = ""
                limit_max_start = -1 ; limit_max_end = -1
                initial_SCposRead_F = [] ; initial_SCposRead_R = []
                # Browse deletion
                for delName in dico_shift_del[key]:
                    start = split(delName,"-")[0]
                    end = split(delName,"-")[1]
                    try: delnbblastF = dicoDel[delName]['scrF']['nbBlast']
                    except: delnbblastF = 0
                    try: delnbblastR = dicoDel[delName]['scrR']['nbBlast']
                    except: delnbblastR = 0
                    # Merge number of Blast
                    nbBlast_F+=delnbblastF
                    nbBlast_R+=delnbblastR
                    # Merge initial SC positions
                    try: initial_SCposRead_F.extend(dicoDel[delName]['scrF']['initial_SCposRead'])
                    except: pass
                    try: initial_SCposRead_R.extend(dicoDel[delName]['scrR']['initial_SCposRead'])
                    except: pass
                    # Find max blast positions for start and end
                    if delnbblastF>maxBlast_F:
                        maxBlast_F = delnbblastF
                        maxBlast_start = start
                        limit_max_start = dicoDel[delName]['scrF']['limit']
                    if delnbblastR>maxBlast_R:
                        maxBlast_R = delnbblastR
                        maxBlast_end = end
                        limit_max_end = dicoDel[delName]['scrR']['limit']
                    # Del entry
                    try : del(dicoDel[delName])
                    except: pass
                dicoDel[maxBlast_start+"-"+maxBlast_end] = { 'scrF':{'nbBlast':nbBlast_F,'limit':limit_max_start,'initial_SCposRead':initial_SCposRead_F},\
                                                             'scrR':{'nbBlast':nbBlast_R,'limit':limit_max_end,'initial_SCposRead':initial_SCposRead_R},\
                                                             'sum_start_sc_fasta':sum_start_sc_fasta, 'sum_end_sc_fasta':sum_end_sc_fasta,\
                                                             'freqF':0.0, 'freqR':0.0 }

        #***** COMPUTE Frequencies *****#
        for delName in dicoDel.keys():
            # FILTERs 
            if dicoInit["bilateral"]==True and (dicoDel[delName]['scrF']['nbBlast']<dicoInit["minblast"] or dicoDel[delName]['scrR']['nbBlast']<dicoInit["minblast"]): del(dicoDel[delName])
            elif dicoInit["bilateral"]==False and dicoDel[delName]['scrF']['nbBlast']<dicoInit["minblast"] and dicoDel[delName]['scrR']['nbBlast']<dicoInit["minblast"]: del(dicoDel[delName])
            # Compute frequency
            else:
                start = int(split(delName,"-")[0])
                end = int(split(delName,"-")[1])
                nb_blast_F = float(dicoDel[delName]['scrF']['nbBlast'])
                nb_blast_R = float(dicoDel[delName]['scrR']['nbBlast'])
                nb_sc_fasta_F = 0 ; set_sc_fasta_pos_F = set()
                nb_sc_fasta_R = 0 ; set_sc_fasta_pos_R = set()
                # Check max total and soft-clipped reads in +-shifted pos
                lst_nb_sc_reads_F = []
                for init_pos in set(dicoDel[delName]['scrF']['initial_SCposRead']):
                    for i in range(init_pos-dicoInit["delShift"],init_pos+dicoInit["delShift"]+1,1):
                        lst_nb_sc_reads_F.append(dicoBam[str(i)]['nb_sc_reads_F'])
                        if not i in set_sc_fasta_pos_F : nb_sc_fasta_F+=dicoBam[str(i)]['nb_sc_fasta_F'] ; set_sc_fasta_pos_F.add(i)
                lst_nb_sc_reads_R = []
                for init_pos in set(dicoDel[delName]['scrR']['initial_SCposRead']):
                    for i in range(init_pos-dicoInit["delShift"],init_pos+dicoInit["delShift"]+1,1):
                        lst_nb_sc_reads_R.append(dicoBam[str(i)]['nb_sc_reads_R'])
                        if not i in set_sc_fasta_pos_R : nb_sc_fasta_R+=dicoBam[str(i)]['nb_sc_fasta_R'] ; set_sc_fasta_pos_R.add(i)
                max_occ_start = max(dicoDel[delName]['scrF']['initial_SCposRead'],key=dicoDel[delName]['scrF']['initial_SCposRead'].count)
                max_occ_end = max(dicoDel[delName]['scrR']['initial_SCposRead'],key=dicoDel[delName]['scrR']['initial_SCposRead'].count)
                nb_reads_F = float(dicoBam[str(max_occ_start)]['nb_reads_F'])
                nb_reads_R = float(dicoBam[str(max_occ_end)]['nb_reads_R'])
                nb_sc_reads_F = float(max(lst_nb_sc_reads_F))
                nb_sc_reads_R = float(max(lst_nb_sc_reads_R))
                dicoDel[delName]['freqF'] = (nb_sc_reads_F/nb_reads_F) * (nb_blast_F/nb_sc_fasta_F) *100.0
                dicoDel[delName]['freqR'] = (nb_sc_reads_R/nb_reads_R) * (nb_blast_R/nb_sc_fasta_R) *100.0
                dicoDel[delName]['depthF'] = nb_reads_F
                dicoDel[delName]['depthR'] = nb_reads_R

        #***** CUMULATIVE Frequency *****#
        dicoCumulFreq = {}
        for delPos in dicoDel.keys():
            start = int(split(delPos,"-")[0])
            end = int(split(delPos,"-")[1])
            if start<end: posRange = range(start,end+1,1)
            else: posRange = range(start,dicoInit["dicoGbk"]['refLength']+1,1) ; posRange.extend(range(1,end+1,1))
            for pos in posRange:
                if pos in dicoCumulFreq: dicoCumulFreq[pos]+=(dicoDel[delPos]['freqF']+dicoDel[delPos]['freqR'])/200.0
                else: dicoCumulFreq[pos] = (dicoDel[delPos]['freqF']+dicoDel[delPos]['freqR'])/200.0

        #***** WRITE .json results files *****#
        JSON = open(pathBLASTjson,'wb')
        JSON.write(json.JSONEncoder().encode(dicoDel))
        JSON.close()
        JSON = open(pathCumuljson,'wb')
        JSON.write(json.JSONEncoder().encode(dicoCumulFreq))
        JSON.close()

    except: 
        exc_type, exc_value, exc_traceback = sys.exc_info()
        lstError.append("DelPredThread \""+titleBam+"\": "+str(exc_value)+" (line "+str(exc_traceback.tb_lineno)+")")

#=================================================================================================================================================#
#=================================================================================================================================================#
#=================================================================================================================================================#

def create_results_table(dicoInit,lstError):
    # Open output files
    pathOutDel = os.path.join(dicoInit['pathOutputDir'],"eKLIPse_deletions.csv")
    pathOutGenes = os.path.join(dicoInit['pathOutputDir'],"eKLIPse_genes.csv")
    OUTDEL = open(pathOutDel,'w')
    OUTGENES = open(pathOutGenes,'w')
    # Write headers
    headerDel = "\"Title\";\"5' breakpoint\";\"3' breakpoint\";\"Freq\";\"Freq For\";\"Freq Rev\";\"5' Blast\";\"3' Blast\";\"5' Depth\";\"3' Depth\";\"Repetition\"\n"
    OUTDEL.write(headerDel)
    headerGenes = "\"Gene\";\"Start\";\"End\";\"Type\""
    for titleBam in dicoInit["dicoBam"].keys(): headerGenes = headerGenes+",\""+titleBam+"\""
    OUTGENES.write(headerGenes+"\n")

    #***** BROWSE input alignments *****#
    dico_max_gene = {}
    for titleBam in dicoInit["dicoBam"].keys():
        # Load json results files
        dicoDel = load_json(os.path.join(dicoInit['pathTmpDir'],titleBam+"_BLAST.json"))
        dicoCumulFreq = load_json(os.path.join(dicoInit['pathTmpDir'],titleBam+"_cumul.json"))

        #***** DELETIONS Results File *****#
        # Sort deletion by deletion start and stop
        dicoSortDel = {}
        for deletion in dicoDel.keys():
            start = int(split(deletion,"-")[0])
            end = int(split(deletion,"-")[1])
            if start in dicoSortDel: dicoSortDel[start][end] = deletion
            else: dicoSortDel[start] = { end:deletion }
        lst_start = dicoSortDel.keys() ; lst_start.sort()
        # Browse sorted deletions
        for start in lst_start:
            lst_end = dicoSortDel[start].keys() ; lst_end.sort()
            for end in lst_end:
                # Search repetition (right)
                repetition_right = "" ; cpt_right = 0
                while dicoInit["dicoGbk"]['refSeq'][start+cpt_right]==dicoInit["dicoGbk"]['refSeq'][end-1+cpt_right]:
                    repetition_right = repetition_right+dicoInit["dicoGbk"]['refSeq'][start+cpt_right]
                    cpt_right+=1
                # Search repetition (right)
                repetition_left = "" ; cpt_left = 0
                while dicoInit["dicoGbk"]['refSeq'][start-2-cpt_left]==dicoInit["dicoGbk"]['refSeq'][end-1-cpt_left]:
                    repetition_left = repetition_left+dicoInit["dicoGbk"]['refSeq'][start-1-cpt_left]
                    cpt_left+=1
                repetition_left = repetition_left[::-1]
                # Keep best one
                if cpt_right==cpt_left==0: repetition = "none"
                elif cpt_right>=cpt_left:
                    repetition = str(start+1)+"-"+repetition_right+"-"+str(start+cpt_right)+" <> "+str(end)+"-"+repetition_right+"-"+str(end-1+cpt_right)
                else:
                    repetition = str(start-cpt_left)+"-"+repetition_left+"-"+str(start-1)+" <> "+str(end+1-cpt_left)+"-"+repetition_left+"-"+str(end) 
                # Write deletion line
                deletion = dicoSortDel[start][end]
                delPercent = ((dicoDel[deletion]['freqF']+dicoDel[deletion]['freqR'])/2.0)
                ToWrite = "\""+titleBam+"\";\""+str(start)+"\";\""+str(end)+"\";\""+str(delPercent).replace(".",",")+\
                          "\";\""+str(dicoDel[deletion]['freqF']).replace(".",",")+"\";\""+str(dicoDel[deletion]['freqR']).replace(".",",")+\
                          "\";\""+str(dicoDel[deletion]['scrF']['nbBlast'])+"\";\""+str(dicoDel[deletion]['scrR']['nbBlast'])+"\";\""+str(dicoDel[delName]['depthF'])+"\";\""+str(dicoDel[delName]['depthR'])+"\";\""+repetition+"\"\n"
                OUTDEL.write(ToWrite)

        #***** GENES max per genes hashtable *****#
        for features in dicoInit["dicoGbk"]['lstGene']:
            try: dico_max_gene[features[0]][titleBam] = 0.0
            except:  dico_max_gene[features[0]] = { titleBam:0.0 }
            for pos in dicoCumulFreq.keys():
                if int(pos)>=features[1] and int(pos)<=features[2]: dico_max_gene[features[0]][titleBam] = max(dico_max_gene[features[0]][titleBam],dicoCumulFreq[pos])
            
    #***** GENES Results File *****#
    for features in dicoInit["dicoGbk"]['lstGene']:
        ToWrite = "\""+features[0]+"\";\""+str(features[1])+"\";\""+str(features[2])+"\";\""+features[3]+"\""
        for titleBam in dicoInit["dicoBam"].keys():
            ToWrite = ToWrite+";\""+str(dico_max_gene[features[0]][titleBam]*100.0).replace(".",",")+"\""
        OUTGENES.write(ToWrite+"\n")

    # Close Files
    OUTDEL.close()
    OUTGENES.close()
