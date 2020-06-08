import sys
import os
import shutil
from string import split
from eKLIPse_fct import *


#=================================================================================================================================================#
#=================================================================================================================================================#
#=================================================================================================================================================#

def circosConf(titleBam,dicoInit,lstError):
    try:
        #***** PATHS *****#
        # Main conf Templates
        pathConfTemplate = os.path.join(dicoInit['pathDataDir'],"circos_template","circos.conf")
        # Copy static conf files
        shutil.copy(os.path.join(dicoInit['pathDataDir'],"circos_template","ideogram.conf"),os.path.join(dicoInit['pathTmpDir'],"ideogram.conf"))
        shutil.copy(os.path.join(dicoInit['pathDataDir'],"circos_template","ticks.conf"),os.path.join(dicoInit['pathTmpDir'],"ticks.conf"))
        shutil.copy(os.path.join(dicoInit['pathDataDir'],"circos_template","bands.conf"),os.path.join(dicoInit['pathTmpDir'],"bands.conf"))
        # Template created from gbk input file
        pathKaryotypeTemplate = os.path.join(dicoInit['pathTmpDir'],titleBam+"_karyotype.colors.mt.txt")
        pathTextBands = os.path.join(dicoInit['pathTmpDir'],titleBam+"_text.bands.txt")
        TEMPLATE_KARYOTYPE = open(pathKaryotypeTemplate,'w')
        TEMPLATE_TEXTBANDS = open(pathTextBands,'w')
        TEMPLATE_KARYOTYPE.write("chr - hsM "+titleBam+" 0 "+str(dicoInit["dicoGbk"]['refLength']-1)+" chrM\n")
        cpt = 1 ; previous_end = 0 ; previous_color = ""
        for gene in dicoInit["dicoGbk"]['lstGene']:
            # Assign color to gene type
            if gene[3]=="trna":
                if previous_color=="orange": color = "dorange"
                else: color = "orange"
            elif gene[3]=="rrna":
                if previous_color=="purple": color = "vdpurple"
                else: color = "purple"
            else:
                if previous_color=="green": color = "vdgreen"
                else: color = "green"
            previous_color = color
            # Create karyotype and bands lines
            if previous_end<gene[1]: TEMPLATE_KARYOTYPE.write("band hsM band"+str(cpt)+" band"+str(cpt)+" "+str(previous_end)+" "+str(gene[1])+" vdgrey\n") ; cpt+=1
            if previous_end>gene[1]: TEMPLATE_KARYOTYPE.write("band hsM band"+str(cpt)+" band"+str(cpt)+" "+str(previous_end+1)+" "+str(gene[2])+" "+color+"\n")
            else: TEMPLATE_KARYOTYPE.write("band hsM band"+str(cpt)+" band"+str(cpt)+" "+str(gene[1])+" "+str(gene[2])+" "+color+"\n")
            TEMPLATE_TEXTBANDS.write("hsM "+str(gene[1])+" "+str(gene[2])+" "+gene[0]+"\n")
            cpt+=1 ; previous_end = gene[2]
        TEMPLATE_KARYOTYPE.write("band hsM band"+str(cpt)+" band"+str(cpt)+" "+str(previous_end+1)+" "+str(dicoInit["dicoGbk"]['refLength']-1)+" vdgrey\n")
        # For human reference display replications origins
        if dicoInit["dicoGbk"]['refName']=="NC_012920":
            TEMPLATE_TEXTBANDS.write("hsM 407 408 *******OH\nhsM 5747 5748 *******OL")
        TEMPLATE_KARYOTYPE.close()
        TEMPLATE_TEXTBANDS.close()

        # CurrenBAM specific files
        pathConf = os.path.join(dicoInit['pathTmpDir'],titleBam+"_circos.conf")
        pathCov = os.path.join(dicoInit["pathTmpDir"],titleBam+".cov")
        pathMeanCov = os.path.join(dicoInit['pathTmpDir'],"mean.cov")
        pathSCcircos = os.path.join(dicoInit['pathTmpDir'],titleBam+".sc")
        pathSCblastCircos = os.path.join(dicoInit['pathTmpDir'],titleBam+"_blastCircos.txt")
        pathSCblastCircos1 = os.path.join(dicoInit['pathTmpDir'],titleBam+"_blastCircos1.txt")
        pathSCblastCircos2 = os.path.join(dicoInit['pathTmpDir'],titleBam+"_blastCircos2.txt")
        pathSCblastCircos3 = os.path.join(dicoInit['pathTmpDir'],titleBam+"_blastCircos3.txt")
        nbLineBlastLink = 0 ; nbLineBlastLink1 = 0 ; nbLineBlastLink2 = 0 ; nbLineBlastLink3 = 0
        pathSCdelPred = os.path.join(dicoInit['pathTmpDir'],titleBam+"_SCdel.txt")
        pathDelCumul = os.path.join(dicoInit['pathTmpDir'],titleBam+"_cumulFreq.txt")

        #***** LOAD JSON files *****#
        dicoBam = load_json(os.path.join(dicoInit['pathTmpDir'],titleBam+"_SC.json"))
        dicoDel = load_json(os.path.join(dicoInit['pathTmpDir'],titleBam+"_BLAST.json"))

        #***** COVERAGE (sample) *****#
        COV = open(pathCov,'w')
        for pos in range(0,dicoInit["dicoGbk"]['refLength']-2,1):
            if str(pos+1) in dicoBam: ToWrite = "hsM\t"+str(pos+1)+"\t"+str(pos+2)+"\t"+str(dicoBam[str(pos+1)]['nb_reads_F']+dicoBam[str(pos+1)]['nb_reads_R'])+"\n"
            else: ToWrite = "hsM\t"+str(pos+1)+"\t"+str(pos+2)+"\t0\n"
            COV.write(ToWrite)
        COV.close()

        #***** SOFTCLIPPING (frequency peak) *****#
        CIRCOS = open(pathSCcircos,'w')
        for pos in range(dicoInit["dicoGbk"]['refLength']-1):
            freqSC = 0.0
            if str(pos+1) in dicoBam:
                if dicoBam[str(pos+1)]['nb_sc_reads_F']>dicoBam[str(pos+1)]['nb_sc_reads_R'] and dicoBam[str(pos+1)]['nb_reads_F']>0: freqSC = float(dicoBam[str(pos+1)]['nb_sc_reads_F'])/float(dicoBam[str(pos+1)]['nb_reads_F'])
                elif dicoBam[str(pos+1)]['nb_sc_reads_F']<dicoBam[str(pos+1)]['nb_sc_reads_R'] and dicoBam[str(pos+1)]['nb_reads_R']>0: freqSC = float(dicoBam[str(pos+1)]['nb_sc_reads_R'])/float(dicoBam[str(pos+1)]['nb_reads_R'])
            ToWrite = "hsM\t"+str(pos)+"\t"+str(pos+1)+"\t"+str(freqSC)+"\n"
            CIRCOS.write(ToWrite)
        CIRCOS.close()

        #***** BLAST Links *****#
        CIRCOS = open(pathSCblastCircos,'w')
        CIRCOS1 = open(pathSCblastCircos1,'w')
        CIRCOS2 = open(pathSCblastCircos2,'w')
        CIRCOS3 = open(pathSCblastCircos3,'w')        
        for key in dicoDel.keys():
            start = int(split(key,"-")[0])
            end = int(split(key,"-")[1])
            meanFreq = (dicoDel[key]['freqF']+dicoDel[key]['freqR'])/2.0
            if meanFreq >= 10:
                CIRCOS.write("hsM\t"+str(dicoDel[key]['scrR']['limit'])+"\t"+str(start)+"\thsM\t"+str(end)+"\t"+str(dicoDel[key]['scrF']['limit'])+"\n")
                nbLineBlastLink+=1
            elif meanFreq >= 1:
                CIRCOS1.write("hsM\t"+str(dicoDel[key]['scrR']['limit'])+"\t"+str(start)+"\thsM\t"+str(end)+"\t"+str(dicoDel[key]['scrF']['limit'])+"\n")
                nbLineBlastLink1+=1
            elif meanFreq >= 0.1:
                CIRCOS2.write("hsM\t"+str(dicoDel[key]['scrR']['limit'])+"\t"+str(start)+"\thsM\t"+str(end)+"\t"+str(dicoDel[key]['scrF']['limit'])+"\n")
                nbLineBlastLink2+=1
            else:
                CIRCOS3.write("hsM\t"+str(dicoDel[key]['scrR']['limit'])+"\t"+str(start)+"\thsM\t"+str(end)+"\t"+str(dicoDel[key]['scrF']['limit'])+"\n")
                nbLineBlastLink3+=1                
        CIRCOS.close()
        CIRCOS1.close()
        CIRCOS2.close()
        CIRCOS3.close()        
       
        #***** CUMULATIVE Frequency *****#
        dicoCumulFreq = load_json(os.path.join(dicoInit['pathTmpDir'],titleBam+"_cumul.json"))
        OUT = open(pathDelCumul,'w')
        for pos in range(dicoInit["dicoGbk"]['refLength']-1):
            if str(pos) in dicoCumulFreq:
                if dicoCumulFreq[str(pos)]>1.0: dicoCumulFreq[str(pos)] = 1.0
                OUT.write("hsM\t"+str(pos)+"\t"+str(pos+1)+"\t"+str(dicoCumulFreq[str(pos)])+"\n")
            else: OUT.write("hsM\t"+str(pos)+"\t"+str(pos+1)+"\t0\n")
        OUT.close()

        #***** CREATE Main Circos Configuration File *****#
        IN = open(pathConfTemplate,'r')
        lst_lines = split(IN.read(),"\n")
        IN.close()
        OUT = open(pathConf,'w')
        for line in lst_lines:
            if line.__contains__("#KARYOTYPE"): new_line = line.replace("#KARYOTYPE",pathKaryotypeTemplate)
            elif line.__contains__("#TEXTBANDS"): new_line = line.replace("#TEXTBANDS",pathTextBands)
            elif line.__contains__("#BLASTLINKS0"): new_line = line.replace("#BLASTLINKS0",pathSCblastCircos)
            elif line.__contains__("#LINKSNUMBER0"): new_line = line.replace("#LINKSNUMBER0",str(nbLineBlastLink)) 
            elif line.__contains__("#BLASTLINKS1"): new_line = line.replace("#BLASTLINKS1",pathSCblastCircos1)
            elif line.__contains__("#LINKSNUMBER1"): new_line = line.replace("#LINKSNUMBER1",str(nbLineBlastLink1))
            elif line.__contains__("#BLASTLINKS2"): new_line = line.replace("#BLASTLINKS2",pathSCblastCircos2)
            elif line.__contains__("#LINKSNUMBER2"): new_line = line.replace("#LINKSNUMBER2",str(nbLineBlastLink2))
            elif line.__contains__("#BLASTLINKS3"): new_line = line.replace("#BLASTLINKS3",pathSCblastCircos3)
            elif line.__contains__("#LINKSNUMBER3"): new_line = line.replace("#LINKSNUMBER3",str(nbLineBlastLink3))           
            elif line.__contains__("#SCFILE"): new_line = line.replace("#SCFILE",pathSCcircos)
            elif line.__contains__("#DATACOV"): new_line = line.replace("#DATACOV",pathCov)
            elif line.__contains__("#DATAMEDIANE"): new_line = line.replace("#DATAMEDIANE",pathMeanCov)
            elif line.__contains__("#DATAFREQ"): new_line = line.replace("#DATAFREQ",pathDelCumul)
            elif line.__contains__("#NAMEDIR"): new_line = line.replace("#NAMEDIR",dicoInit['pathOutputDir'])
            elif line.__contains__("#NAMEFILE"): new_line = line.replace("#NAMEFILE","eKLIPse_"+titleBam+".png")
            else: new_line = line
            OUT.write(new_line+"\n")
        OUT.close()
    except: 
        exc_type, exc_value, exc_traceback = sys.exc_info()
        lstError.append("ConfThread \""+titleBam+"\": "+str(exc_value)+" (line "+str(exc_traceback.tb_lineno)+")")

#=================================================================================================================================================#
#=================================================================================================================================================#
#=================================================================================================================================================#

def circosPlot(titleBam,dicoInit,lstError):
    try:
        pathConf = dicoInit['pathTmpDir']+"/"+titleBam+"_circos.conf"
        pathLog = os.path.join(dicoInit["pathTmpDir"],titleBam+"_circos.log")
        pathPng = os.path.join(dicoInit['pathOutputDir'],"eKLIPse_"+titleBam+".png")
        cmd = dicoInit['pathCircosBin']+" -nosvg -conf "+pathConf+" > "+pathLog+" 2>&1"
        os.system(cmd)
        if not os.path.isfile(pathPng):
            lstError.append("PlotThread \""+titleBam+"\": any output .png created (check log file)")
    except: 
        exc_type, exc_value, exc_traceback = sys.exc_info()
        lstError.append("PlotThread \""+titleBam+"\": "+str(exc_value)+" (line "+str(exc_traceback.tb_lineno)+")")