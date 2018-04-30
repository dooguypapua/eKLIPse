import os
import shutil
import uuid
import time
import tabulate
import pybam
from Bio import SeqIO
from eKLIPse_fct import *


#=================================================================================================================================================#
#=================================================================================================================================================#
#=================================================================================================================================================#

def arg_manager(argv,pathRootDir,boolColor,boolQT,spinner):
    struuid = str(uuid.uuid4())[:8] # analysis identifier
    #***** Initialization Dictionnary *****#
    dicoInit = {
                'startTime':time.time(),\
                'pathRootDir':os.path.dirname(os.path.realpath(__file__)),\
                'pathDataDir':os.path.join(os.path.dirname(os.path.realpath(__file__)),"data"),\
                'pathCircosTemplateDir':os.path.join(os.path.dirname(os.path.realpath(__file__)),"data","circos_template"),\
                'pathInput':"",\
                'pathOutputDir':os.path.join(os.getcwd(),"eKLIPse_"+struuid),\
                'pathGbkRef':"",\
                'pathTmpDir':os.path.join("/tmp","eKLIPse_"+struuid),\
                'pathBlastN':"blastn",\
                'pathMakeblastdb':"makeblastdb",\
                'pathCircosBin':"circos",\
                'pathSamtools':"samtools",\
                'dicoBam':{},\
                'dicoGbk':{},\
                'windows':False,\
                'lstTitleBam':[],\
                'minQ':20, 'minlen':100, 'SCsize':25, 'MappedPart':20,\
                'downCov':500000, 'delShift':5,\
                'mitosize':1000, 'minblast':1, 'bilateral':True,\
                'blastIdThreshold':80, 'blastCovThreshold':70,'blastGapOpen':0, 'blastGapExt':2,\
                'nbThread':1, 'boolColor':boolColor,'boolQT':boolQT,\
                }

    #***** Check Required Arguments *****#
    if len(argv)==1 or argv[1]=="-h" or argv[1]=="--help": manual_display([],dicoInit,spinner)
    if not "-in" in argv and not "--test" in argv: manual_display(["Missing required \"-in\" argument"],dicoInit,spinner)
    if not "-ref" in argv and not "--test" in argv: manual_display(["Missing required \"-ref\" argument"],dicoInit,spinner)
    # Executables path (local for Windows  / $PATH or arguments for linux)
    path_exe_src = os.path.join(os.path.dirname(os.path.realpath(__file__)),"src")
    if os.path.isdir(path_exe_src):
        dicoInit['pathBlastN'] = os.path.join(path_exe_src,"blast-2.6.0+","blastn.exe")
        dicoInit['pathMakeblastdb'] = os.path.join(path_exe_src,"blast-2.6.0+","makeblastdb.exe")
        dicoInit['pathCircosBin'] = os.path.join(path_exe_src,"circos-0.69-6","bin","circos.exe")
        dicoInit['pathSamtools'] = os.path.join(path_exe_src,"samtools.exe")
        dicoInit['windows'] = True
    else:
        if not "-blastn" in argv and any(os.access(os.path.join(path, "blastn"), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))==False: manual_display(["eKLIPse required \"blastn\" in PATH or \"-blastn\" defined"],dicoInit,spinner)
        if not "-makeblastdb" in argv and any(os.access(os.path.join(path, "makeblastdb"), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))==False: manual_display(["eKLIPse required \"makeblastdb\" in PATH or \"-makeblastdb\" defined"],dicoInit,spinner)
        if not "-circos" in argv and any(os.access(os.path.join(path, "circos"), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))==False: manual_display(["eKLIPse required \"circos\" in PATH or \"-circos\" defined"],dicoInit,spinner)
        if not "-samtools" in argv and any(os.access(os.path.join(path, "samtools"), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))==False: manual_display(["eKLIPse required \"samtools\" in PATH or \"-samtools\" defined"],dicoInit,spinner)

    #**** Read Passed Arguments *****#
    prev_arg = ""
    lstErrorDisplay = []
    for i in range(1,len(argv),2):
        if prev_arg=="--": i-=1
        prev_arg = ""
        # Input file with BAMs path
        if argv[i]=="-in":
            if not os.path.isfile(argv[i+1]): lstErrorDisplay.append("Input file \""+argv[i+1]+"\" not found")
            else:
                dicoInit['pathInput'] = argv[i+1]
                IN = open(argv[i+1],'r') ; lstLine = split(IN.read(),"\n") ; IN.close()
                for j in range(len(lstLine)):
                    if lstLine[j]!="" and not lstLine[j][0]=="#":
                        try:
                            pathBam = split(lstLine[j],"\t")[0]
                            titleBam = split(lstLine[j],"\t")[1]
                            if not os.path.isfile(pathBam): lstErrorDisplay.append("Input file at line"+str(j+1)+": \""+pathBam+"\" not found")
                            if titleBam in dicoInit['lstTitleBam']: lstErrorDisplay.append("Input file at line "+str(j+1)+": \""+titleBam+"\" name already used")
                            dicoInit['lstTitleBam'].append(titleBam)
                            dicoInit["dicoBam"][titleBam] = {'path':pathBam, 'refName':"", 'nbReads':0, 'path_downsampling':""}
                        except: lstErrorDisplay.append("Input file at line "+str(j+1)+": bad syntax")
        # Reference GBK file
        elif argv[i]=="-ref":
            if not os.path.isfile(argv[i+1]): lstErrorDisplay.append("Reference GBK file \""+argv[i+1]+"\" not found")
            else:
                dicoInit['pathGbkRef'] = argv[i+1]
                for gb_record in SeqIO.parse(open(argv[i+1],"r"), "genbank") :
                    dicoInit["dicoGbk"]['refName'] = gb_record.name
                    dicoInit["dicoGbk"]['refDescr'] = gb_record.description
                    dicoInit["dicoGbk"]['refSeq'] = str(gb_record.seq)
                    dicoInit["dicoGbk"]['refLength'] = len(gb_record.seq)
                    dicoInit["dicoGbk"]['lstGene'] = []
                    for feature in gb_record.features:
                        if feature.type=="gene":
                            if feature.qualifiers['nomenclature'][0].lower().__contains__("trna"): geneType = "trna"
                            elif feature.qualifiers['nomenclature'][0].lower().__contains__("rna"): geneType = "rrna"
                            else: geneType = "protein"
                            dicoInit["dicoGbk"]['lstGene'].append([feature.qualifiers['gene'][0],int(feature.location.start)+1,int(feature.location.end),geneType])
            if len(dicoInit["dicoGbk"])==0: lstErrorDisplay.append("Cannot Read Reference GBK file \""+argv[i+1]+"\"")
        elif argv[i]=="--test":
            prev_arg = "--"
            dicoInit['lstTitleBam'] = ["test_Illumina","test_proton"]
            dicoInit["dicoBam"] = {"test_Illumina":{'path':dicoInit["pathDataDir"]+"/test_illumina.bam", 'refName':"", 'nbReads':0, 'path_downsampling':""},"test_proton":{'path':dicoInit["pathDataDir"]+"/test_proton.bam", 'refName':"", 'nbReads':0, 'path_downsampling':""}}
            dicoInit['pathGbkRef'] = dicoInit["pathDataDir"]+"/NC_012920.1.gb"
            for gb_record in SeqIO.parse(open(dicoInit['pathGbkRef'],"r"), "genbank") :
                dicoInit["dicoGbk"]['refName'] = gb_record.name
                dicoInit["dicoGbk"]['refDescr'] = gb_record.description
                dicoInit["dicoGbk"]['refSeq'] = str(gb_record.seq)
                dicoInit["dicoGbk"]['refLength'] = len(gb_record.seq)
                dicoInit["dicoGbk"]['lstGene'] = []
                for feature in gb_record.features:
                    if feature.type=="gene":
                        if feature.qualifiers['nomenclature'][0].lower().__contains__("trna"): geneType = "trna"
                        elif feature.qualifiers['nomenclature'][0].lower().__contains__("rna"): geneType = "rrna"
                        else: geneType = "protein"
                        dicoInit["dicoGbk"]['lstGene'].append([feature.qualifiers['gene'][0],int(feature.location.start)+1,int(feature.location.end),geneType])
        # Results Directory
        elif argv[i]=="-out":
            if boolQT==True: dicoInit['pathOutputDir'] = argv[i+1] # uuid already done
            else: dicoInit['pathOutputDir'] = os.path.join(argv[i+1],"eKLIPse_"+struuid)
        # Temporary Directory
        elif argv[i]=="-tmp": # Temporary folder
            if boolQT==True: dicoInit['pathTmpDir'] = argv[i+1] # uuid already done
            else: dicoInit['pathTmpDir'] = os.path.join(argv[i+1],"eKLIPse_"+struuid)
        # Min Read quality
        elif argv[i]=="-minq": # The minimum read quality to consider
            if string_to_num(argv[i+1])!="int": lstErrorDisplay.append("\"-minQ\" integer expected")
            else: dicoInit["minQ"] = int(argv[i+1])
        # Min Read length
        elif argv[i]=="-minlen":
            if string_to_num(argv[i+1])!="int": lstErrorDisplay.append("\"-minlen\" integer expected")
            else: dicoInit["minlen"] = int(argv[i+1])
        # SoftClipping analysis threshold     
        elif argv[i]=="-scsize": # The minimum length of softclip to keep
            if string_to_num(argv[i+1])!="int": lstErrorDisplay.append("\"-SC_size\" integer expected")
            else: dicoInit["SCsize"] = int(argv[i+1])
        # SoftClipping analysis threshold     
        elif argv[i]=="-mapsize": # The size of mapped read sequence send to blast
            if string_to_num(argv[i+1])!="int": lstErrorDisplay.append("\"-MappedPart\" integer expected")
            else: dicoInit["MappedPart"] = int(argv[i+1])
        # Maximum reads coverage
        elif argv[i]=="-downcov": # 
            if string_to_num(argv[i+1])!="int": lstErrorDisplay.append("\"-downCov\" integer expected")
            else: dicoInit["downCov"] = int(argv[i+1])
        # Deletion position shift
        elif argv[i]=="-shift": # 
            if string_to_num(argv[i+1])!="int": lstErrorDisplay.append("\"-shift\" integer expected")
            else: dicoInit["delShift"] = int(argv[i+1])
        # Minimal length for deleted mito
        elif argv[i]=="-mitosize": # 
            if string_to_num(argv[i+1])!="int": lstErrorDisplay.append("\"-mitosize\" integer expected")
            else: dicoInit["mitosize"] = int(argv[i+1])
        # Minimal number of BLAST per breakpoint
        elif argv[i]=="-minblast": # 
            if string_to_num(argv[i+1])!="int": lstErrorDisplay.append("\"-minblast\" integer expected")
            else: dicoInit["minblast"] = int(argv[i+1])
        # Filter non-bilateral BLAST deletions
        elif argv[i]=="-bilateral": # 
            if argv[i+1].lower()!="true" and argv[i+1].lower()!="false": lstErrorDisplay.append("\"-bilateral\" boolean expected")
            else:
                if argv[i+1].lower()=="true": dicoInit["bilateral"] = True
                else: dicoInit["bilateral"] = False
        # BlastN analysis threshold
        elif argv[i]=="-id": # The minimum identity to keep a blast
            if string_to_num(argv[i+1])!="int": lstErrorDisplay.append("\"-blast_id\" integer expected")
            else: dicoInit["blastIdThreshold"] = int(argv[i+1])        
        elif argv[i]=="-cov": # The minimum depth between the query and the subject in the blast
            if string_to_num(argv[i+1])!="int": lstErrorDisplay.append("\"-blast_cov\" integer expected")
            else: dicoInit["blastCovThreshold"] = int(argv[i+1])
        elif argv[i]=="-gapopen": # Cost to open a gap
            if string_to_num(argv[i+1])!="int": lstErrorDisplay.append("\"-blast_gapopen\" integer expected")
            else: dicoInit["blastGapOpen"] = int(argv[i+1])        
        elif argv[i]=="-gapext": # Cost to extend a gap
            if string_to_num(argv[i+1])!="int": lstErrorDisplay.append("\"-blast_gapext\" integer expected")
            else: dicoInit["blastGapExt"] = int(argv[i+1])           
        # Number of threads
        elif argv[i]=="-thread":
            if string_to_num(argv[i+1])!="int": lstErrorDisplay.append("\"-nb_threads\" integer expected")
            else: dicoInit["nbThread"] = int(argv[i+1])
        # Programs
        elif argv[i]=="-blastn":
            if not os.path.isfile(argv[i+1]): lstErrorDisplay.append("\"-blastn\" binary not found")
            else: dicoInit['pathBlastN'] = argv[i+1]
        elif argv[i]=="-makeblastdb":
            if not os.path.isfile(argv[i+1]): lstErrorDisplay.append("\"-makeblastdb\" binary not found")
            else: dicoInit['pathMakeblastdb'] = argv[i+1]            
        elif argv[i]=="-circos":
            if not os.path.isfile(argv[i+1]): lstErrorDisplay.append("\"-circos\" binary not found")
            else: dicoInit['pathCircosBin'] = argv[i+1]
        elif argv[i]=="-samtools":
            if not os.path.isfile(argv[i+1]): lstErrorDisplay.append("\"-samtools\" binary not found")
            else: dicoInit['pathSamtools'] = argv[i+1]
        elif argv[i]=="--qtgui" or argv[i]=="--nocolor": prev_arg = "--"
        else: lstErrorDisplay.append("\""+argv[i]+"\" unrecognized parameters")

    #***** CHECK GAP OPEN & EXTEND *****#
    boolGap = False
    if dicoInit["blastGapOpen"]>=2 and dicoInit["blastGapExt"]>=2: boolGap = True
    if dicoInit["blastGapOpen"]==1 and dicoInit["blastGapExt"]==2: boolGap = True
    if dicoInit["blastGapOpen"]==0 and dicoInit["blastGapExt"]==2: boolGap = True
    if dicoInit["blastGapOpen"]==2 and dicoInit["blastGapExt"]==1: boolGap = True
    if dicoInit["blastGapOpen"]==1 and dicoInit["blastGapExt"]==1: boolGap = True
    if boolGap==False: lstErrorDisplay.append("Gap existence and extension "+str(dicoInit["blastGapOpen"])+" and "+str(dicoInit["blastGapExt"])+" are not supported. Please use (>=2 and >=2) or (1 and 2) or (0 and 2) or (2 and 1) or (1 and 1)")

    #***** CREATE DIRECTORIES *****#
    if not os.path.isdir(dicoInit['pathOutputDir']):
        try: os.makedirs(dicoInit['pathOutputDir'])
        except: lstErrorDisplay.append("Unable to create output folder \""+dicoInit['pathOutputDir']+"\"")
    if not os.path.isdir(dicoInit['pathTmpDir']):
        try: os.makedirs(dicoInit['pathTmpDir'])
        except: lstErrorDisplay.append("Unable to create temporary folder \""+dicoInit['pathTmpDir']+"\"")
    
    #***** CHECK errors & BAMs *****#
    # Return error(s) if found
    if len(lstErrorDisplay)!=0: manual_display(lstErrorDisplay,dicoInit,spinner)
    # Else check BAMs
    else:
        # Copy circos conf files
        headerBool = False
        for titleBam in dicoInit["dicoBam"].keys():
            # Retrieve header
            path_idxstats = os.path.join(dicoInit['pathTmpDir'],"idxstats.txt")
            if dicoInit['windows']==False: cmd_idxstats = dicoInit['pathSamtools']+" idxstats "+dicoInit['dicoBam'][titleBam]['path']+" 1>"+path_idxstats+" 2>/dev/null"
            else: cmd_idxstats = dicoInit['pathSamtools']+" idxstats "+dicoInit['dicoBam'][titleBam]['path']+" 1>"+path_idxstats+" 2> nul"
            os.system(cmd_idxstats)
            IN = open(path_idxstats,'r') ; lst_lines = split(IN.read(),"\n") ; IN.close()
            for line in lst_lines:
                split_line = split(line,"\t")
                if len(split_line)>1:
                    refName = split_line[0]
                    refSize = int(split_line[1])
                    nb_reads = int(split_line[2])
                    if refSize==dicoInit["dicoGbk"]['refLength']:
                        headerBool = True
                        dicoInit["dicoBam"][titleBam]['refName'] = refName
                        dicoInit["dicoBam"][titleBam]['nbReads'] = nb_reads
                        break
            if not headerBool: lstErrorDisplay.append("Any mitochondrial header found in file \""+dicoInit["dicoBam"][titleBam]['path']+"\"")
            if dicoInit["dicoBam"][titleBam]['refName']=="": lstErrorDisplay.append("Unable to find correspunding reference size in header \""+dicoInit["dicoBam"][titleBam]['path']+"\"")    
            
        # Return BAM error(s) if found
        if len(lstErrorDisplay)!=0: manual_display(lstErrorDisplay,dicoInit,spinner)
    
    #***** WRITE REFERENCE Fasta and blastdb *****#
    path_fasta_out = os.path.join(dicoInit['pathTmpDir'],"ref.fa")
    OUT = open(path_fasta_out,'w')
    OUT.write(">"+dicoInit["dicoGbk"]['refName']+"\n"+dicoInit["dicoGbk"]['refSeq']+"\n")
    OUT.close()
    if dicoInit['windows']==False: os.system(dicoInit['pathMakeblastdb']+" -dbtype nucl -in "+path_fasta_out+" > /dev/null 2>&1")
    else: os.system(dicoInit['pathMakeblastdb']+" -dbtype nucl -in "+path_fasta_out+" > nul 2>&1")

    # Return dicoInit
    return dicoInit

#=================================================================================================================================================#
#=================================================================================================================================================#
#=================================================================================================================================================#

def manual_display(lst_error,dicoInit,spinner):
    if dicoInit['boolQT']==False: spinner.stop()
    printcolor("\n\n  USAGE: python eKLIPse -in <FILE with Alignment paths> -ref <GBK reference> [OPTIONS]\n\n","bold yellow",dicoInit['boolColor'])
    printcolor("  OPTIONS:\n","bold white",dicoInit['boolColor'])
    printcolor("    -out          <str>  : Output directory path                  [current]\n","white",dicoInit['boolColor'])
    printcolor("    -tmp          <str>  : Temporary directory path               [/tmp]\n","white",dicoInit['boolColor'])
    printcolor("    -scsize       <int>  : Soft-clipping minimal length           [25]\n","white",dicoInit['boolColor'])
    printcolor("    -mapsize      <int>  : Upstream mapping length                [20]\n","white",dicoInit['boolColor'])
    printcolor("    -downcov      <int>  : Downsampling reads number              [500000] (0=disable)\n","white",dicoInit['boolColor'])
    printcolor("    -minq         <int>  : Read quality threshold                 [20]\n","white",dicoInit['boolColor'])
    printcolor("    -minlen       <int>  : Read length threshold                  [100]\n","white",dicoInit['boolColor'])
    printcolor("    -shift        <int>  : Breakpoint BLAST shift length          [5]\n","white",dicoInit['boolColor'])
    printcolor("    -minblast     <int>  : Minimal number of BLAST per breakpoint [1]\n","white",dicoInit['boolColor'])
    printcolor("    -bilateral    <bool> : Filter non-bilateral BLAST deletions   [True]\n","white",dicoInit['boolColor'])
    printcolor("    -mitosize     <int>  : Minimal resulting mitochondria size    [1000]\n","white",dicoInit['boolColor'])
    printcolor("    -id           <int>  : BLAST %identity threshold              [80]\n","white",dicoInit['boolColor'])
    printcolor("    -cov          <int>  : BLAST %coverage threshold              [70]\n","white",dicoInit['boolColor'])
    printcolor("    -gapopen      <int>  : BLAST cost to open a gap               [0:proton, 5:illumina]\n","white",dicoInit['boolColor'])
    printcolor("    -gapext       <int>  : BLAST cost to extend a gap             [2]\n","white",dicoInit['boolColor'])
    printcolor("    -thread       <int>  : Number of thread to use                [2]\n","white",dicoInit['boolColor'])
    printcolor("    -samtools     <str>  : samtools bin path                      [$PATH]\n","white",dicoInit['boolColor'])   
    printcolor("    -blastn       <str>  : BLASTn bin path                        [$PATH]\n","white",dicoInit['boolColor'])
    printcolor("    -makeblastdb  <str>  : makeblastdb bin path                   [$PATH]\n","white",dicoInit['boolColor'])    
    printcolor("    -circos       <str>  : circos bin path                        [$PATH]\n","white",dicoInit['boolColor'])
    printcolor("    --test               : eKLIPse test\n","white",dicoInit['boolColor'])
    printcolor("    --nocolor            : Disable output colors\n\n","white",dicoInit['boolColor'])
    printerror(lst_error)
    try: shutil.rmtree(dicoInit['pathOutputDir']) ; shutil.rmtree(dicoInit['pathTmpDir'])
    except: pass
    sys.exit()   



#=================================================================================================================================================#
#=================================================================================================================================================#
#=================================================================================================================================================#

def config_display(dicoInit):
    print "\n"
    format_time = time.strftime('%d/%m/%y %H:%M:%S',time.localtime(dicoInit['startTime']))
    printcolor("    Start time           : ","bold white",dicoInit['boolColor'])
    printcolor(format_time+"\n","white",dicoInit['boolColor']) ; time.sleep(0.1)
    printcolor("    Input file path      : ","bold white",dicoInit['boolColor'])
    printcolor(dicoInit['pathInput']+"\n","white",dicoInit['boolColor']) ; time.sleep(0.1)
    printcolor("    Reference            : ","bold white",dicoInit['boolColor'])
    printcolor(dicoInit['pathGbkRef']+"\n","white",dicoInit['boolColor']) ; time.sleep(0.1)
    printcolor("    Output folder        : ","bold white",dicoInit['boolColor'])
    printcolor(dicoInit['pathOutputDir']+"\n","white",dicoInit['boolColor']) ; time.sleep(0.1)
    printcolor("    Temporary folder     : ","bold white",dicoInit['boolColor'])
    printcolor(dicoInit['pathTmpDir']+"\n","white",dicoInit['boolColor']) ; time.sleep(0.1)
    printcolor("    Read threshold       : ","bold white",dicoInit['boolColor'])
    printcolor("minQ="+str(dicoInit['minQ'])+" | minlen="+str(dicoInit['minlen'])+"\n","white",dicoInit['boolColor']) ; time.sleep(0.1)
    printcolor("    SC threshold         : ","bold white",dicoInit['boolColor'])
    printcolor("SCsize="+str(dicoInit['SCsize'])+" | MappedPart="+str(dicoInit["MappedPart"])+"\n","white",dicoInit['boolColor']) ; time.sleep(0.1)
    printcolor("    Deletion shift       : ","bold white",dicoInit['boolColor'])
    printcolor(str(dicoInit['delShift'])+"\n","white",dicoInit['boolColor']) ; time.sleep(0.1)
    printcolor("    Min mito size        : ","bold white",dicoInit['boolColor'])
    printcolor(str(dicoInit['mitosize'])+"\n","white",dicoInit['boolColor']) ; time.sleep(0.1)
    printcolor("    Min breakpoint BLAST : ","bold white",dicoInit['boolColor'])
    printcolor(str(dicoInit['minblast'])+"\n","white",dicoInit['boolColor']) ; time.sleep(0.1)
    printcolor("    Filter non-bilateral : ","bold white",dicoInit['boolColor'])
    printcolor(str(dicoInit['bilateral'])+"\n","white",dicoInit['boolColor']) ; time.sleep(0.1)
    printcolor("    BLAST thresholds     : ","bold white",dicoInit['boolColor'])
    printcolor("id="+str(dicoInit['blastIdThreshold'])+" | cov="+str(dicoInit['blastCovThreshold'])+" | gapopen="+str(dicoInit['blastGapOpen'])+" | gapext="+str(dicoInit['blastGapExt'])+"\n","white",dicoInit['boolColor']) ; time.sleep(0.1)
    printcolor("    Downsampling         : ","bold white",dicoInit['boolColor'])
    if dicoInit['downCov']==0 : printcolor("none\n","white",dicoInit['boolColor']) ; time.sleep(0.1)
    else : printcolor(str(dicoInit['downCov'])+"\n","white",dicoInit['boolColor']) ; time.sleep(0.1)
    printcolor("    Threads number       : ","bold white",dicoInit['boolColor'])
    printcolor(str(dicoInit['nbThread'])+"\n","white",dicoInit['boolColor']) ; time.sleep(0.1)
    bamTable = tabulate_bam(dicoInit['dicoBam'],dicoInit['lstTitleBam'],).replace("\n","\n      ")
    printcolor("    Input alignments     : ","bold white",dicoInit['boolColor'])
    printcolor(str(len(dicoInit['dicoBam']))+"\n","white",dicoInit['boolColor'])
    printcolor("      "+bamTable+"\n","white",dicoInit['boolColor']) ; time.sleep(0.1)

#=================================================================================================================================================#
#=================================================================================================================================================#
#=================================================================================================================================================#

def tabulate_bam(dicoBam,lstTitleBam):
    bamContent = [] ; bamHeader = ["Name","BAM header","Nb Reads"]
    for titleBam in lstTitleBam:
        bamContent.append([titleBam,dicoBam[titleBam]['refName'],dicoBam[titleBam]['nbReads']])
    return str(tabulate.tabulate(bamContent,bamHeader,tablefmt="psql",stralign="center",numalign="center").encode('utf-8'))

#=================================================================================================================================================#
#=================================================================================================================================================#
#=================================================================================================================================================#