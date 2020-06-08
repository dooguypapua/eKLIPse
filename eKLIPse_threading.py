import sys
import tqdm
import threading
from eKLIPse_fct import *
from eKLIPse_sc import *
from eKLIPse_circos import *




def launch_threads_cli(dicoInit,description,target_fct,thread_nameformat):
    lstThread = []
    lstError = []
    lstThreadcopy = []
    PrevNbFinishThread = 0
    # One thread per BAM
    for titleBam in dicoInit["dicoBam"].keys():
        if target_fct=="Alignment_downsampling":
            if dicoInit['dicoBam'][titleBam]['nbReads']>dicoInit['downCov']: lstThread.append(threading.Thread(target=Alignment_downsampling, args=(titleBam,dicoInit,lstError), name=thread_nameformat+titleBam))
            else: dicoInit['dicoBam'][titleBam]['path_downsampling'] = dicoInit['dicoBam'][titleBam]['path']
        elif target_fct=="Read_alignment": lstThread.append(threading.Thread(target=Read_alignment, args=(titleBam,dicoInit,lstError), name=thread_nameformat+titleBam))
        elif target_fct=="SC_blast": lstThread.append(threading.Thread(target=SC_blast, args=(titleBam,dicoInit,lstError), name=thread_nameformat+titleBam))
        elif target_fct=="deletionPrediction": lstThread.append(threading.Thread(target=deletionPrediction, args=(titleBam,dicoInit,lstError), name=thread_nameformat+titleBam))
        elif target_fct=="circosConf": lstThread.append(threading.Thread(target=circosConf, args=(titleBam,dicoInit,lstError), name=thread_nameformat+titleBam))
        elif target_fct=="circosPlot": lstThread.append(threading.Thread(target=circosPlot, args=(titleBam,dicoInit,lstError), name=thread_nameformat+titleBam))
    # Tasks variables
    NbTasks = len(lstThread)
    lstThreadcopy = lstThread
    # Launch threads
    if NbTasks>0:
        printcolor("\n  "+description,'bold yellow',dicoInit['boolColor'])
        with tqdm.tqdm(total=NbTasks,ncols=80,desc="        <Running> ",position=1,leave=False) as pbar: # tqdm loop
            # Continue while some tasks are availables
            while (len(lstThread)!=0 or thread_nameformat in str(threading.enumerate())) and len(lstError)==0:
                nbThreadActive = str(threading.enumerate()).count(thread_nameformat)
                if nbThreadActive<dicoInit['nbThread']:
                    try: thread = lstThread.pop() ; thread.start() ; nbThreadActive+=1
                    except: pass
                NbFinishThread = NbTasks-len(lstThread)-nbThreadActive
                if PrevNbFinishThread!=NbFinishThread:
                    pbar.update(NbFinishThread-PrevNbFinishThread)
                    PrevNbFinishThread = NbFinishThread
            pbar.close()
        for thread in lstThreadcopy:
            try : thread.join()
            except : pass
        printcolor("\b"*57,"white",dicoInit['boolColor'])
        printerror(lstError)        
        printcolor("[OK]\n","white",dicoInit['boolColor'])


def launch_threads_qt(dicoInit,description,target_fct,thread_nameformat,qtprocesslabel):
    lstThread = []
    lstError = []
    lstThreadcopy = []
    lst_progress = []
    PrevNbFinishThread = 0
    # One thread per BAM
    for titleBam in dicoInit["dicoBam"].keys():
        if target_fct=="Alignment_downsampling":
            if dicoInit['dicoBam'][titleBam]['nbReads']>dicoInit['downCov']: lstThread.append(threading.Thread(target=Alignment_downsampling, args=(titleBam,dicoInit,lstError), name=thread_nameformat+titleBam))
            else: dicoInit['dicoBam'][titleBam]['path_downsampling'] = dicoInit['dicoBam'][titleBam]['path']
        elif target_fct=="Read_alignment": lstThread.append(threading.Thread(target=Read_alignment, args=(titleBam,dicoInit,lstError), name=thread_nameformat+titleBam))
        elif target_fct=="SC_blast": lstThread.append(threading.Thread(target=SC_blast, args=(titleBam,dicoInit,lstError), name=thread_nameformat+titleBam))
        elif target_fct=="deletionPrediction": lstThread.append(threading.Thread(target=deletionPrediction, args=(titleBam,dicoInit,lstError), name=thread_nameformat+titleBam))
        elif target_fct=="circosConf": lstThread.append(threading.Thread(target=circosConf, args=(titleBam,dicoInit,lstError), name=thread_nameformat+titleBam))
        elif target_fct=="circosPlot": lstThread.append(threading.Thread(target=circosPlot, args=(titleBam,dicoInit,lstError), name=thread_nameformat+titleBam))
    # Tasks variables
    NbTasks = len(lstThread)
    lstThreadcopy = lstThread
    # Launch threads
    if NbTasks>0:
        sys.stdout.write(description+"\n")
        sys.stdout.flush()
        # Continue while some tasks are availables
        while (len(lstThread)!=0 or thread_nameformat in str(threading.enumerate())) and len(lstError)==0:
            nbThreadActive = str(threading.enumerate()).count(thread_nameformat)
            if nbThreadActive<dicoInit['nbThread']:
                try: thread = lstThread.pop() ; thread.start() ; nbThreadActive+=1
                except: pass
            NbFinishThread = NbTasks-len(lstThread)-nbThreadActive
            if not NbFinishThread in lst_progress:
                sys.stdout.write(qtprocesslabel+str(int((float(NbFinishThread)*100.0)/float(NbTasks)))+"\n")
                sys.stdout.flush()
                lst_progress.append(NbFinishThread)
            if PrevNbFinishThread!=NbFinishThread: PrevNbFinishThread = NbFinishThread
        for thread in lstThreadcopy:
            try : thread.join()
            except : pass
        printerror(lstError)
        sys.stdout.write(qtprocesslabel+"100\n")
        sys.stdout.flush()
