# -*- coding: utf-8 -*-
import os
import sys
import shutil
import time
from spinner import *
from eKLIPse_init import *
from eKLIPse_fct import *
from eKLIPse_threading import *


#***** HEAD DISPLAY *****#
boolColor = True
boolQT = False
lstError = []
if "--nocolor" in sys.argv: boolColor = False
if "--qtgui" in sys.argv: boolQT = True ; boolColor = False
if boolQT==False:
    os.system('cls' if os.name == 'nt' else 'clear')
    printcolor("\n   ___              __   __   ___ \n  |__  |__/ |    | |__) /__` |__  \n  |___ |  \ |___ | |    .__/ |___ \n\n","bold yellow",boolColor)


#*********************************************************************
#---------------------------------------------------------------------
#  [ 0 ] -  INITIALIZATION
#---------------------------------------------------------------------
#*********************************************************************
if "-in" in sys.argv and "-ref" in sys.argv: printcolor("\n  Initialization ","bold yellow",boolColor)
# Create Initialization Hastable
spinner = Spinner(0.4)
if boolQT==False: spinner.start()
dicoInit = arg_manager(sys.argv,os.path.dirname(os.path.realpath(__file__)),boolColor,boolQT,spinner)
if boolQT==False: spinner.stop()
# Display Configuration
config_display(dicoInit)



#*********************************************************************
#---------------------------------------------------------------------
#  [ 1 ] -  DOWNSAMPLING
#---------------------------------------------------------------------
#*********************************************************************
if dicoInit["downCov"]!=0:
    if dicoInit["boolQT"]==False: launch_threads_cli(dicoInit,"BAM Downsampling    ","Alignment_downsampling","DownsamplingThread_")
    else: launch_threads_qt(dicoInit,"BAM Downsampling    ","Alignment_downsampling","DownsamplingThread_","#DC")



#*********************************************************************
#---------------------------------------------------------------------
#  [ 2 ] -  COMPUTING
#---------------------------------------------------------------------
#*********************************************************************
#***** Read alignments *****#
if dicoInit["boolQT"]==False: launch_threads_cli(dicoInit,"Read alignments  ","Read_alignment","ReadThread_")
else: launch_threads_qt(dicoInit,"Read alignments  ","Read_alignment","ReadThread_","#RA")

#***** Mean coverage *****#
printcolor("\n  Compute mean depth   ","bold yellow",boolColor)
mean_coverage(dicoInit,lstError)
printerror(lstError)
printcolor("[OK]\n","white",boolColor)


#***** BLAST SC sequence *****#
if dicoInit["boolQT"]==False: launch_threads_cli(dicoInit,"Blast SC sequences ","SC_blast","ThreadSCblast_")
else: launch_threads_qt(dicoInit,"Blast SC sequences ","SC_blast","ThreadSCblast_","#BL")



#*********************************************************************
#---------------------------------------------------------------------
#  [ 3 ] -  DELETION PREDICTION & EXCEL
#---------------------------------------------------------------------
#*********************************************************************
#***** Search deletion *****#
if dicoInit["boolQT"]==False: launch_threads_cli(dicoInit,"Search deletion  ","deletionPrediction","ThreadDelPred_")
else: launch_threads_qt(dicoInit,"Search deletion  ","deletionPrediction","ThreadDelPred_","#DE")

#***** Make .csv *****#
printcolor("\n  Make results table   ","bold yellow",boolColor)
create_results_table(dicoInit,lstError)
printerror(lstError)
printcolor("[OK]\n","white",boolColor)



#*********************************************************************
#---------------------------------------------------------------------
#  [ 4 ] -  CREATE CIRCOS PLOTS
#---------------------------------------------------------------------
#*********************************************************************
#***** Make the CIRCOS configuration file *****#
if dicoInit["boolQT"]==False: launch_threads_cli(dicoInit,"Make circos .conf   ","circosConf","ThreadCircosConf_")
else: launch_threads_qt(dicoInit,"Make circos .conf   ","circosConf","ThreadCircosConf_","#CO")

#***** Launch CIRCOS *****#
if dicoInit["boolQT"]==False: launch_threads_cli(dicoInit,"Create circos plots ","circosPlot","ThreadCircosPlot_")
else: launch_threads_qt(dicoInit,"Create circos plots ","circosPlot","ThreadCircosPlot_","#PL")



#*********************************************************************
#---------------------------------------------------------------------
#  [ * ] -  CLEANING & DISPLAY
#---------------------------------------------------------------------
#*********************************************************************
# shutil.rmtree(dicoInit['pathTmpDir'])
printcolor("\n  eKLIPse Completed\n",'bold green',boolColor)
printcolor("    End time           : ","bold green",boolColor)
printcolor(time.strftime('%d/%m/%y %H:%M:%S',time.localtime(time.time()))+"\n","green",boolColor)
(t_min, t_sec) = divmod(time.time()-dicoInit['startTime'],60)
(t_hour,t_min) = divmod(t_min,60) 
printcolor("    Run time           : ","bold green",boolColor)
printcolor("{}h:{}m:{}s".format(int(t_hour),int(t_min),int(t_sec))+"\n\n","green",boolColor)