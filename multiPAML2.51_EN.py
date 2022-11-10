# -*- coding: utf-8 -*-


##. multiPAML2.5 is useful to run multiple PAML analyses given the MSA (.fas) and the corresponding trees (.tre) tagged in the clade to test


##. REQUIREMENTS: a MSA file in fasta format for each gene, one or multiple tree files (with .tre extension) for each gene, 'branch-site.ctl' and
##.                     'branch-site_Ho.ctl' control files and this multiPAML2.5 python script file must be located in the same folder


##. EXPECTED OUTPUT: multiPAML2.5 will create the directories structure and manage necesary file to run positive selection tests for all given genes 
##.                     testing one or many branches/clades for each gene, depending on the tagged trees given for each gene. multiPAML2.5 will crate
##.                     a directory for each gene and inside it a sub-directory for each branch/clade to be tested for that gene.
##.                     Output files for each run will be parsed and most relevant results (such as log likelihood values, LRT estimator value, resulting
##.                     p-value and the table of estimated parametres) for each test will be summarized in a general 'results.txt' output file.
##.                     Positive results will be additionaly gathered in another ouput file for easy access ('PosSel-tests.txt')


##.  WARNINGS: this version of multiPAML2 (v2.5) is prepared to run positive selection. If you previously want to run multiple M0 (one-ratio test)
##.                in order to optimice branch lengthes under GY94, please use multiPAML2.5-M0 script before
##.
##.                - this script asssumens codeml binary is included in the $PATH and can be invoked from any location
##.
##.                - MSA files must bear the exact same extension as the one set in the filesExt variable (default: .fas) 
##.
##.                - in order to be properly assigned to the corresponding gene, tree files must START with the exact same name of the MSA file. After
##.                     starting with their corresponding MSA file name, tree files should vary in order to show differences between them (e.g., diffent
##.                     tagged branches or clades)
##.
##.                - sub-directories of the same gene directory, corresponding to different tested branches/clades for that gene, will be named according
##.                     to the corresponding tree file name. Specifically they will be named with MSA file name and the last part or the corresponding
##.                     tree file (everything that appears in the tree file name after the LAST underscore).
##.
##.                     Recommendation: in order to easyly distinguish among results add the name of the tagged branch/clade at the END of the tree file
##.                     name AFTER an UNDERSCORE.


##. Created by Dr. Francisco Pisciottano (fpisciottano@gmail.com) at Biology and Experimental Medicine Institute -
##. National Scientific and Technical Research Council (IBYME-CONICET), Ciudad Autónoma de Buenos Aires, Argentina

## Welcome


print "\nHi! Are we going to run some PAML analyses?\n\n  Alright! Let's start! \n  Importing modules and loading tools...\n"


## modules and tools
import os
import shutil
from scipy import stats


def openFile(filetoopen):
        """takes a file name and returns its content in a string"""
        f = open(filetoopen)
        out = f.read()
        f.close()
        return out

        
def save(nombre, content):
        """saves content to a named file
        content must be a string"""
        g = open(nombre, 'w')
        g.write(content)
        g.close

    
def listFiles(ext):
    """lists the files of a given extension in the current directory"""
    filesList = list()
    for i in os.listdir(os.getcwd()):
        if i.endswith(ext):
            filesList.append(i)
    return filesList


def listTrees(geneName):
    """returns a list of .tre files in the current directory that correspond to the given gene"""
    treesList = list()
    for i in os.listdir(os.getcwd()):
        if i.endswith('.tre') and i.startswith(geneName):
            treesList.append(i)
    return treesList


def setPAMLtest(ctlfile, seqfile, treefile, geneName):
    ctlCont = openFile(ctlfile)
    ctlCont = ctlCont.replace("seqfile = .phy", "seqfile = " + seqfile)
    ctlCont = ctlCont.replace("treefile = .tre", "treefile = " + treefile)
    ending = str()
    if ctlfile == "M0.ctl":
        ending = "-M0"
    elif ctlfile == "branch-site.ctl":
        ending = "-BS"
    elif ctlfile == "branch-site_Ho.ctl":
        ending = "-Ho"
    ctlCont = ctlCont.replace("outfile = mlc", "outfile = mlc_" + geneName + ending)
    save(ctlfile, ctlCont)


stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)



print "  Modulues and tools correctly imported and loaded. We can begin! \n"


## ########

filesExt = '.fas'  #  .fas is the default extension for MSA files

# control files names
M0ctlfile = "M0.ctl"
BSctlfile = "branch-site.ctl"
Hoctlfile = "branch-site_Ho.ctl"


# creo una lista con los nombres de todos los archivos de la extensión dada
filesList = list()
filesList = listFiles(filesExt)

if len(filesList) == 0:
    print "ERROR! No MSA files found with the given extension in current location"
    raise Exception ("No MSA files found with the given extension")
elif len(filesList) > 0:
    print "\n > %i detected MSA files. Starting to analyse them...\n" % len(filesList)


# defino el directorio desde el cual voy a correr
main = os.getcwd()



###########

alpha = 0.05

outfile = "results.txt"


# creo el output para ir metiendo los resultados de likelihoods
likelihoodOutputs = str()
seleccionados = str()

#conciseOutput = "Gene Symbol" + "," + "p value" + "," + "q value (FDR corrected)" + "," + "bool (p < %s)" % str(alpha)
conciseOutput = "Gene Symbol" + "," + "p-value" + "\n"

# creo el errorLog
errorlog = 'Detected errors: \n'


divisor = '-'*60 + '\n'
genesAnalizados = list()
pvals = list()
qvals = list()
PSgenes = list()



for seqfile in filesList:
    # inicializo/reseteo las variables de nombre e índice de árboles para este gen
    name = str()
    t = 0
    
    # obtengo el nombre del gen a partir del archivo del alineamiento
    name = seqfile[:seqfile.find(filesExt)]
    print "-- gene %s --" % name


    # busco si hay árboles que correspondan a ese gen
    if len(listTrees(name)) == 0:
        print "ERROR! No tree files detected for gene %s. Please check tree file names" % name
        raise Exception ("No tree files detected for gene %s" % name)

    else:
        # creo la carpeta para ese gen y guardo su ruta en carpetaGen
        os.mkdir(name)
        carpetaGen = "./" + name

        # empiezo a iterar por los árboles para ese gen haciendo un test con cada uno
        while t < len(listTrees(name)):
            inputTree_name = listTrees(name)[t]
            rama = inputTree_name[inputTree_name.rfind('_')+1:inputTree_name.find('.')]

            # inicializo las variables de salida para este test
            mlcBS = None
            lnL_line = str()
            lnL = str()
            mlcHo = None
            lnL0_line = str()
            lnL0 = str()
            LRT = float()
            p = float()
            estimatesTable = str()
            
            # creo las subcarpetas del análisis de esa rama/subárbol, sus subcarpetas para las hipótesis Halt y Ho
            # y guardo las rutas de estas carpetas
            carpetaRama = carpetaGen + "/" + name + "_%s" %rama
            carpetaBS = carpetaRama + "/" + name + "-BS"
            carpetaHo = carpetaRama + "/" + name + "-Ho"
            os.mkdir(carpetaRama)
            os.mkdir(carpetaBS)
            os.mkdir(carpetaHo)

            # copio los archivos necesario desde la carpeta main a las subcarpetas
            shutil.copy(seqfile, carpetaBS)
            shutil.copy(seqfile, carpetaHo)

            shutil.copy(BSctlfile, carpetaBS)
            shutil.copy(Hoctlfile, carpetaHo)

            shutil.copy(inputTree_name, carpetaBS)
            shutil.copy(inputTree_name, carpetaHo)


            # me muevo a la carpeta BS, seteo y corro el analisis branch-site para este gen
            os.chdir(carpetaBS)
            setPAMLtest(BSctlfile, seqfile, inputTree_name, name)
            print "running positive selection hypothesis for %s gene for branch %s" % (name, rama)
            os.system("codeml branch-site.ctl")
            
            # parseo el archivo de salida 'mlc-BS' y extraigo el valor de likelihood 'lnL'
            mlcBS = openFile("mlc_" + name + "-BS") # (!!) linea delicada!
            lnL_line = mlcBS[mlcBS.find("lnL"):mlcBS.find("+0.000000")]
            lnL = lnL_line[lnL_line.find("-"):].strip()
            likelihoodOutputs += name + "_%s" %rama + '\n\n' + "lnL = " + lnL + '\n'

            # extraigo también del archivo mlc-BS la tabla de parámentros estimados (proporciones y omegas de las clases de sitios)
            estimatesTable = mlcBS[mlcBS.find("site class             0        1       2a       2b"):mlcBS.find("Empirical Bayes ")-7]

            
            # vuelvo al directorio main
            os.chdir(main)


            # me muevo a la carpeta Ho, seteo y corro la hipótesis nula para este gen
            os.chdir(carpetaHo)
            setPAMLtest(Hoctlfile, seqfile, inputTree_name, name)
            print "running null hypothesis for %s gene for branch %s \n" % (name, rama)
            os.system("codeml branch-site_Ho.ctl")
            
            # parseo el archivo de salida 'mlc-Ho' y extraigo el valor de likelihood 'lnL0'
            mlcHo = openFile("mlc_" + name + "-Ho") # (!!) linea complicada!
            lnL0_line = mlcHo[mlcHo.find("lnL"):mlcHo.find("+0.000000")]
            lnL0 = lnL0_line[lnL0_line.find("-"):].strip()
            likelihoodOutputs += "lnL0 = " + lnL0 + '\n'


            #calculo el LRT y el p-valor de este test
            LRT = 2*(float(lnL)-float(lnL0))
            p = stats.chisqprob(LRT,1)
            pvals.append(p)
            likelihoodOutputs += "LRT = " + str(LRT) + '\n' + "p = " + str(p) + '\n\n' + estimatesTable + '\n\n' + divisor


            # asiento el gen analizado y si resultó seleccionado y vuelvo al directorio main
            genesAnalizados.append(name + "-" + rama)

            if p < 0.05:
                PSgenes.append(name + "-" + rama)
                
            os.chdir(main)

            
            # salvo los resultados en cada vuelta (para no perder el avance si se interrumpe)
            save("errorLog", errorlog)
            save(outfile, likelihoodOutputs)
            
            # sumo +1 al índice de árboles de este gen para pasar al siguiente
            t = t+1


###corrección x testeo múltiple
##FDRout = multipletests(pvals,alpha=alpha,method="fdr_bh")
##qvals = FDRout[1]
##posSel_bool = FDRout[0]

##for n in range(len(genesAnalizados)):
##    conciseOutput += genesAnalizados[n] + "," + str(pvals[n]) + "," + str(qvals[n]) + "," + str(posSel_bool[n])
##    if posSel_bool[n] == True:
##        PSgenes.append(genesAnalizados[n])        
##    conciseOutput += '\n'

    

seleccionados = '\n'.join(PSgenes)


print "\nReady! Analyses finished!\n"
print "Analysed genes: " + str(len(genesAnalizados))
likelihoodOutputs += "\nAnalysed genes: " + str(len(genesAnalizados)) + '\n'

print "Positive Selected genes (p < %s) : " % str(alpha) + str(len(PSgenes))
likelihoodOutputs += "Positive Selected genes (p < %s) : " % str(alpha) + str(len(PSgenes)) + '\n'
save("PosSel-tests.txt", seleccionados)

save(outfile, likelihoodOutputs)











