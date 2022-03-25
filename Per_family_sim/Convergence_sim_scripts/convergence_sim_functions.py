import numpy as np
import glob
import os
import re
import math
import subprocess

#makes rose input files

def make_infiles(treefileA, treefileB, hmmfile):

    ## IN and OUT file names are hard-coded here
    
    global infileAname 
    infileAname = 'convroseinA'
    global infileBname
    infileBname = 'convroseinB'
    global outfileAname
    outfileAname = 'convroseoutA'
    global outfileBname 
    outfileBname = 'convroseoutB'
    
    ## MAKE INPUT FILES
    treeA = open(treefileA, 'r').read()
    treeB = open(treefileB, 'r').read()

    speccountA = 0
    speconly = re.split('\(|\)|,|:',treeA) 
    for i in range(0, len(speconly)): 
        if speconly[i].isalpha() == True:
            speccountA += 1
            
    speccountB = 0
    speconly = re.split('\(|\)|,|:',treeB) 
    for i in range(0, len(speconly)): 
        if speconly[i].isalpha() == True:
            speccountB += 1
                        
    ## MAKE ROSE INPUT FILE FOR SPECIES GROUP A

    os.system('hmmemit ' + hmmfile + '> hmmseq.fa')
    seqfile = open('hmmseq.fa', 'r').readlines()[1:]
    startseq = ''
    for i in range(0, len(seqfile)): 
        startseq += re.sub('\n', '', seqfile[i])

    roseinfile = open(infileAname, 'w')

    roseinfile.write('%include Supp_Files/protein-defaults')
    roseinfile.write('\n')
    roseinfile.write('SequenceNum = ' + str(speccountA))
    roseinfile.write('\n')
    roseinfile.write('ChooseFromLeaves = True')
    roseinfile.write('\n')
    roseinfile.write('TheTree = ' + treeA)
    roseinfile.write('\n')
    roseinfile.write('TheSequence = ' + '"' + startseq + '"')
    roseinfile.write('\n')
    roseinfile.write('OutputFilebase = "' + outfileAname + '"')
    roseinfile.write('\n')

    roseinfile.close()

    ## MAKE ROSE INPUT FILE FOR SPECIES GROUP B
    ## SAME HMM, DIFFERENT SEQUENCE

    os.system('hmmemit ' + hmmfile + '> hmmseq.fa')
    seqfile = open('hmmseq.fa', 'r').readlines()[1:]
    startseq = ''
    for i in range(0, len(seqfile)): 
        startseq += re.sub('\n', '', seqfile[i])

    roseinfile = open(infileBname, 'w')

    roseinfile.write('%include Supp_Files/protein-defaults')
    roseinfile.write('\n')
    roseinfile.write('SequenceNum = ' + str(speccountB))
    roseinfile.write('\n')
    roseinfile.write('ChooseFromLeaves = True')
    roseinfile.write('\n')
    roseinfile.write('TheTree = ' + treeB)
    roseinfile.write('\n')
    roseinfile.write('TheSequence = ' + '"' + startseq + '"')
    roseinfile.write('\n')
    roseinfile.write('OutputFilebase = "' + outfileBname + '"')
    roseinfile.write('\n')

    roseinfile.close()


## run ROSE and PROML

def run_sim(roseinA, roseinB, roseoutA, roseoutB):

    ## RUN
    
    roseoutseqsA = roseoutA + '.fas' # stem only in input
    roseoutseqsB = roseoutB + '.fas' # stem only in input

    # run ROSE, A and B

    runcmd = 'rose ' + roseinA
    os.system(runcmd)
    
    runcmd = 'rose ' + roseinB
    os.system(runcmd)

    ## HOMOLOGY HYPOTHESIS

    # combine sequence files (Separate bc convergent simulation)

    # remove existing one from previous run

    os.system('rm convsim_homseqs.fas')

    combinecmd = 'cat ' + roseoutseqsA + ' >> convsim_homseqs.fas'
    os.system(combinecmd)
    combinecmd = 'cat ' + roseoutseqsB + ' >> convsim_homseqs.fas'
    os.system(combinecmd)    

    # align output sequences 

    aligncmd = 'clustalo -i ' + 'convsim_homseqs.fas' + ' -o ' + 'convsim_homseqs' + '_ali --outfmt=phy --force'
    os.system(aligncmd)

    # find alignment length (for normalization)
    homalidata = np.genfromtxt('convsim_homseqs' + '_ali', dtype=str, delimiter='\n')
    homalilen = float(str.split(homalidata[0])[1])


    # make phylip input file
    phyin = open('phy_convsim_homseqs_in' ,'w')
    phyin.write('convsim_homseqs' + '_ali')
    phyin.write('\n')
    phyin.write('y')
    phyin.write('\n')
    phyin.close()

    # clear any existing PHYLIP output files
    os.system('rm outtree')
    os.system('rm outfile')

    # run phylip proml
    phycmd = 'proml < phy_convsim_homseqs_in'
    os.system(phycmd)

    global hompromlout
    hompromlout = 'convsim_homseqs' + '_ali' + '_proml_outfile'
    
    # rename proml output files
    move1cmd = 'mv outfile ' + hompromlout
    move2cmd = 'mv outtree ' + 'convsim_homseqs' + '_ali' + '_proml_outtree'

    os.system(move1cmd)
    os.system(move2cmd)

    ## CONVERGENCE HYPOTHESIS 

    # use existing two files with different partitions of sequences

    # align two sequence files

    aligncmd1 = 'clustalo -i ' + roseoutseqsA + ' -o ' + roseoutseqsA + '_ali --outfmt=phy --force'
    aligncmd2 = 'clustalo -i ' + roseoutseqsB + ' -o ' + roseoutseqsB + '_ali --outfmt=phy --force'

    os.system(aligncmd1)
    os.system(aligncmd2)

    # find alignment lengths (for normalization)

    convaliAdata = np.genfromtxt(roseoutseqsA + '_ali', dtype=str, delimiter='\n')
    convaliAlen = float(str.split(convaliAdata[0])[1])

    convaliBdata = np.genfromtxt(roseoutseqsB + '_ali', dtype=str, delimiter='\n')
    convaliBlen = float(str.split(convaliBdata[0])[1])

    # proml with phylip
    # remove any existing outfiles

    os.system('rm outtree')
    os.system('rm outfile')

    # make phylip input file 1
    phyin = open('phy_convsim_convseqsA_in' ,'w')
    phyin.write(roseoutseqsA + '_ali')
    phyin.write('\n')
    phyin.write('y')
    phyin.write('\n')
    phyin.close()

    # make phylip input file 2
    phyin = open('phy_convsim_convseqsB_in' ,'w')
    phyin.write(roseoutseqsB + '_ali')
    phyin.write('\n')
    phyin.write('y')
    phyin.write('\n')
    phyin.close()
    
    global convApromlout
    convApromlout = roseoutseqsA + '_ali' + '_proml_outfile'
    global convBpromlout
    convBpromlout = roseoutseqsB + '_ali' + '_proml_outfile'
    
    # run proml
    phycmd = 'proml < phy_convsim_convseqsA_in'
    os.system(phycmd)
    movecmd = 'mv outfile ' + convApromlout
    os.system(movecmd)
    movecmd = 'mv outtree ' + roseoutseqsA + '_ali' + '_proml_outtree'
    os.system(movecmd)


    phycmd = 'proml < phy_convsim_convseqsB_in'
    os.system(phycmd)
    movecmd = 'mv outfile ' + convBpromlout
    os.system(movecmd) 
    movecmd = 'mv outtree ' + roseoutseqsB + '_ali' + '_proml_outtree'
    os.system(movecmd)

    return homalilen, convaliAlen, convaliBlen



## returns log likelihood ratio of convergence over homology

def find_likelihood(proml_outfile):
    outfile = open(proml_outfile, 'r').readlines()
    for i in range(0, len(outfile)): 
        if 'Ln Likelihood' in outfile[i]: 
            likelihood = float(str.split(outfile[i])[3]) 
            return likelihood

        
def find_normalized_likelihood(proml_outfile, alilen):
    outfile = open(proml_outfile, 'r').readlines()
    for i in range(0, len(outfile)): 
        if 'Ln Likelihood' in outfile[i]: 
            likelihood = float(str.split(outfile[i])[3]) 
            normlikelihood = likelihood/alilen
            return likelihood, normlikelihood
                       

def calc_likelihood_ratio(conv1_proml_outfile, conv2_proml_outfile, hom_proml_outfile, conv1len, conv2len, homlen):
    homlikelihood, homnormlikelihood = find_normalized_likelihood(hom_proml_outfile, homlen)
    conv1likelihood, conv1normlikelihood = find_normalized_likelihood(conv1_proml_outfile, conv1len)
    conv2likelihood, conv2normlikelihood = find_normalized_likelihood(conv2_proml_outfile, conv2len)
    likelihood_ratio = conv1likelihood+conv2likelihood-homlikelihood
    normlikelihood_ratio = conv1normlikelihood+conv2normlikelihood-homnormlikelihood
    return likelihood_ratio, normlikelihood_ratio

## everything

def conv_simulation(numreps, treeA, treeB, hmm):
    global likelihoodratios
    likelihoodratios = []
    global normlikelihoodratios
    normlikelihoodratios = []
    for i in range(0, numreps): 
        make_infiles(treeA, treeB, hmm)
        homalilen, convAalilen, convBalilen = run_sim(infileAname, infileBname, outfileAname, outfileBname)
        lr, nlr = calc_likelihood_ratio(convApromlout, convBpromlout, hompromlout, convAalilen, convBalilen, homalilen)
        print('iteration ', i)
        likelihoodratios.append(lr)
        normlikelihoodratios.append(nlr)
        print('likelihoods:', lr, nlr)
        print('lengths: hom, conva, convb', homalilen, convAalilen, convBalilen)
    print(likelihoodratios)
    return likelihoodratios, normlikelihoodratios