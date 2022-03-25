import os
import numpy as np
import re
import subprocess
import math
import glob

# make rose input file

def make_infile(treefile, hmmfile): # infilename, outfilename):

    global convroseinfile 
    convroseinfile = 'convrosein'
    
    global convroseoutfile 
    convroseoutfile = 'convroseout'
    
    tree = open(treefile, 'r').read()

    speccount = 0
    speconly = re.split('\(|\)|,|:',tree) 
    for i in range(0, len(speconly)): 
        if speconly[i].isalpha() == True:
            speccount += 1

    os.system('hmmemit ' + hmmfile + '> hmmseq.fa')
    seqfile = open('hmmseq.fa', 'r').readlines()[1:]
    startseq = ''
    for i in range(0, len(seqfile)): 
        startseq += re.sub('\n', '', seqfile[i])

    ## MAKE ROSE INPUT FILE
    roseinfile = open(convroseinfile, 'w')

    roseinfile.write('%include Supp_Files/protein-defaults')
    roseinfile.write('\n')
    roseinfile.write('SequenceNum = ' + str(speccount))
    roseinfile.write('\n')
    roseinfile.write('ChooseFromLeaves = True')
    roseinfile.write('\n')
    roseinfile.write('TheTree = ' + tree)
    roseinfile.write('\n')
    roseinfile.write('TheSequence = ' + '"' + startseq + '"')
    roseinfile.write('\n')
    roseinfile.write('OutputFilebase = "' + convroseoutfile + '"')
    roseinfile.write('\n')

    roseinfile.close()



def run_sim(specieslist1, specieslist2):

    ## RUN
    
    roseoutfile = convroseoutfile + '.fas' # stem only in input
    
    # run ROSE

    runcmd = 'rose ' + convroseinfile

    os.system(runcmd)

    ## HOMOLOGY HYPOTHESIS

    # align output sequences 

    aligncmd = 'clustalo -i ' + roseoutfile + ' -o ' + roseoutfile + '_hom_ali --outfmt=phy --force'
    os.system(aligncmd)

    # find alignment length (for normalization)
    homalidata = np.genfromtxt(roseoutfile + '_hom_ali', dtype=str, delimiter='\n')
    homalilen = float(str.split(homalidata[0])[1])


    # make phylip input file
    phyin = open('phy_homsim_homin' ,'w')
    phyin.write(roseoutfile + '_hom_ali')
    phyin.write('\n')
    phyin.write('y')
    phyin.write('\n')
    phyin.close()

    # clear any existing PHYLIP output files
    os.system('rm outtree')
    os.system('rm outfile')

    # run phylip proml
    phycmd = 'proml < phy_homsim_homin'
    os.system(phycmd)
    
    global homsim_hompromlout 
    homsim_hompromlout = roseoutfile + '_hom_proml_outfile'

    # rename proml output files
    move1cmd = 'mv outfile ' + homsim_hompromlout
    move2cmd = 'mv outtree ' + roseoutfile + '_hom_proml_outtree'

    os.system(move1cmd)
    os.system(move2cmd)

    ## CONVERGENCE HYPOTHESIS 

    # make two files with different partitions of sequences

    # index so can use sfetch
    indexcmd = 'esl-sfetch --index ' + roseoutfile 

    os.system(indexcmd)

    # use sfetch

    extcmd1 = 'esl-sfetch -f ' + roseoutfile + ' ' + specieslist1 + ' > ' + roseoutfile + '_conv_A'
    extcmd2 = 'esl-sfetch -f ' + roseoutfile + ' ' + specieslist2 + ' > ' + roseoutfile + '_conv_B'

    os.system(extcmd1)
    os.system(extcmd2)

    # align two sequence files

    aligncmd1 = 'clustalo -i ' + roseoutfile + '_conv_A -o ' + roseoutfile + '_conv_A_ali --outfmt=phy --force'
    aligncmd2 = 'clustalo -i ' + roseoutfile + '_conv_B -o ' + roseoutfile + '_conv_B_ali --outfmt=phy --force'

    os.system(aligncmd1)
    os.system(aligncmd2)

    # find alignment lengths (for normalization)

    convaliAdata = np.genfromtxt(roseoutfile + '_conv_A_ali', dtype=str, delimiter='\n')
    convaliAlen = float(str.split(convaliAdata[0])[1])

    convaliBdata = np.genfromtxt(roseoutfile + '_conv_B_ali', dtype=str, delimiter='\n')
    convaliBlen = float(str.split(convaliBdata[0])[1])


    # proml with phylip
    # remove any existing outfiles

    os.system('rm outtree')
    os.system('rm outfile')

    # make phylip input file 1
    phyin = open('phy_convA_in' ,'w')
    phyin.write(roseoutfile + '_conv_A_ali')
    phyin.write('\n')
    phyin.write('y')
    phyin.write('\n')
    phyin.close()

    # make phylip input file 2
    phyin = open('phy_convB_in' ,'w')
    phyin.write(roseoutfile + '_conv_B_ali')
    phyin.write('\n')
    phyin.write('y')
    phyin.write('\n')
    phyin.close()

    # run proml
    
    global homsim_convA_promlout
    homsim_convA_promlout = roseoutfile + '_conv_A_proml_outfile'
    
    
    phycmd = 'proml < phy_convA_in'
    os.system(phycmd)
    movecmd = 'mv outfile ' + homsim_convA_promlout
    os.system(movecmd)
    movecmd = 'mv outtree ' + roseoutfile + '_conv_A_proml_outtree'
    os.system(movecmd)

    global homsim_convB_promlout
    homsim_convB_promlout = roseoutfile + '_conv_B_proml_outfile'
    

    phycmd = 'proml < phy_convB_in'
    os.system(phycmd)
    movecmd = 'mv outfile ' + homsim_convB_promlout
    os.system(movecmd) 
    movecmd = 'mv outtree ' + roseoutfile + '_conv_B_proml_outtree'
    os.system(movecmd)

    return homalilen, convaliAlen, convaliBlen

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
                       

def calc_likelihood_ratio(conv1len, conv2len, homlen):
    homlikelihood, homnormlikelihood = find_normalized_likelihood(homsim_hompromlout, homlen)
    conv1likelihood, conv1normlikelihood = find_normalized_likelihood(homsim_convA_promlout, conv1len)
    conv2likelihood, conv2normlikelihood = find_normalized_likelihood(homsim_convB_promlout, conv2len)
    likelihood_ratio = conv1likelihood+conv2likelihood-homlikelihood
    normlikelihood_ratio = conv1normlikelihood+conv2normlikelihood-homnormlikelihood
    return likelihood_ratio, normlikelihood_ratio

## all


def hom_simulation(numreps, tree, hmm, clade1, clade2):
    global likelihoodratios
    likelihoodratios = []
    global normlikelihoodratios
    normlikelihoodratios = []
    for i in range(0, numreps): 
        make_infile(tree, hmm)
        homalilen, convAalilen, convBalilen = run_sim(clade1, clade2)
        lr, nlr = calc_likelihood_ratio(convAalilen, convBalilen, homalilen)
        print('iteration ', i)
        likelihoodratios.append(lr)
        normlikelihoodratios.append(nlr)
        print('likelihoods:', lr, nlr)
        print('lengths: hom, conva, convb', homalilen, convAalilen, convBalilen)
    print(likelihoodratios)
    return likelihoodratios, normlikelihoodratios
