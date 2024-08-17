#!/usr/bin/env python
import sys
import getopt
import os, glob
from collections import defaultdict
import re
import numpy as np
from scipy.optimize import curve_fit
from scipy.misc import factorial
import scipy.stats as ss
import scipy.optimize as so

from scipy.stats import nbinom
import matplotlib.pyplot as plt
#import statsmodels
from statsmodels.discrete.discrete_model import NegativeBinomial

#import scipy.optimize as so
#-i C:\Users\falkh\Dropbox\Tec2\SNPs\T2/FB_reGT_ali_30_25.vcf -o tee


usage_message = '''
vcf2cons.py - Falk.Hildebrand@gmail.com
rewrite of vcf2cons.pl, including more advanced error model
'''
#precompile regex
regLen = re.compile('L=(\d+)=')
regAZ = re.compile('^([A-Za-z.])(,[A-Za-z])*$')
regAF = re.compile('AF=(-?\d[\.\d+]?)')
regDP = re.compile('DP=(\d+)')
regQA = re.compile('QA=(\d+)')
regINDEL = re.compile('INDEL')

if sys.version > '3':
    long = int

def likelihood_f((n,p,loc), x, neg=-1):
    n=np.round(n) #by definition, it should be an integer
    loc=np.round(loc)
    return neg*(np.log(ss.nbinom.pmf(x, n, p, loc))).sum()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg
class Options():
    def __init__(self):
        self.host_mem = 0.9
class indexCls():
    def __init__(self,ins):
        tfa = ins.split(':')
        self.AOi = -1
        if ("AO" in tfa):
            self.AOi = tfa.index("AO")
        self.ROi = tfa.index("RO")
        self.QAi = -1
        if "QA" in tfa:
            self.QAi = tfa.index("QA")
        self.DPi = tfa.index("DP")
        #print ins,self.DPi, self.ROi
        #exit(33)

class ChromStats():
    def __init__(self, name,idx=0,leng=0):
        self.name=name
        self.idx=idx
        self.reset(leng)
        self.depthStat = defaultdict(lambda: defaultdict(lambda: 0))


    def reset(self,leng):
        # empty ints
        self.seq = ['n']*leng;  self.qual = 0;
        self.gaps = [];  self.spos = []
        self.lastPos=0
        self.bcnts = 0; self.ccnts = 0;   self.lbcnts = 0;   self.lccnts = 0

    def extendSeq(self,last_pos):
        if (last_pos > len(self.seq)):
            self.seq += "n" * (last_pos - len(self.seq) - 1)

#core routine
    def evalSNP(self,ts,t,ref,alt,altidx,freq,chPos,Cidx,_D,_d,_Q):
        #GT:DP:AD:RO:QR:AO:QA:GL
        #DPi = 1;AOi=5;ROi=3;QAi=6; #indexes for important information
        tsa=ts.split(':')
        if (len(t[3]) == 1 ): #and regINDEL.search(ts) == None
            # a SNP or reference
            q = 0
            b = ref
            dep = 0
            #regDPres = regDP.search(t[7])
            #if (regDPres != None):
            #    dep = int(regDPres.group(1))
                # print t[7] + " " + str(dep) + "\n";exit(4)
            if tsa[Cidx.DPi]is not ".":
                dep = int(tsa[Cidx.DPi])
            if (dep < 1):
                b = 'N'
            if ("," in alt or len(alt) > 1):
                print t + "\n"
                exit(3)
            altset=False
            if freq > 0:
                lfreq=0
                if Cidx.AOi >=0:
                    AO=0;RO=0;
                    tsAOs = tsa[Cidx.AOi].split(',')
                    if len(tsAOs)>altidx and tsAOs[altidx] is not ".":
                        AO = float(tsAOs[altidx])
                    if tsa[Cidx.ROi] is not ".":
                        RO = float(tsa[Cidx.ROi])
                    #print AO,RO
                    if (AO>0 or RO > 0):
                        lfreq = AO/(RO+AO)
                if (lfreq > 0.501 and alt != '.' and dep > 1):
                    #print tsa[Cidx.AOi], tsa[Cidx.ROi], lfreq,dep;
                    b = alt
                    altset=True
                    self.ccnts += 1
                    self.lccnts += 1
                    self.spos.append(t[1])
                    self.depthStat[dep]['alt'] += 1
                    self.depthStat[dep]['norm'] -= 1
                    QA=0
                    #if len(tsa) < Cidx.QAi:
                    #    print ts; exit(55)
                    if tsa[Cidx.QAi] is not ".":
                        QA = tsa[Cidx.QAi]
                    if (QA > _Q and dep >= _d and dep < _D):
                        b = b.upper()
                    else:
                        b = b.lower()
            if self.lastPos != chPos or altset:
                if self.lastPos != chPos:
                    self.depthStat[dep]['norm'] += 1
                    self.bcnts += 1
                    self.lbcnts += 1
                #print len(self.seq),chPos
                self.seq[chPos-1] = b
        elif (t[4] not in '.'):
            print t + "\n"
            exit(37)
            self.gaps.append(t[1], len(t[3]))

        self.lastPos = chPos

    def depthStats(self, oh):
        matrix = []
        dsStr = ''

        for dep in sorted(self.depthStat.keys()):
            dsStr += str(dep) + "\t"
            cntsD = 0
            if ('alt' in self.depthStat[dep]):
                dsStr += str(self.depthStat[dep]['alt']) + "\t"
                cntsD += self.depthStat[dep]['alt']
            else:
                dsStr += "0\t"
            if ('norm' in self.depthStat[dep]):
                dsStr += str(self.depthStat[dep]['norm']) + "\n"
                cntsD += self.depthStat[dep]['norm']
            else:
                dsStr += "0\n"

            matrix.append([dep, cntsD])

        oh.write(dsStr)

        matrix = np.array(matrix)
        matrix_sum = np.sum(matrix[:, 1])
        cutoffs = (int(matrix_sum * 0.025 + 0.5), int(matrix_sum * 0.975 + 0.5))
        cs = np.cumsum(matrix[:, 1])
        print cs
        print cutoffs
        lowDep = 1
        if np.sum(cs < cutoffs[0]) != 0:
            lowDep = np.max(matrix[cs < cutoffs[0], 0])
        highDep = 1
        if np.sum(cs > cutoffs[1]) != 0:
            highDep = np.min(matrix[cs > cutoffs[1], 0])

        print lowDep, highDep, cutoffs, "\n"

        # fit with curve_fit
        # parameters, cov_matrix = curve_fit(poisson, matrix[:,0], matrix[:,1])
        # print parameters, cov_matrix

        # Use negative binomial instead..
        # plt.plot(range(0, 30000), ss.nbinom.pmf(range(0, 30000), n=3, p=1.0 / 300, loc=0), 'g-')
        # bins = plt.hist(all_hits, 100, normed=True, alpha=0.8)
        n, p =  0.000105, 3.86
        #n, p = 0.4, 0.4
        rv = nbinom(n, p)
        matFrac =  matrix[:, 1] #matrix_sum
        #res = NegativeBinomial(np.ones(len(matrix[:, 1])), matrix[:, 1])
        #res.fit_regularized()
        #https://stackoverflow.com/questions/29338152/using-scipy-optimize-to-fit-data-to-negative-binomial
        xopt, fopt, iterations, funcalls, warn = so.fmin(likelihood_f, (3, 1.0 / 300, 0), args=(matrix[:,1], -1),
                                                         full_output=True, disp=False)
        print 'optimal solution: r=%f, p=%f, loc=%f' % tuple(xopt)
        print 'log likelihood = %f' % (fopt)
        print 'ran %d iterations with %d function calls' % (iterations, funcalls)
        print 'warning code %d' % (warn)

        x=matrix[:,0]
        fig, ax = plt.subplots(1, 1)
        ax.plot(x, matrix[:,1], 'bo', ms=8, label='real data', alpha=0.5)
        ax.vlines(x, 0, rv.pmf(x) * matrix_sum, colors='k', linestyles='-', lw=1, label='nbiom pmf')
        ax.legend(loc='best', frameon=False)
        #rf.cdf(x)
        plt.show()
        nbinom.cdf(x, n, p)

        print "Done depth analysis\n"


#general options for program
opt = Options()
cp = 0


def parse_opt(argv):
    try:
        opts, args = getopt.getopt(argv, "hm:o:i:n:",
                                   ["help",
                                    "out-dir=",
                                    "out-name="
                                    ])
    except getopt.error as msg:
        raise Usage('vcf2cons.py\n' + str(msg))
    if len(opts) == 0:
        raise Usage('vcf2cons.py\n' + usage_message)

    global opt
    opt.temp_dir = ""
    opt.out_name = "consens"
    for option, value in opts:
        if option in ("-h", "--help"):
            print('vcf2cons.py\n' + usage_message)
            exit(0)
        elif option in ("-o", "--out-dir"):
            if value is "-":
                opt.out_dir = value
            else:
                opt.out_dir = value + "/"
        elif option in ("-n", "--out-name"):
            opt.out_name = value
        elif option in ("-i"):
            opt.inF = value
        else:
            raise Usage("Invalid option %s", option)

    print opt.out_dir + "\n"
    if opt.temp_dir == "":
        opt.temp_dir = opt.out_dir + "tmp/"
    if not os.path.exists(opt.out_dir):
        os.makedirs(opt.out_dir)
# v2q_post_process(last_chr, last_pos, CStats[sm], 5, outhandle)
# v2q_post_process(last_chr, last_pos, seq, qual, gaps, 5, lbcnts, lccnts, spos, outhandle)
#def v2q_post_process(chr, seq, qual, gaps, l, reports, replaces, pos, outhandle):
def v2q_post_process(chr, last_pos, CS, outhandle):

    #	for g in gaps:
    #		beg=0;
    #		if  g[0] > l:
    #			beg = g[0]-l
    #		end = g[0]+g[1]+l
    #		if (end > len(seq)):
    #			end = len(seq)
    #print str(CS.idx)+"\n"
    outhandle[CS.idx-9].write( ">" + chr + " COV=" + str(CS.lbcnts) +
            " REPL=" + str(CS.lccnts) + " POS=" + ",".join(CS.spos) + "\n" + ''.join(CS.seq) + "\n")


# poisson function, parameter lamb is the fit parameter
def poisson(k, lamb):
    return (lamb ** k / factorial(k)) * np.exp(-lamb)



def main(argv=None):
    # return 3
    print sys.argv

    try:

        parse_opt(sys.argv[1:])
        CStats = {}
        # arrays
        nsmpls=1
        last_chr = ""
        last_pos = 0;  last_chr = ""
        lcnt = 0;  chromL = 0
        last_pos = 0
        _Q = 20;  _d = 1;   _D = 30000
        samples=[]

        outhandle = [sys.stdout]


        if opt.inF is "-":
            data = sys.stdin.readlines()
        else:
            filO = open(opt.inF, "r")
            data = filO.readlines()
            filO.close()
        #if opt.out_dir is not "-":


        print "easy", len(data), "lines\n"
        for line in data:
            lcnt += 1
            if line.startswith("##"):
                continue
            line=line.rstrip()
            t = line.split('\t')
            if line.startswith("#"):
                samples = t[9:] #each sample get an outfile
                nsmpls = len(samples)
                if opt.out_dir is "-":
                    if len(samples) > 1:
                        print "print to Stdout requested, but more than one sample in vcf.. aborting\n"
                        exit (543)
                else:
                    it=0
                    for sm in samples:
                        if it ==0:
                            outhandle[it] = open(opt.out_dir + opt.out_name + "."+sm + ".fna", "w")
                        else:
                            outhandle.append(open(opt.out_dir + opt.out_name + "." + sm + ".fna", "w"))
                        CStats[sm]=ChromStats(sm,it+9,0)
                        it+=1
                print "Using "+str(len(samples))+" samples\n"
                continue
            # print t[0]+t[2]+"\n"+line+"\n"
            # exit(0)
            if (last_chr == ""):#first round..
                last_chr = t[0]
                chromL = int(regLen.search(last_chr).group(1))  # ;  chromL = _; #}#die "$chromL\n";}
                for sm in samples:
                    CStats[sm].reset(chromL)

            #write out to file the complete chromosome
            if (last_chr != t[0]):
                last_chr = t[0]
                if chromL > last_pos:
                    last_pos = chromL + 1
                chromL = int(regLen.search(last_chr).group(1))
                if (last_chr > 0):
                    for sm in CStats:
                        #CStats[sm].extendSeq(last_pos)
                        #v2q_post_process(last_chr, last_pos, seq, qual, gaps, 5, lbcnts, lccnts, spos, outhandle)
                        v2q_post_process(last_chr, last_pos, CStats[sm], outhandle)
                        CStats[sm].reset(chromL)
                #reset variables
                last_pos=0


            t1 = int(t[1])
            #if (t1 - last_pos) > 1:
                #fill all seqs with gaps
                #for sm in CStats:
                    #CStats[sm].seq += 'n' * (t1 - last_pos - 1)

            if (t1 - last_pos < 0):
                #something entirely wrong with the input
                print ("[vcf2cons] unsorted input\n" + t[1] + "-" + str(last_pos) + "\non line " + str(lcnt) + "\n")
                exit(26)
            ref = t[3]
            altp = t[4]
            alta = altp.split(',')
            regAZres = regAZ.search(altp)
            if regAZres is None: #weird alt entry
                print "Alt is weird: ",altp;exit(83)
                #alt = regAZres.group(1)
            freq=0.0
            regAFres = regAF.search(t[7])
            if (regAFres != None):
                #is an alternate SNP
                freq = float(regAFres.group(1))

            Cidx = indexCls(t[8])
            altidx=0
            for alt in alta:
                #now go over each sample, inserting ref or alt allele
                for sm in CStats:
                    #print t[CStats[sm].idx];  exit(33)
                    CStats[sm].evalSNP(t[CStats[sm].idx],t,ref,alt,altidx,freq,t1,Cidx,_D,_d,_Q)
                altidx += 1
            last_pos = t1
            #end of main loop

        #just finish up the last chromosome...
        if (chromL > last_pos):
            last_pos = chromL + 1
        for sm in CStats:
            #print 'dd'
            v2q_post_process(last_chr, last_pos, CStats[sm], outhandle)
            #CStats[sm].extendSeq(last_pos)
            # last time write out derrived sequence
            #v2q_post_process(last_chr, last_pos, seq, qual, gaps, 5, lbcnts, lccnts, spos, outhandle)

        if opt.out_dir is not "-":
            for i in range(len(outhandle)):
                outhandle[i].close()

        for sm in CStats:
            outDepStat = sys.stderr
            if (opt.out_dir is not '-'):
                outDepStat=open(opt.out_dir + opt.out_name + "."+CStats[sm].name + ".depStats", "w")
            CStats[sm].depthStats(outDepStat)
            if (opt.out_dir is not '-'):
                outDepStat.close()

    except Usage as err:
        print(sys.argv[0].split("/")[-1] + ": " + str(err.msg))



    return 0

if __name__ == "__main__":
    sys.exit(main())
