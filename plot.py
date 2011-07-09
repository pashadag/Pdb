import sys
import numpy as np
#from pylab import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math

width = 1.5 
curFile = 0
linetypes = ('k-v', 'r-s', 'g-^', 'b-o', 'c-D', 'm-*', 'k-', 'k-', 'k-', 'k-', 'k-', 'k-')
filelist = ()
desc = ()
gentype = ""
runtype = sys.argv[1]

legendloc = 4 #bottom right #4 bottom right, 7 center right, 1 upper right
if runtype == "human1" :
    filelist = [ "../082010/data/runh3.contigs", "../082010/data/runh2.contigs", "../082010/data/runh1.contigs", "../082010/data/runh0.contigs" ]
	#filelist = [ "../082010/data/runh1.contigs", "../082010/data/runh2.contigs", "../082010/data/runh3.contigs", "../082010/data/runh0.contigs" ]
    desc = [ "IS=5000", "IS=1000", "IS=200", "non-paired"]
    gentype = "chr22"
elif runtype == "ecoli1" :
    filelist = [ "../082010/data/run14.contigs", "../082010/data/run1.contigs", "../082010/data/run0.contigs" ] 
    desc = [ "IS=2000", "IS=200", "non-paired" ]
    gentype = "ecoli"
elif runtype == "ecoli2" :
    filelist = [ "data/runf4d.contigs", "data/runf4b.contigs", "data/runf4c.contigs" ]
    desc = [ "len=1000\n(non-paired)", "len=20", "len=10"]
    gentype = "ecoli"
    legendloc = 7 #centerright
elif runtype == "human2" :
    gentype = "chr22"
    legendloc = 7 #centerright
    filelist = [ "data/runf4h+.contigs", "data/runf4l.contigs", "data/runf4i.contigs", "data/runf4e.contigs", "data/runf4f.contigs" ]
    desc = [ "len=1000\n(non-paired)", "len=300", "len=100", "len=50", "len=20"]
    #data/chr22filt.1000.contigs.lengths" "data/runf4g.contigs"(l=10)
    #filelist = [ "data/chr22filt.1000.contigs.lengths", "data/runf4i.contigs", "data/runf4e.contigs", "data/runf4f.contigs" ]
    #desc = [ "len=1000\n(non-paired)", "len=100", "len=50", "len=20"]
    #filelist = [ "data/chr22filt.1000.contigs.lengths", "data/runn4a.contigs", "data/runn4d.contigs", "data/runn4c.contigs", "data/runn4b.contigs" ]
    #desc = [ "l=1000", "IS=5k", "IS=4k", "IS=3k", "IS=2k" ]
elif runtype == "ecoli3" :
    filelist = [ "data/runf4a.contigs", "data/runf2e.contigs", "data/runf2f.contigs", "data/runf2h.contigs", "../082010/data/run0.contigs" ]
    desc = [ r'$\Delta=0$', r'$\Delta=20$', r'$\Delta=40$', r'$\Delta=200$', "non-\npaired" ]
    gentype = "ecoli"
elif runtype == "human3" :
    filelist = [ "../082010/data/runh2.contigs", "data/runf2a.contigs", "data/runf2b.contigs", "data/runf2c.contigs", "../082010/data/runh0.contigs" ]
    desc = [ r'$\Delta=0$', r'$\Delta=5$', r'$\Delta=20$', r'$\Delta=40$', "non-\npaired" ]
    gentype = "chr22"
elif runtype == "eulerecoliold1" :
    filelist = [ "../082010/data/run14.contigs", "../082010/data/run1.contigs" ]
    desc = [ "IS=2000", "IS=200" ]
    gentype = "ecoli"
elif runtype == "eulerecoliold2" :
    filelist = [ "data/runf4b.contigs", "data/runf4c.contigs" ]
    desc = [ "len=20", "len=10"]
    gentype = "ecoli"
    legendloc = 7 #centerright
elif runtype == "eulerecoliold3" :
    filelist = [ "data/runf4a.contigs", "data/runf2e.contigs", "data/runf2f.contigs", "data/runf2h.contigs" ]
    desc = [ r'$\Delta=0$', r'$\Delta=20$', r'$\Delta=40$', r'$\Delta=200$'  ]
    gentype = "ecoli"
elif runtype == "eulerecoli1" :
    filelist = [ "data/euler/runeco1a.lengths", "data/euler/runeco1b.lengths" ]
    desc = [ 'IS=2000', "IS=200" ]
    gentype = "ecoli"
elif runtype == "eulerecoli2" :
    filelist = [ "data/euler/runeco2a.lengths", "data/euler/runeco2b.lengths" ]
    desc = [ 'len=20', "len=10" ]
    #linetypes = ('r-s', 'g-^', 'b-o', 'c-D', 'm-*', 'k-', 'k-', 'k-', 'k-', 'k-', 'k-')
    gentype = "ecoli"
elif runtype == "eulerecoli3" :
    filelist = [ "data/euler/runeco1z.lengths", "data/euler/runeco3b.lengths", "data/euler/runeco3c.lengths", "data/euler/runeco3d.lengths" ]
    desc = [ r'$\Delta=0$', r'$\Delta=20$', r'$\Delta=40$', r'$\Delta=200$' ]
    gentype = "ecoli"
elif runtype == "eulerhuman1" :
    filelist = [ "data/euler/runhum1a.lengths", "data/euler/runhum1b.lengths", "data/euler/runhum1c.lengths" ]
    desc = [ "IS=5000", "IS=1000", "IS=200"]
    gentype = "chr22"
elif runtype == "eulerhuman3" :
    filelist = [ "data/euler/runhum1b.lengths", "data/euler/runhum3b.lengths", "data/euler/runhum3c.lengths", "data/euler/runhum3d.lengths" ]
    desc = [ r'$\Delta=0$', r'$\Delta=5$', r'$\Delta=20$', r'$\Delta=40$' ]
    gentype = "chr22"
else :
    gentype = runtype
    filelist = sys.argv[2:]
    #desc = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'e', 'f', 'g', 'h', 'i')  
    desc = sys.argv[2:]

#ecolimain
#filelist = [ "data/runf3f.contigs", "data/runf3e.contigs", "data/runf3d.contigs" ]
#desc = [ "IS=5000", "IS=1000", "IS=200"]
#filelist = [ "data/runf3i.contigs", "data/runf3h.contigs", "data/runf3g.contigs", "data/rune100.contigs" ]
#desc = [ "IS=5000", "IS=1000", "IS=200", "non-paired"]


#left .2 right .8
plt.subplots_adjust(left=.15, bottom=.2, right=0.85, top=.9, wspace=None, hspace=None)


for fileName in filelist:
    print fileName
    f = open(fileName, "r")
    x = []
    y  = []
    s = []
    cumula = 0
    count = 0
    lbl = ""
    linumber = 0;
   
    for line in f:
        if linumber < 1: #number of lines in header
            lbl = line
            linumber = linumber +1
            continue
        

        x.append(int(line))
        y.append(count)
        count +=1
    
    x.sort(reverse = True)
    for i in x:
        cumula += i
        s.append(cumula)

	#width += 0.31 
    plt.plot(y[1:100],s[1:100],linetypes[curFile], markersize=10, markevery=10 + 2*curFile, label= desc[curFile], linewidth = width )
	#plt.plot(y[1:100],s[1:100],label= lbl, linewidth = width)
	#plt.plot(y[1:100],s[1:100],label= fileName, linewidth = width)
	#plt.plot(y,s,label= fileName, linewidth = width)
    
    curFile = curFile + 1;


# matplotlib.text.Text instances
#plt.axhline(0, color='black', lw=4)
#plt.title('Effect of Insert Length')

if runtype == "ecoli3" or runtype == "human3" or runtype == "eulerecoli3" or runtype == "eulerhuman3" :
    plt.xlabel('Contigs', fontsize=30 )

tickfontsize = 20
if gentype == "chr22" :
    plt.yticks((0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000, 35000000), ("0", "5M", "10M", "15M", "20M", "25M", "30M", "35M"), size=tickfontsize)
elif gentype == "ecoli" :
    plt.yticks((0, 1000000, 2000000, 3000000, 4000000, 5000000), ("0", "1M", "2M", "3M", "4M", "5M"), size=tickfontsize)
    plt.ylabel('Cumulative Length', fontsize=30)
else :
    plt.yticks(size=tickfontsize)

plt.xticks(size=tickfontsize)

if runtype == "human2" :
    leg = plt.legend(loc=legendloc, bbox_to_anchor = (0.95, 0.6)) #4 bottom right, 7 center right, 1 upper right
else :
    leg = plt.legend(loc=legendloc) #4 bottom right, 7 center right, 1 upper right

for t in leg.get_texts():
	t.set_fontsize(tickfontsize)    # the legend text fontsize
#ltext = leg.get_texts() # all the text.Text instance in the legend
#set(ltext, fontsize=35) # the legend text fontsize

#plt.figure(figsize=(17, 8))
plt.savefig(runtype + ".pdf")
#plot.show()
