# -*- coding: utf-8 -*-
"""
Created on Fri May 14 15:46:06 2021

@author: Alexia
"""
from scipy import stats
from fitter import Fitter
import pandas as pd
import numpy as np
import os
from os import listdir
from os.path import isfile, join
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import pylab as plt
import pylab as pzt
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pylab as plt
import bioinfokit
from bioinfokit.visuz import cluster
import matplotlib.pyplot as plt
import pylab as plt
from fpdf import FPDF
arr = os.listdir()
####  all scripts functions
def create_pdf(val1,val2,val3,val4):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size = 25)
    
    # create a cell
    pdf.cell(200, 15, txt = "Plotting Report",
    		ln = 1, align = 'C')
    # add another cell
    pdf.set_font("Arial", size = 20)
    pdf.cell(200, 15, txt = "The script has produced the following plots: ",
    		ln = 2, align = 'C')
    pdf.set_font("Arial", size = 15)
    # add another cell
    pdf.cell(200, 10, txt = "1) PCA plot",
    		ln = 4, align = 'C')
    # add another cell
    m='Total explained variance :'+val1+' %'
    pdf.cell(200, 10, txt = m,
    		ln = 5, align = 'C')
    n='Variance of each component respectively:'+val2
    pdf.cell(200, 10, txt = n,
    		ln = 6, align = 'C')   
    pdf.cell(200, 10, txt = "2) Box plot ",
    		ln = 10, align = 'C')
    pdf.cell(200, 10, txt = "3) Cluster plot",
    		ln = 11, align = 'C')
    pdf.cell(200, 10, txt = "4) Correlation plot",
    		ln = 12, align = 'C')
    pdf.cell(200, 10, txt = "5) Histogram on each sample",
    		ln = 13, align = 'C')
    pdf.set_font("Arial", size = 10)
    z='There were a total of '+val3+' samples analysed under '+val4+' different conditions.'
    pdf.cell(200, 15, txt = z,
    		ln = 13, align = 'C')
    # save the pdf with name .pdf
    pdf.output("report.pdf")
def detect_distribution(data):
    f = Fitter(data)
    f.fit()
    f.summary()
def readfiles(l):
    #open all files and make them in a form that is readable
    if l[-3:]=='csv':
        #print('Opening csv file ',l)
        name=pd.read_csv(l,sep='\t') 
        namenp=np.array(name)       
    if l[-3:]=='lsx':
        #print('Opening excel file',l)
        name=pd.read_excel(l) 
        namenp=np.array(name)         
    return (name,namenp)
def pca(fil):
    print('PCA analysis on file: ',fil)
    g1=fil
    def readfiles(g1):
        # check for the names
        m=g1
        if m[-3:]=='lsx':
            name = pd.read_excel(io=m)
            name=name.fillna(0)
            namenp=np.array(name)
        if m[-3:]=='txt':
        #open all files and make them in a form that is readable
            name=pd.read_csv(m,sep='\t')
            namenp=np.array(name)
        return(name,namenp)
    ########################################################################
    print('.')
    print('.')
    print('The file you chose has the following columns:')
    print('.')
    print('.')
    m=readfiles(g1)[0].columns.values
    counter=0 
    for x in m :  
        print(counter,x)
        counter+=1  
    #########################################################################
    m=readfiles(g1)[0] 
    keepfrom = input(""" --> Type the columns you want to keep from:""")
    keepuntil= input(""" --> Type the columns you want to keep until(+1):""")  
    fr=int(keepfrom)
    un=int(keepuntil)
    i=m.iloc[:,fr:un]
    m=i.columns.values
    print('you have kept the following columns:')
    counter=0
    for x in m :  
        print(counter,x)
        counter+=1  
    dropfrom = input(""" --> Type the columns you want to drop from:""")
    dropuntil= input(""" --> Type the columns you want to drop until(+1):""")
    dfr=int(dropfrom)
    dun=int(dropuntil)
    ix=i.drop(i.columns[dfr:dun], axis=1)
    m=ix.columns.values
    print('you have kept the following columns:')
    count=0
    for x in m:
        count+=1
        print(count, x)
    numberList = []
    print(numberList)
    print("\n")
    for i in m:
        print("what is the index of",i, "[**the values must be numerical] :")
        item = int(input())
        numberList.append(item)
    print("User List is ", numberList)
    dirname=input("Name the folder that you want to save your plots:")
    os.mkdir(dirname)
    os.chdir(dirname)
    colum=[]
    for x in m:
        colum.append(x)
        count+=1
        print(count, x)
        plt.grid(b=True, which='both', color='lightgrey', linestyle='-', alpha=1)
        plt.hist(ix[x], color = np.random.random(3), edgecolor = 'black',bins = int(100))
        # seaborn histogram
        sns.distplot(ix[x], hist=True, kde=False, 
                     bins=int(100), color = 'blue',
                     hist_kws={'edgecolor':'black'})
        # Add labels
        l='Histogram of '+x
        plt.title(l)
        plt.xlabel('Values')
        plt.ylabel('Number of counts with the specif value')
        m='Histogram of '+x
        plt.savefig(m)
        plt.show()
        plt.clf()
    print(colum)  
    ax = plt.axes()
    ax.set_facecolor('whitesmoke')
    sns.boxplot(data=ix[colum],palette="Blues")  
    plt.savefig('Boxplot')
    plt.show()
    plt.clf()
    f, ax = plt.subplots(figsize=(10, 8))
    corr = ix.corr()
    bs=sns.heatmap(corr, mask=np.zeros_like(corr, dtype=np.bool), cmap='Blues',square=True, ax=ax,annot=True)
    bs.set_title('Correlation Matrix')
    plt.savefig('Correlation_matrix')
    plt.show()
    plt.clf()
    g = sns.clustermap(ix,cmap='mako')
    plt.title("Cluster-map")
    plt.savefig('Cluster')
    plt.show()
    plt.clf()
    print('You are done you can go and check your file ',dirname,' countaining all the polts!')
    ##############################################################################
    #                           MANUALLY                                         #
    ##############################################################################
    ls=numberList
    lst=pd.DataFrame(numberList) 
    ini=ix.T
    sns.set_style("white") 
    # Run The PCA
    pca = PCA(n_components=3)
    pca.fit(ini)
    # Store results of PCA in a data frame
    result=pd.DataFrame(pca.transform(ini), columns=['PCA%i' % i for i in range(3)], index=ini.index)
    X_pca=pca.transform(ini) 
    total_var = pca.explained_variance_ratio_.sum() * 100
    ex_variance=np.var(X_pca,axis=0)
    ex_variance_ratio = ex_variance/np.sum(ex_variance)
   # print (ex_variance_ratio) 
    al=len(ls)
    t=len(np.unique(ls))
    create_pdf(str(total_var),str(ex_variance_ratio),str(al),str(t))
   # print(total_var)
    # Plot initialisation
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    scatter=ax.scatter(result['PCA0'], result['PCA1'], result['PCA2'],alpha=0.7,c=ls, cmap="tab20", s=10)
    # make simple, bare axis lines through space:
    xAxisLine = ((min(result['PCA0']), max(result['PCA0'])), (0, 0), (0,0))
    ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r')
    yAxisLine = ((0, 0), (min(result['PCA1']), max(result['PCA1'])), (0,0))
    ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r')
    zAxisLine = ((0, 0), (0,0), (min(result['PCA2']), max(result['PCA2'])))
    ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r')
    # label the axes
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_zlabel("PC3")
    ax.set_title('pca')
    plt.legend()
    labels = np.unique(ls)
    handles = [plt.Line2D([],[],marker="o", ls="", 
                          color=scatter.cmap(scatter.norm(yi))) for yi in labels]
    plt.legend(handles, labels,title='Conditions',handlelength=0.4,fontsize='xx-small')
    f=g1[:-4]+'_pca.pdf'
    plt.savefig(f, bbox_inches='tight')
    plt.show()
print(  """
   ====================================================
  |   These are the files in your working directory:  |
  ===================================================== 
  
  """ )
for x in arr:
    if x[-3:]=='csv' or x[-3:]=='lsx':
        print(x)
counter=0
g1 = input("""Which file do you want to open :""")
counter+=1
readfiles(g1)
#more = input("""Are there more files you  would like to open y/n:""")
#if more[0]=='y':
 #   counter+=1
  #  global g2
   # g2 = input("""Which file do you want to open :""")
    #readfiles(g2) 
    #mor = input("""Are there more files you  would like to open y/n:""")
    #if mor[0]=='y':
        #counter+=1
        #global g3
        #g3 = input("""Which file do you want to open :""")
        #readfiles(g3) 
# keep all my files open 
if counter==1:
    global file1,file1np
    file1=readfiles(g1)[0]
    file1np=readfiles(g1)[1]
    #print('the first file has the following columns: ')
    #for col in file1.columns:
     #print(col)
elif  counter==2:
    file1=readfiles(g1)[0]
    file1np=readfiles(g1)[1]
    print('the first file has the following columns: ')
    for col in file1.columns:
        print(col)
    global file2,file2np
    file2=readfiles(g2)[0]
    file2np=readfiles(g2)[1]
    print('the second file has the following columns: ')
    for col in file2.columns:
        print(col)
elif  counter==3:
    file1=readfiles(g1)[0]
    file1np=readfiles(g1)[1]
    print('the first file has the following columns: ')
    for col in file1.columns:
        print(col)
    file2=readfiles(g2)[0]
    file2np=readfiles(g2)[1]
    print('the second file has the following columns: ')
    for col in file2.columns:
        print(col)
    global file3,file3np
    file3=readfiles(g3)[0]
    file3np=readfiles(g3)[1]
    print('the third file has the following columns: ')
    for col in file3.columns:
        print(col)
ploting=input("Would you like to continue with ploting y/n:")   
if ploting=='y':  
    pca(g1)
else :
    print('''Actually..... I don't really care...
          I will run this anyway !!''')
    print('''
          What are you gonna do about it ?
         ----------------------------------
          \                                                .
           \                                               .
            \                                              .    
                      __--- __                             .
                    _-       /--______                     .
               __--( /     \ )XXXXXXXXXXX\v.               .
             .-XXX(   O   O  )XXXXXXXXXXXXXXX-             .
            /XXX(       U     )        XXXXXXX\            .
          /XXXXX(              )--_  XXXXXXXXXXX\          .
         /XXXXX/ (      O     )   XXXXXX   \XXXXX\         .
         XXXXX/   /            XXXXXX   \__ \XXXXX         .
         XXXXXX__/          XXXXXX         \__---->        .
 ---___  XXX__/          XXXXXX      \__         /         .
   \-  --__/   ___/\  XXXXXX            /  ___--/=         .
    \-\    ___/    XXXXXX              '--- XXXXXX         .
       \-\/XXX\ XXXXXX                      /XXXXX         .
         \XXXXXXXXX   \                    /XXXXX/         .
          \XXXXXX      >                 _/XXXXX/          .
            \XXXXX--__/              __-- XXXX/            .
             -XXXXXXXX---------------  XXXXXX-             .
                \XXXXXXXXXXXXXXXXXXXXXXXXXX/               .
                  ""VXXXXXXXXXXXXXXXXXXV""                 .
          
          ''')
    pca(g1)

#detect_distribution(m)



