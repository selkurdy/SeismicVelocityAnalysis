# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 22:16:11 2016

@author: Sami Elkurdy

python svelan.py --datafilename testtvr.txt --datatype tvr --datacols 2 3 4 5 --listindata

python svelan.py --datafilename testtvr.txt --datatype tvr --datacols 2 3 4 5 --listquad

python svelan.py --datafilename testzva.txt --datatype zvav --datacols 4 5 0 3 --listindata

python svelan.py --datafilename testzva.txt --datatype zvav --datacols 4 5 0 3 --listquad

python svelan.py --datafilename tvrms.txt --datacols 2 3 4 5 --plot map --datatype tvr --listindata >t

python svelan.py --datafilename tvrms.txt --datacols 2 3 4 5 --plot every --datatype tvr --listindata --pickevery 200 >t
~ New version reading in segy
python svelansegy.py bokor_intvel.sgy --datatype segyzvi > t

~ for segy input:
    run program to generate zt output then run again using datatype zt to generate quad coefs
*******
~ step 1
python svelansegy.py bokor_intvel.sgy --datatype segyzvi
~ step 2
python svelansegy.py bokor_intvel_zt.txt --datatype zt
python svelansegy.py P7300_999_FINAL_STACKING_VEL_zvr.txt --datatype tvr
python svelansegy.py P7300_999_FINAL_STACKING_VEL_zvr.txt --datatype zvav
python svelansegy.py P7300_999_FINAL_STACKING_VEL_zvr.txt --datatype tvr --listquad


"""
import sys, os.path
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
# from scipy  import interpolate
# import scipy.stats as sts
# import pandas as pd
import segyio


def calctvr(t,vr):
#t is in t2w ms
    # t0=0.0
    # z0=0.0
    # allrec=[]
    # oneline=()
    reclen=len(t)

    z= []
    vi= []
    vav=[]
    vav.append(vr[0])
    t0sec=t[0]/1000.0
    z.append(t0sec * vav[0]/2.0)
    vi.append(vav[0])
    for i in range(1,reclen):
#        print"tvr t: %10.3f  vrms: %10.0f " %(tvr[i,0],tvr[i,1])
#        print"tvr t: %10.3f  vrms: %10.0f " %(t1wall[i],vrall[i])
        t1sec=t[i]/1000.0
        dt2w=t1sec-t0sec
        dv=vr[i]**2 * t1sec -vr[i-1]**2 * t0sec
        vi.append( np.sqrt(np.fabs(dv)/(dt2w)))
        dz= dt2w/2.0 * vi[i]
        z.append( z[i-1] +dz)
        vav.append(z[i]/t1sec * 2.0)
        t0sec=t1sec
    vi = np.array(vi)
    vav = np.array(vav)
    z = np.array(z)
    return vi,vav,z

def residuals(p,y,x):
    err = y -peval(x,p)
    return err

def peval(t,p):
    p[2]=0.0
    z = (t*p[1] + (t**2) * p[0])
    return z


def recordout(fp,rec,x,y,t,z):
    for i in range(t.size):
        print("%10d %12.2f  %12.2f  %10.0f  %10.0f " %(rec,x,y,t[i],z[i]),file=fp)


def pdlistqcoef(fp,id,x,y,tzqc,tzqclstsq,bsqa0,bsqa1,bsqa2,bslsqa0,bslsqa1,bslsqa2):
    """
    listing quadratic boostrap bs pvals 10 50 90 for each quadratic coeff
    """
    print("%10d  %12.2f  %12.2f  %10.3f  %10.3f  %10.3f  %10.3f  %10.3f  %10.3f  %10.3f  %10.3f  %10.3f  %10.3f  %10.3f  %10.3f %10.3f  %10.3f  %10.3f   %10.3f  %10.3f  %10.3f  %10.3f  %10.3f  %10.3f   %10.3f  %10.3f  %10.3f " %
    (id,x,y,
        tzqc[2],tzqc[1],tzqc[0],
        tzqclstsq[2],tzqclstsq[1],tzqclstsq[0],
        bsqa0[0],bsqa1[0],bsqa2[0],
        bsqa0[1],bsqa1[1],bsqa2[1],
        bsqa0[2],bsqa1[2],bsqa2[2],
        bslsqa0[0],bslsqa1[0],bslsqa2[0],
        bslsqa0[1],bslsqa1[1],bslsqa2[1],
        bslsqa0[2],bslsqa1[2],bslsqa2[2]),file=fp)



def listbyrecordzva(fp,x,y,z,v,listquad,pvals=[10,50,90],dp=2,bsiter=1000):
    xr=[]
    yr=[]
    rcounter=[]
    x0 = x[0]
    y0 = y[0]
    xr.append(x[0])
    yr.append(y[0])
    recordcounter=0
    rcounter.append(recordcounter)

    trecord=[]
    zrecord=[]

    if listquad:
        print('ID   X   Y   QA0   QA1   QA2   QLSQA0    QLSQA1    QLSQA2  ',end='',file=fp)
        print('QA0BSP{} QA1BSP{} QA2BSP{}  '.format(pvals[0],pvals[0],pvals[0]),end='',file=fp)
        print('QA0BSP{} QA1BSP{} QA2BSP{}  '.format(pvals[1],pvals[1],pvals[1]) ,end='',file=fp)
        print('QA0BSP{} QA1BSP{} QA2BSP{}  '.format(pvals[2],pvals[2],pvals[2]),end='',file=fp)
        print('QLSQA0BSP{}   QLSQA1BSP{}   QLSQA2BSP{} '.format(pvals[0],pvals[0],pvals[0]),end='',file=fp)
        print('QLSQA0BSP{}   QLSQA1BSP{}   QLSQA2BSP{} '.format(pvals[1],pvals[1],pvals[1]),end='',file=fp)
        print('QLSQA0BSP{}   QLSQA1BSP{}   QLSQA2BSP{} '.format(pvals[2],pvals[2],pvals[2]),file=fp)

        # print('ID   X   Y   QA0   QA1   QA2   QLSQA0    QLSQA1    QLSQA2   QA0BSP{} QA1BSP{} QA2BSP{}  QA0BSP{} QA1BSP{} QA2BSP{}  QA0BSP{} QA1BSP{} QA2BSP{}   QLSQA0BSP{}   QLSQA1BSP{}   QLSQA2BSP{}    QLSQA0BSP{}   QLSQA1BSP{}   QLSQA2BSP{}    QLSQA0BSP{}   QLSQA1BSP{}   QLSQA2BSP{}  '.
        #     format(pvals[0],pvals[0],pvals[0],
        #         pvals[1],pvals[1],pvals[1],
        #         pvals[2],pvals[2],pvals[2],
        #         pvals[0],pvals[0],pvals[0],
        #         pvals[1],pvals[1],pvals[1],
        #         pvals[2],pvals[2],pvals[2]),file=fp)
    for i in range(x.size):
        if x[i] == x0 and y[i] == y0:
            trecord.append( z[i]/v[i])
            zrecord.append(z[i])
        else:
            xr.append(x[i])
            yr.append(y[i])
            recordcounter +=1
            rcounter.append(recordcounter)
            tr=np.array(trecord)
            tr2ws = tr*2000.0
            zr=np.array(zrecord)
            tzqcoef = np.polyfit(tr,zr,2)
            p0=np.array(tzqcoef)
            plsq = leastsq(residuals,p0,args=(zr,tr), maxfev=2000)
            bsqa0,bsqa1,bsqa2,bslsqa0,bslsqa1,bslsqa2 =bsonly(tr,zr,n=bsiter,dp=dp,pvals=pvals)
            if listquad:
                pdlistqcoef(fp,recordcounter,x0,y0,tzqcoef,plsq[0],bsqa0,bsqa1,bsqa2,bslsqa0,bslsqa1,bslsqa2)
            else:
                recordout(fp,recordcounter,x0,y0,tr2ws,zr)
            x0=x[i]
            y0=y[i]
            trecord=[]
            zrecord=[]
            trecord.append( z[i]/v[i])
            zrecord.append(z[i])

    xr.append(x[i])
    yr.append(y[i])
    recordcounter +=1
    rcounter.append(recordcounter)
    tr=np.array(trecord)
    tr2ws = tr*2000.0
    zr=np.array(zrecord)
    tzqcoef = np.polyfit(tr,zr,2)
    p0=np.array(tzqcoef)
    plsq = leastsq(residuals,p0,args=(zr,tr), maxfev=2000)
    bsqa0,bsqa1,bsqa2,bslsqa0,bslsqa1,bslsqa2 =bsonly(tr,zr,n=bsiter,dp=dp,pvals=pvals)
    if listquad:
        pdlistqcoef(fp,recordcounter,x0,y0,tzqcoef,plsq[0],bsqa0,bsqa1,bsqa2,bslsqa0,bslsqa1,bslsqa2)
    else:
        recordout(fp,recordcounter,x0,y0,tr2ws,zr)




    xa=np.array(xr)
    ya=np.array(yr)
    rca=np.array(rcounter)
    return xa,ya,rca






def listbyrecordtvr(fp,x,y,t,vr,listquad,pvals=[10,50,90],dp=2,bsiter=1000):
    xr=[]
    yr=[]
    rcounter =[]
    x0 = x[0]
    y0 = y[0]
    recordcounter=0
    xr.append(x[0])
    yr.append(y[0])
    rcounter.append(recordcounter)

    trecord=[]
    vrrecord=[]
    if listquad:
        print('ID   X   Y   QA0   QA1   QA2   QLSQA0    QLSQA1    QLSQA2  ',end='',file=fp)
        print('QA0BSP{} QA1BSP{} QA2BSP{}  '.format(pvals[0],pvals[0],pvals[0]),end='',file=fp)
        print('QA0BSP{} QA1BSP{} QA2BSP{}  '.format(pvals[1],pvals[1],pvals[1]) ,end='',file=fp)
        print('QA0BSP{} QA1BSP{} QA2BSP{}  '.format(pvals[2],pvals[2],pvals[2]),end='',file=fp)
        print('QLSQA0BSP{}   QLSQA1BSP{}   QLSQA2BSP{} '.format(pvals[0],pvals[0],pvals[0]),end='',file=fp)
        print('QLSQA0BSP{}   QLSQA1BSP{}   QLSQA2BSP{} '.format(pvals[1],pvals[1],pvals[1]),end='',file=fp)
        print('QLSQA0BSP{}   QLSQA1BSP{}   QLSQA2BSP{} '.format(pvals[2],pvals[2],pvals[2]),file=fp)
        # print('ID   X   Y   QA0   QA1   QA2   QLSQA0     QLSQA1    QLSQA2   QA0BSP{} QA1BSP{} QA2BSP{}  QA0BSP{} QA1BSP{} QA2BSP{}  QA0BSP{} QA1BSP{} QA2BSP{} '.
        #     format(pvals[0],pvals[0],pvals[0],pvals[1],pvals[1],pvals[1],pvals[2],pvals[2],pvals[2]),file=fp)
    for i in range(x.size):
        if x[i] == x0 and y[i] == y0:
            trecord.append( t[i])
            vrrecord.append(vr[i])
        else:
            xr.append(x[i])
            yr.append(y[i])
            recordcounter +=1
            rcounter.append(recordcounter)
            vi,vav,z = calctvr(trecord,vrrecord)
            vrrecord=np.array(vrrecord)
            t1w=np.array(trecord)/ 2000.0
            # trec1w = trec/2000.0
            tzqcoef = np.polyfit(t1w,z,2)
            p0=np.array(tzqcoef)
            plsq = leastsq(residuals,p0,args=(z,t1w), maxfev=2000)
            bsqa0,bsqa1,bsqa2,bslsqa0,bslsqa1,bslsqa2 =bsonly(t1w,z,n=bsiter,dp=dp,pvals=pvals)
            if listquad:
                pdlistqcoef(fp,recordcounter,x0,y0,tzqcoef,plsq[0],bsqa0,bsqa1,bsqa2,bslsqa0,bslsqa1,bslsqa2)
            else:
                recordout(fp,recordcounter,x0,y0,t1w * 2000.0,z)
            x0=x[i]
            y0=y[i]
            trecord=[]
            vrrecord=[]
            trecord.append( t[i])
            vrrecord.append(vr[i])

    xr.append(x[i])
    yr.append(y[i])
    recordcounter +=1
    rcounter.append(recordcounter)
    # trec=np.array(trecord)
    t1w=np.array(trecord)/ 2000.0
    # trec1w = trec/2000.0
    vi,vav,z = calctvr(trecord,vrrecord)
    vrrecord=np.array(vrrecord)
    tzqcoef = np.polyfit(t1w,z,2)
    p0=np.array(tzqcoef)
    plsq = leastsq(residuals,p0,args=(z,t1w), maxfev=2000)
    bsqa0,bsqa1,bsqa2,bslsqa0,bslsqa1,bslsqa2 =bsonly(t1w,z,n=bsiter,dp=dp,pvals=pvals)
    if listquad:
        pdlistqcoef(fp,recordcounter,x0,y0,tzqcoef,plsq[0],bsqa0,bsqa1,bsqa2,bslsqa0,bslsqa1,bslsqa2)
    else:
        recordout(fp,recordcounter,x0,y0,t1w * 2000.0 ,z)




    xa=np.array(xr)
    ya=np.array(yr)
    rca=np.array(rcounter)
    return xa,ya,rca


def listbyrecordzt(fp,x,y,z,t,listquad,pvals=[10,50,90],dp=2,bsiter=1000):
    """
    t is already t2w in ms

    """


    xr=[]
    yr=[]
    rcounter=[]
    x0 = x[0]
    y0 = y[0]
    xr.append(x[0])
    yr.append(y[0])
    recordcounter=0
    rcounter.append(recordcounter)

    trecord=[]
    zrecord=[]
    if listquad:
        print('ID   X   Y   QA0   QA1   QA2   QLSQA0    QLSQA1    QLSQA2  ',end='',file=fp)
        print('QA0BSP{} QA1BSP{} QA2BSP{}  '.format(pvals[0],pvals[0],pvals[0]),end='',file=fp)
        print('QA0BSP{} QA1BSP{} QA2BSP{}  '.format(pvals[1],pvals[1],pvals[1]) ,end='',file=fp)
        print('QA0BSP{} QA1BSP{} QA2BSP{}  '.format(pvals[2],pvals[2],pvals[2]),end='',file=fp)
        print('QLSQA0BSP{}   QLSQA1BSP{}   QLSQA2BSP{} '.format(pvals[0],pvals[0],pvals[0]),end='',file=fp)
        print('QLSQA0BSP{}   QLSQA1BSP{}   QLSQA2BSP{} '.format(pvals[1],pvals[1],pvals[1]),end='',file=fp)
        print('QLSQA0BSP{}   QLSQA1BSP{}   QLSQA2BSP{} '.format(pvals[2],pvals[2],pvals[2]),file=fp)


        # print('ID   X   Y   QA0   QA1   QA2   QLSQA0     QLSQA1    QLSQA2   QA0BSP{} QA1BSP{} QA2BSP{}  QA0BSP{} QA1BSP{} QA2BSP{}  QA0BSP{} QA1BSP{} QA2BSP{} '.
        #     format(pvals[0],pvals[0],pvals[0],pvals[1],pvals[1],pvals[1],pvals[2],pvals[2],pvals[2]),file=fp)
    for i in range(x.size):
        if x[i] == x0 and y[i] == y0:
            trecord.append( t[i])
            zrecord.append(z[i])
        else:
            xr.append(x[i])
            yr.append(y[i])
            recordcounter +=1
            rcounter.append(recordcounter)
            tr=np.array(trecord)/2000.0 # convert it to t1w in sec
            tr2ws = tr*2000.0
            zr=np.array(zrecord)
            tzqcoef = np.polyfit(tr,zr,2)
            p0=np.array(tzqcoef)
            plsq = leastsq(residuals,p0,args=(zr,tr), maxfev=2000)
            bsqa0,bsqa1,bsqa2,bslsqa0,bslsqa1,bslsqa2 =bsonly(tr,zr,n=bsiter,dp=dp,pvals=pvals)
            if listquad:
                pdlistqcoef(fp,recordcounter,x0,y0,tzqcoef,plsq[0],bsqa0,bsqa1,bsqa2,bslsqa0,bslsqa1,bslsqa2)
            else:
                recordout(fp,recordcounter,x0,y0,tr2ws,zr)
            x0=x[i]
            y0=y[i]
            trecord=[]
            zrecord=[]
            trecord.append( t[i])
            zrecord.append(z[i])

    xr.append(x[i])
    yr.append(y[i])
    rcounter.append(recordcounter)
    tr=np.array(trecord)/2000.0 # convert it to t1w in sec
    tr2ws = tr*2000.0
    zr=np.array(zrecord)
    tzqcoef = np.polyfit(tr,zr,2)
    p0=np.array(tzqcoef)
    plsq = leastsq(residuals,p0,args=(zr,tr), maxfev=2000)
    bsqa0,bsqa1,bsqa2,bslsqa0,bslsqa1,bslsqa2 =bsonly(tr,zr,n=bsiter,dp=dp,pvals=pvals)
    if listquad:
        pdlistqcoef(fp,recordcounter,x0,y0,tzqcoef,plsq[0],bsqa0,bsqa1,bsqa2,bslsqa0,bslsqa1,bslsqa2)
    else:
        recordout(fp,recordcounter,x0,y0,tr2ws,zr)




    xa=np.array(xr)
    ya=np.array(yr)
    rca=np.array(rcounter)

    return xa,ya,rca









def plotallrecordszt(za,ta,listquad,plottitle = 'All Wells'):
    #make sure ta is t1w
    plt.gca().invert_yaxis()
    plt.axis([min(ta),max(ta),max(za),min(za)])
    plt.scatter(ta,za)

    if listquad:
        qc=np.polyfit(ta,za,2)
        ti=np.linspace(min(ta),max(ta))
        zi=np.polyval(qc,ti)
        plt.plot(ti,zi,c='r')

        p0=np.array(qc)
        plsq = leastsq(residuals,p0,args=(za,ta), maxfev=2000)
        zilsq= np.polyval(plsq[0],ti)
        plt.plot(ti,zilsq,c='g',lw=3)
        plt.annotate('z=%5.0f+%6.0f*t1w+%6.0f*t1w**2' % (qc[2],qc[1],qc[0]),\
            xy=(ta[4],za[4]),xytext=(0.5,0.85),\
            textcoords='figure fraction')

        plt.annotate('z=%5.0f+%6.0f*t1w+%6.0f*t1w**2' % (plsq[0][2],plsq[0][1],plsq[0][0]),\
            xy=(ta[4],za[4]),xytext=(0.5,0.80),\
            textcoords='figure fraction')


    #plt.title('All Wells ')
    plt.title(plottitle)
    plt.xlabel('T1W')
    plt.ylabel('Depth')
#    plt.text(0.5,2000,'z=%5.0f+%6.0f*t1w+%6.0f*t1w**2' % (qc[2],qc[1],qc[0]))
    plt.show()
    if listquad:
        print('%5.0f   %6.0f   %6.0f Quadratic' % (qc[2],qc[1],qc[0]))
        print('%5.0f   %6.0f   %6.0f Least Squares' % (plsq[0][2],plsq[0][1],plsq[0][0]))


def plotmap(x,y):
    plt.scatter(x,y,alpha= .8,s=10)
    plt.show()


def datain(fname,nheader,cols,zvmultiplier):
    xyzv = np.genfromtxt(fname,skip_header=nheader,usecols=cols)
    z =xyzv[:,2] * zvmultiplier[0]
    v = xyzv[:,3] * zvmultiplier[1]
    return xyzv[:,0],xyzv[:,1],z,v

def listindata(x,y,z,v,reduceby):
    recordcounter = 0
    rcounter=[]
    xr= []
    yr = []
    zr = []
    vr = []
    rcounter.append(recordcounter)
    x0,y0 = x[0],y[0]
    xr.append(x[0])
    yr.append(y[0])
    for i in range(x.size):
        if x[i] == x0 and y[i] == y0:
            if recordcounter % reduceby == 0 :
                # print('%10d  %12.2f  %12.2f  %10.0f  %10.0f'%(recordcounter,x[i],y[i],z[i],v[i]))
                zr.append(z[i])
                vr.append(v[i])
        else:
            x0,y0 = x[i] , y[i]
            recordcounter+=1
            if recordcounter % reduceby == 0 :
                # print('%10d  %12.2f  %12.2f  %10.0f  %10.0f'%(recordcounter,x[i],y[i],z[i],v[i]))
                zr.append(z[i])
                vr.append(v[i])
                rcounter.append(recordcounter)
                xr.append(x[i])
                yr.append(y[i])


        if recordcounter % reduceby == 0 :
            rcounter.append(recordcounter)
            zr.append(z[i])
            vr.append(v[i])
            xr.append(x[i])
            yr.append(y[i])


    return np.array(xr),np.array(yr),np.array(zr),np.array(vr),np.array(rcounter)



def process_zvi(tr,dz):
    """
    dz is vertical sampling in depth
    convert t to t2w
    output list of z t2w in ms
    """
    z = np.array(tr)
    t = np.array(tr)
    z[0] = 0
    t[0] = 0
    for i in range(1,tr.size):
        z[i] = i * dz
        t[i] =  t[i-1] +(dz / tr[i]) * 2000.0  # to make t2w in ms
    return z,t


def process_zva(tr,dz):
    """
    dz is vertical sampling in depth
    convert t to t2w
    output list of z t2w in ms
    """
    z = np.array(tr)
    t = np.array(tr)
    z[0] = 0
    t[0] = 0
    for i in range(1,tr.size):
        z[i] = i * dz
        t[i] = z[i] / tr[i] * 2000.0
    return z,t



def process_tvi(tr,dz):
    """
    dz is vertical sampling in t2w
    output list of z t2w in ms
    """

    z = np.array(tr)
    t = np.array(tr)
    dzsec = dz / 2000.0  #dz in sec t1w
    z[0] = 0
    t[0] = 0
    for i in range(1,tr.size):
        z[i] = z[i-1] + dzsec * tr[i]
        t[i] =  i * dz #  in t2w ms
    return z,t




def process_tvr(tr,dz):
    t = np.array(tr)
    vr = np.array(tr)
    for i in range(1,tr.size):
        t[i] = i * dz
        vr[i] =  tr[i]
    return t,vr


def list_zt(fp,z,t,x,y,trn,every = 100):
    for i in range(1,z.size):
        if i % every == 0:
            print('{:12.2f}   {:12.2f}   {:6.0f}   {:6.0f}'.format(x,y,z[i],t[i]),file=fp)



def list_tvr(fp,t,vr,x,y,trn,every = 100):
    for i in range(1,t.size):
        if i % every == 0:
            print('{:12.2f}   {:12.2f}   {:6.0f}   {:6.0f}'.format(x,y,t[i],vr[i]),file=fp)



def draw_bs_pairs_quadratic(x, y, size=1):
    """Perform pairs bootstrap for quadratic."""

    # Set up array of indices to sample from: inds
    inds = np.arange(len(x))

    # Initialize replicates: bs_slope_reps, bs_intercept_reps
    bs_a0_reps = np.empty(size)
    bs_a1_reps = np.empty(size)
    bs_a2_reps = np.empty(size)
    bs_a0lsq_reps = np.empty(size)
    bs_a1lsq_reps = np.empty(size)
    bs_a2lsq_reps = np.empty(size)

    # Generate replicates
    for i in range(size):
        bs_inds = np.random.choice(inds, len(inds))
        bs_x, bs_y = x[bs_inds], y[bs_inds]
        bs_a2_reps[i], bs_a1_reps[i], bs_a0_reps[i] = np.polyfit(bs_x, bs_y, 2)
        p0=np.array([bs_a2_reps[i],bs_a1_reps[i],bs_a0_reps[i]])
        plsq = leastsq(residuals,p0,args=(bs_y, bs_x), maxfev=2000)
        bs_a2lsq_reps[i] = plsq[0][0]
        bs_a1lsq_reps[i] = plsq[0][1]
        bs_a0lsq_reps[i] = plsq[0][2]

    return bs_a2_reps,bs_a1_reps, bs_a0_reps,bs_a2lsq_reps,bs_a1lsq_reps, bs_a0lsq_reps


def bsonly(x,y,n=1000,dp=2,pvals=[10,50,90]):
    """
    x should be t1w
    y should be z
    """


    bs_a2_reps,bs_a1_reps, bs_a0_reps,bs_a2lsq_reps,bs_a1lsq_reps, bs_a0lsq_reps = draw_bs_pairs_quadratic(x, y, size=n)

    #x is predicted, y is actual
    # print('Bootstrap: a0  {:.{prec}f} a1 {:.{prec}f}  a2 {:.{prec}f}'.format(icpt,s,prec=dp))
    bsa0 = np.array(pvals)
    bsa1 = np.array(pvals)
    bsa2 = np.array(pvals)
    bslsqa0 = np.array(pvals)
    bslsqa1 = np.array(pvals)
    bslsqa2 = np.array(pvals)



    for i in range(len(pvals)):
        bsa0[i]=np.percentile(bs_a0_reps, pvals[i])
        bsa1[i]=np.percentile(bs_a1_reps, pvals[i])
        bsa2[i]=np.percentile(bs_a2_reps, pvals[i])

        bslsqa0[i]=np.percentile(bs_a0lsq_reps, pvals[i])
        bslsqa1[i]=np.percentile(bs_a1lsq_reps, pvals[i])
        bslsqa2[i]=np.percentile(bs_a2lsq_reps, pvals[i])

        # print('Bootstrap: P{:<d} a0  {:.2f}  a1 {:.2f}  a2 {:.2f}'.format(pvals[i],bsa0[i],bsa1[i],bsa2[i]))
    return  bsa0,bsa1,bsa2,bslsqa0,bslsqa1,bslsqa2





def getcommandline():
    parser = argparse.ArgumentParser(description = 'Analyse Seismic Velocity segy or flat Files' )
    parser.add_argument('datafilename', help ='Either segy or Flat file with x y z[t] v[av,rms]')
    parser.add_argument('--datatype', choices=['segyzvi','segyzva','segytvr','segytvi','zvav','tvr','zt'],default='tvr',
            help='Data input options segy z vint or segy t2w vrms  or segy t2w vint or flatfile z vav or flatfile t2w in ms vrms or z t2w from segy')
    parser.add_argument('--datacols',nargs=4,type=int,default=[0,1,2,3],
            help='flat file x y z v column numbers. dfv= 0 1 2 3 ')
    parser.add_argument('--xhdr',type=int,default=181,help='xcoord header.default=181')
    parser.add_argument('--yhdr',type=int,default=185,help='ycoord header. default=185')
    parser.add_argument('--xyscalerhdr',type=int,default=71,help='hdr of xy scaler to divide by.default=71')

    parser.add_argument('--zvmultiplier',nargs=2,type=float,default=[1.0,1.0],help='z v multipliers. dfv= 1 1')
    parser.add_argument('--nheader',type=int,default=1,help='flat file header lines to skip. dfv=1')
    parser.add_argument('--minmaxz',type=float,nargs=2,default=[0.0,100000.0],help='Min Max of z[or t]. default= full data')
    parser.add_argument('--pickevery',type=int,default=50,help='reduce data by.dfv=50')
    parser.add_argument('--vertsample',type=int,default=100,help='vertical sampling. default = list every 100 samples')
    parser.add_argument('--listquad',action='store_true',default=False,help='List quadratic coeff. dfv= list zt databy record')
    parser.add_argument('--pvals',type=float,nargs=3,default=[10,50,90],
            help='3 values of P to compute corresponding Quadratic coefs, dfv= 10 50 90')
    parser.add_argument('--decimalplaces',type=int,default=2,help='decimal places. dfv=2')
    parser.add_argument('--bs_iter',type=int,default=1000,help='bootstrap iteration number.dfv=1000')
    parser.add_argument('--plot',choices=['map','all'],default='map',help='Plotting besides listing')

    result= parser.parse_args()
    if not result.datafilename :
        parser.print_help()
        exit()
    else:
        return result

def main():
    cmdl=getcommandline()
    in_fname = cmdl.datafilename
    dirsplit,fextsplit= os.path.split(cmdl.datafilename)
    fname,fextn= os.path.splitext(fextsplit)
    if cmdl.datatype == 'segyzvi':
        outfname = fname + '_zt.txt'
        with open(outfname,'w') as fp:

            with segyio.open( in_fname, "r" ,strict=False) as srcp:
                dz = srcp.bin[segyio.BinField.Interval]/1000.0
                alltraces = len(srcp.trace)
                for trnum,tr in enumerate(srcp.trace):
                    if trnum % cmdl.pickevery == 0 :
                        print('Processing trace # {} of {}'.format(trnum,alltraces),file=sys.stderr)
                        xysc = np.fabs(srcp.header[trnum][cmdl.xyscalerhdr])
                        xc = srcp.header[trnum][cmdl.xhdr]/ xysc
                        yc = srcp.header[trnum][cmdl.yhdr]/ xysc
                        z,t = process_zvi(tr,dz)
                        list_zt(fp,z,t,xc,yc,trnum,every = cmdl.vertsample)
                        #  from segy  list depth  time file
                print('Successfully generated {}'.format(outfname))

    elif cmdl.datatype == 'segyzva':
        outfname = fname + '_zt.txt'
        with open(outfname,'w') as fp:

            with segyio.open( in_fname, "r" ,strict=False) as srcp:
                dz = srcp.bin[segyio.BinField.Interval]/1000.0
                for trnum,tr in enumerate(srcp.trace):
                    if trnum % cmdl.pickevery == 0 :
                        print('Processing trace # {}'.format(trnum),file=sys.stderr)
                        xysc = np.fabs(srcp.header[trnum][cmdl.xyscalerhdr])
                        xc = srcp.header[trnum][cmdl.xhdr]/ xysc
                        yc = srcp.header[trnum][cmdl.yhdr]/ xysc
                        z,t = process_zva(tr,dz)
                        list_zt(fp,z,t,xc,yc,trnum,every = cmdl.vertsample)
                        #  from segy  list depth  time file
                print('Successfully generated {}'.format(outfname))



    elif cmdl.datatype == 'segytvr':
        outfname = fname + '_tvr.txt'
        with open(outfname,'w') as fp:

            with segyio.open( in_fname, "r" ,strict=False) as srcp:
                dz = srcp.bin[segyio.BinField.Interval]/1000.0
                for trnum,tr in enumerate(srcp.trace):
                    if trnum % cmdl.pickevery == 0 :
                        print('Processing trace # {}'.format(trnum),file=sys.stderr)
                        xysc = np.fabs(srcp.header[trnum][cmdl.xyscalerhdr])
                        xc = srcp.header[trnum][cmdl.xhdr]/ xysc
                        yc = srcp.header[trnum][cmdl.yhdr]/ xysc
                        t,vr = process_tvr(tr,dz)
                        list_tvr(fp,t,vr,xc,yc,trnum,every = cmdl.vertsample)
                        # the only one that outputs t2w vrms file
                print('Successfully generated {}'.format(outfname))


    elif cmdl.datatype == 'segytvi':
        outfname = fname + '_zt.txt'
        with open(outfname,'w') as fp:

            with segyio.open( in_fname, "r" ,strict=False) as srcp:
                dz = srcp.bin[segyio.BinField.Interval]/1000.0
                for trnum,tr in enumerate(srcp.trace):
                    if trnum % cmdl.pickevery == 0 :
                        print('Processing trace # {}'.format(trnum),file=sys.stderr)
                        xysc = np.fabs(srcp.header[trnum][cmdl.xyscalerhdr])
                        xc = srcp.header[trnum][cmdl.xhdr]/ xysc
                        yc = srcp.header[trnum][cmdl.yhdr]/ xysc
                        z,t = process_tvi(tr,dz)
                        list_zt(fp,z,t,xc,yc,trnum,every = cmdl.vertsample)
                        #  from segy  list depth  time file
                print('Successfully generated {}'.format(outfname))

    else:

        x,y,z,v =datain(cmdl.datafilename,cmdl.nheader,cmdl.datacols,cmdl.zvmultiplier)
        if cmdl.minmaxz[0] > z[0]:
            v = v[z>= cmdl.minmaxz[0]]
            x = x[z>= cmdl.minmaxz[0]]
            y = y[z>= cmdl.minmaxz[0]]
            z = z[z>= cmdl.minmaxz[0]]
        if cmdl.minmaxz[1] < z[-1]:
            v = v[z< cmdl.minmaxz[1]]
            x = x[z< cmdl.minmaxz[1]]
            y = y[z< cmdl.minmaxz[1]]
            z = z[z< cmdl.minmaxz[1]]


        if cmdl.datatype =='zvav':
            if cmdl.listquad:
                outfname = fname + '_qcf.txt'
            else:
                outfname = fname + '_t2wz.txt'
            with open(outfname,'w') as fp:
                xr,yr,rn=listbyrecordzva(fp,x,y,z,v,cmdl.listquad,cmdl.pvals,cmdl.bs_iter,cmdl.decimalplaces)
            print('Successfully generated {}'.format(outfname))
            if cmdl.plot == 'all':
                plt.scatter(z,v)
                plt.xlabel('Z')
                plt.ylabel('VAV')
                plt.show()
            elif cmdl.plot == 'map':
                plt.scatter(xr,yr,alpha= .8,s=10)
                plt.show()
        elif cmdl.datatype =='tvr': #z is t2w
            if cmdl.listquad:
                outfname = fname + '_qcf.txt'
            else:
                outfname = fname + '_t2wz.txt'
            with open(outfname,'w') as fp:
                xr,yr,rn=listbyrecordtvr(fp,x,y,z,v,cmdl.listquad,cmdl.pvals,cmdl.bs_iter,cmdl.decimalplaces)
            print('Successfully generated {}'.format(outfname))
            if cmdl.plot == 'all':
                plt.scatter(z,v)
                plt.xlabel('T2W ms')
                plt.ylabel('VRMS')
                plt.show()
            elif cmdl.plot == 'map':
                plt.scatter(xr,yr,alpha= .8,s=10)
                plt.show()

        elif cmdl.datatype =='zt': #zt output from segy file zvi
            if cmdl.listquad:
                outfname = fname + '_qcf.txt'
            else:
                outfname = fname + '_t2wz.txt'
            with open(outfname,'w') as fp:
                xr,yr,rn=listbyrecordzt(fp,x,y,z,v,cmdl.listquad,cmdl.pvals,cmdl.bs_iter,cmdl.decimalplaces)
            print('Successfully generated {}'.format(outfname))
            if cmdl.plot == 'all':
                plt.scatter(z,v)
                plt.xlabel('Z')
                plt.ylabel('T2W ms')
                plt.show()
            elif cmdl.plot == 'map':
                plt.scatter(xr,yr,alpha= .8,s=10)
                plt.show()


if __name__ == '__main__':
    main()
