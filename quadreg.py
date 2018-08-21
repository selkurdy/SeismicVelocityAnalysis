# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 23:22:19 2017

@author: Sami Elkurdy
Quadratic coefficients regular listing for import into Landmark or Paradigm

F:\sekData\Baronia\velocities>python quadreg.py --quadfilename quad.lst --regula
r 100 5100 200  > tdregular1.lst
F:\sekData\Baram\Velocities\Seismic Velocities>python quadreg.py --quadfilename
baramquad.lst --regular 100 5100 200 >tdregular.lst


**reading in 2 surfer ascii grids of a1 a2 null value has to be used
F:\sekData\Baronia\velocities>python quadreg.py --quadfilename t --quadtype surf
er2grid --surfer2gridnames a1a_trimmed.grd a2a_trimmed.grd --regular 100 5100 20
0 --nullvalue 1.70141e38

**reading in one xylist of a0 a1 a2 in the same file creating a xplot of t2w vs z
F:\sekData\Baronia\velocities>python quadreg.py --quadfilename quad_trimmed_300.
txt --quadtype xy1list --xy1listcols 1 2 6 7 8 --regular 100 5100 200 --plot >tr
egular300.txt

**out regular slices
F:\sekData\Baronia\velocities>python quadreg.py --quadtype xya2lists --xya2filen
ames a1_trimmed.dat a2_trimmed.dat --regular 100 5100 200 --regularslice
**Feb0516: functionality to read only 2 flat files for coefs and 1 with all coefs have been
properly tested. Surfer grid does not work!! reading three files has not been tested yet

**Feb0517: added depth conversion of horizons: Assumption:
    all dta has the same xy, i.e use gridding then export to xyz of a1, a2, and
    horizon time or depth to guarantee all having the same xy. I you need interpolation
    use ztqr.py ztqmc.py
    I filter input data to be less than 10000

python quadreg.py --quadtype xya1list --xya1filename P7300_999_FINAL_STACKING_VEL_zvr_qcf.txt --horizon --horfilename depth2300.txt --hordatacols 1 2 3
python quadreg.py --quadtype xya1list --xya1filename P7300_999_FINAL_STACKING_VEL_zvr_qcf.txt --regularslice


"""
import os.path
import argparse
import numpy as np
import matplotlib.pyplot as plt


def qcoefin(fname,qccol,nheader): #used for all three coefs in one file
    xyqc=np.genfromtxt(fname,usecols=qccol,skip_header =nheader)
    return xyqc[:,0],xyqc[:,1],xyqc[:,2],xyqc[:,3],xyqc[:,4]

def acoefin(fname,xyacols,nheader): #used for single coef per file
    xya=np.genfromtxt(fname,usecols=xyacols,skip_header=nheader)
    #filter surfer null values by taking all less than 10000, arbitrary!!
    xya = xya[xya[:,2]<10000.0]
    #xya = xya[~xya[:,2]==  missing]
    return xya[:,0],xya[:,1],xya[:,2]



def quadeval(cf,xi):
    z=cf[0]+ cf[1]*xi + cf[2] * xi * xi
    return z

def quadroots(cf0,cf1,cf2,xi):
    a0 =cf0- xi
    a1 = cf1
    a2 =cf2

    num=(a1 *a1) - 4.0 * a2 * a0
    if num >= 0.0:
        r1=(-a1 + np.sqrt(np.fabs(num)))/(2.0 * a2)
        r2 = (-a1 - np.sqrt(np.fabs(num)))/(2.0 * a2)
    else:
        r1= -a1/ (2.0 *a2)
        r2= np.sqrt(np.fabs(num))/(2.0 *a2)
    return r1,r2



def zconvlist(x,y,tregular,qc0,qc1,qc2,recid,toplot):
    t=np.arange(tregular[0],tregular[1],tregular[2])
    t1w = t /2000.0
    print('ID   X   Y   T2W   Z   VAV   QUADLSQA0     QUADLSQA1    QUADLSQA2')
    for j in range(x.size):
        zl=[]
        #recid +=j
        for i in range(t.size):
            zc=quadeval([qc0[j],qc1[j],qc2[j]],t1w[i])
            zl.append(zc)
            vav= zc/t1w[i]
            print("%10.0f  %12.2f  %12.2f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f" %\
            (recid+j,x[j],y[j],t[i],zc,vav,qc0[j],qc1[j],qc2[j]))
        if toplot:
            zi=np.array(zl)
            plt.scatter(t,zi,s=3)
    if toplot:
        plt.xlabel('T2W')
        plt.ylabel('DEPTH')
        plt.show()

def tconvlist(x,y,zregular,qc0,qc1,qc2,recid,toplot):
    z=np.arange(zregular[0],zregular[1],zregular[2])
    print('ID   X   Y   T2W   Z   VAV   QUADLSQA0     QUADLSQA1    QUADLSQA2')
    for j in range(x.size):
        #recid =+j
        tl =[]
        for i in range(z.size):
            t0,t1=quadroots(qc0[j],qc1[j],qc2[j],z[i])
            tl.append(t0)
            vav=z[i]/t0
            t2w= t0 * 2000.0
            print("%10.0f  %12.2f  %12.2f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  " %\
            (recid+j,x[j],y[j],t2w,z[i],vav,qc0[j],qc1[j],qc2[j]))
        if toplot:
            ti= np.array(tl)
            plt.scatter(ti,z,s=3)
    if toplot:
        plt.xlabel('T2W')
        plt.ylabel('DEPTH')
        plt.show()


def zconvslice(x,y,tregular,qc0,qc1,qc2,recid):
    t=np.arange(tregular[0],tregular[1],tregular[2])
    t1w = t /2000.0
    for i in range(t.size):
        #recid +=j
        slicefname = "time"+"%-.0f"% t[i] +'.txt'
        f = open(slicefname,'w')
        f.write('ID   X   Y   T2W   Z   VAV   QUADLSQA0     QUADLSQA1    QUADLSQA2 \n')
        for j in range(x.size):
            zc=quadeval([qc0[j],qc1[j],qc2[j]],t1w[i])
            vav= zc/t1w[i]
            f.write("%10.0f  %12.2f  %12.2f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n" %\
            (recid+j,x[j],y[j],t[i],zc,vav,qc0[j],qc1[j],qc2[j]))
        f.close()

def tconvslice(x,y,zregular,qc0,qc1,qc2,recid):
    z=np.arange(zregular[0],zregular[1],zregular[2])
    for i in range(z.size):
        #recid =+j
        slicefname = "depth"+"%-.0f"% z[i] +'.txt'
        f = open(slicefname,'w')
        f.write('ID   X   Y   T2W   Z   VAV   QUADLSQA0     QUADLSQA1    QUADLSQA2 \n')
        for j in range(x.size):
            t0,t1=quadroots(qc0[j],qc1[j],qc2[j],z[i])
            vav=z[i]/t0
            t2w= t0 * 2000.0
            f.write("%10.0f  %12.2f  %12.2f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n" %\
            (recid+j,x[j],y[j],t2w,z[i],vav,qc0[j],qc1[j],qc2[j]))
        f.close()


#depth convert horizons
def zconvhor(outfname,x,y,t,a0,a1,a2):
    with open(outfname,'w') as fp:
        t1w =t/2000.0
        print('X   Y   T2W   Z   VAV   QUADLSQA0     QUADLSQA1    QUADLSQA2',file=fp)
        for i in range(x.size):
            z=quadeval([a0[i],a1[i],a2[i]],t1w[i])
            vav=z/t1w[i]
            print("%12.2f  %12.2f  %10.0f  %10.0f   %10.0f  %10.0f  %10.0f  %10.0f " %\
             (x[i], y[i],t[i], z ,vav,a0[i],a1[i],a2[i]),file=fp)


#time convert horizons
def tconvhor(outfname,x,y,z,a0,a1,a2):
    with open(outfname,'w') as fp:
        print('X   Y   T2W   Z   VAV   QUADLSQA0     QUADLSQA1    QUADLSQA2',file=fp)
        for i in range(x.size):
            t0,t1=quadroots(a0[i],a1[i],a2[i],z[i])

            vav=z[i]/t0
            print("%12.2f  %12.2f  %10.0f  %10.0f   %10.0f  %10.0f  %10.0f  %10.0f " %\
             (x[i], y[i], t0 *2000.0,z[i],vav,a0[i],a1[i],a2[i]),file=fp)




def getcommandline():
    parser= argparse.ArgumentParser(description='Quadratic coefficients at regular increments and domain conversion of horizons')
    parser.add_argument('--quadtype',\
                        choices=['xya1list','xya2lists','xya3lists'],\
                                default='xya1list',help='quad coef input file type. dfv= xya1list')
    parser.add_argument('--xya1filename',help='Quadratic listing with x y a0 a1 a2 columns')
    parser.add_argument('--xya1listcols',type=int,nargs=5,default=[1,2,6,7,8],\
    help=' x y a0 a1 a2 data column numbers. dfv=1 2 6 7 8')
    parser.add_argument('--headerlinesa',type=int,default=1,help='Coef header lines to skip.default=1')
    parser.add_argument('--horizon',default=False,action='store_true',\
                        help='Domain convert hor not regular. dfv=false, i.e. regular')
    parser.add_argument('--headerlinesh',type=int,default=1,help='Horizon header lines to skip.default=1')
    parser.add_argument('--d2t',default=False,action='store_true',\
                        help='depth to time conversion, dfv=time to depth conversion')
    parser.add_argument('--horfilename',help='Horizons file name')
    parser.add_argument('--hordatacols',nargs=3,type=int,default=[0,1,2],\
                        help='x y t[z] column numbers. xy has to be the same as coefs. No interpolation. default= 0 1 2')
    parser.add_argument('--horscaleshift',nargs=2,type=float,default=[1.0,0.0],\
                        help='horizon multiplier then shift.dfv= 1.0 0.0')
    parser.add_argument('--startid',type=int,default=0,help='Start record id default=0')
    parser.add_argument('--regular',nargs=3,type=float,default=(100,4000,200),\
    help='Regular interpolation, start end increment.default = 100 4000 200')
    parser.add_argument('--regularz',action='store_true',default=False,\
    help='dfv=False, i.e. regularization in t2w')
    parser.add_argument('--regularslice',default=False,action='store_true',help='Generate slices not vertical functions')
    parser.add_argument('--plot',action='store_true',default=False,help='Plot td pairs ')
    parser.add_argument('--xya2filenames',nargs=2,help='a1 and a2 coef file names')
    parser.add_argument('--xya2listcols',nargs=6 ,type=int, default=(0,1,2,0,1,2),\
    help='x y a1 x y a2 in seperate files. dfv = 0 1 2 0 1 2')

    parser.add_argument('--xya3filenames',nargs=3,help='a0  a1 and a2 coef file names')
    parser.add_argument('--xya3listcols',nargs=9 ,type=int, default=(0,1,2,0,1,2,0,1,2),\
    help='x y a0 x y a1 x y a2 in seperate files. dfv = 0 1 2 0 1 2 0 1 2')

    result=parser.parse_args()


    if (not result.xya1filename and not result.xya2filenames )  :
        parser.print_help()
        exit()
    else:
        return result


def main():
    cmdl= getcommandline()
    if cmdl.horizon:

        xh,yh,th = acoefin(cmdl.horfilename,cmdl.hordatacols,cmdl.headerlinesh)
        th *= cmdl.horscaleshift[0]
        th += cmdl.horscaleshift[1]
        dirsplit,fextsplit= os.path.split(cmdl.horfilename)
        fname,fextn= os.path.splitext(fextsplit)
        if cmdl.d2t:
            horconvfile = os.path.join(dirsplit,fname) +"_z2t.txt"
        else:
            horconvfile = os.path.join(dirsplit,fname) +"_t2z.txt"
    if cmdl.quadtype =='xya1list':
        x,y,qc0,qc1,qc2 = qcoefin(cmdl.xya1filename,cmdl.xya1listcols,cmdl.headerlinesa)
        #i.e. time convert and z is regular
        if cmdl.horizon:
            #horizon depth to time conversion
            if cmdl.d2t:
                tconvhor(horconvfile,x,y,th,qc0,qc1,qc2)
            #horizon depth to time conversion
            else:
                zconvhor(horconvfile,x,y,th,qc0,qc1,qc2)

        elif cmdl.regularz:
            if cmdl.regularslice:
                tconvslice(x,y,cmdl.regular,qc0,qc1,qc2,cmdl.startid)
            else:
                tconvlist(x,y,cmdl.regular,qc0,qc1,qc2,cmdl.startid,cmdl.plot)
        #i.e. depth convert and t is regular
        else:
            if cmdl.regularslice:
                zconvslice(x,y,cmdl.regular,qc0,qc1,qc2,cmdl.startid)
            else:
                zconvlist(x,y,cmdl.regular,qc0,qc1,qc2,cmdl.startid,cmdl.plot)



#2 coef file listings representing a1 a2
    elif cmdl.quadtype== 'xya2lists':
        x1,y1,a1=acoefin(cmdl.xya2filenames[0],cmdl.xya2listcols[:3],cmdl.headerlinesa)
        x2,y2,a2=acoefin(cmdl.xya2filenames[1],cmdl.xya2listcols[3:],cmdl.headerlinesa)
        a0= np.zeros_like(x1) #generate a array of zeros of same size as x1
        if cmdl.horizon:
            #horizon depth to time conversion
            if cmdl.d2t:
                tconvhor(horconvfile,x,y,th,qc0,qc1,qc2)
            #horizon depth to time conversion
            else:
                zconvhor(horconvfile,x,y,th,qc0,qc1,qc2)
        elif cmdl.regularz:
            if cmdl.regularslice:
                tconvslice(x1,y1,cmdl.regular,a0,a1,a2,cmdl.startid)
            else:
                tconvlist(x1,y1,cmdl.regular,a0,a1,a2,cmdl.startid,cmdl.plot)
        else:
            if cmdl.regularslice:
                zconvslice(x1,y1,cmdl.regular,a0,a1,a2,cmdl.startid)
            else:
                zconvlist(x1,y1,cmdl.regular,a0,a1,a2,cmdl.startid,cmdl.plot)
#3 coef files representing a0 a1 a2
    elif cmdl.quadtype== 'xya3lists':
        x0,y0,a0=acoefin(cmdl.xya3filenames[0],cmdl.xya3cols[:3],cmdl.headerlinesa)
        x1,y1,a1=acoefin(cmdl.xya3filenames[1],cmdl.xya3cols[3:6],cmdl.headerlinesa)
        x2,y2,a2=acoefin(cmdl.xya3filenames[2],cmdl.xya3cols[6:8],cmdl.headerlinesa)
        if cmdl.horizon:
            #horizon depth to time conversion
            if cmdl.d2t:
                tconvhor(horconvfile,x,y,th,qc0,qc1,qc2)
            #horizon depth to time conversion
            else:
                zconvhor(horconvfile,x,y,th,qc0,qc1,qc2)
        if cmdl.regularz:
            if cmdl.regularslice:
                tconvslice(x0,y0,cmdl.regular,a0,a1,a2,cmdl.startid)
            else:
                tconvlist(x0,y0,cmdl.regular,a0,a1,a2,cmdl.startid,cmdl.plot)
        else:
            if cmdl.regularslice:
                zconvslice(x0,y0,cmdl.regular,a0,a1,a2,cmdl.startid)
            else:
                zconvlist(x0,y0,cmdl.regular,a0,a1,a2,cmdl.startid,cmdl.plot)





if __name__ == '__main__':
    main()

