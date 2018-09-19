import os
import numpy
from sys import exit
import numpy.ma as ma

'''
Author: Prof. P. Lewis, UCL
email: p.lewis@ucl.ac.uk

'''

class Kernels:
    '''
    Linear kernel models
    '''
    def __init__(self,vza,sza,raa,critical=1,RossHS=True,RecipFlag=True,HB=2.0,BR=1.0,MODISSPARSE=True,MODISDENSE=False,RossType='Thick',normalise=1,normalize=0,LiType='Transit',doIntegrals=True,BSAangles=[],nbar=0.0):
        '''
        The class creator sets up the kernels for some angle set. Default Li is MODISSPARSE parameter set
        The kernels are accessible from:
                self.Isotropic
                self.Ross
                self.Li
        The angles are accesible from:
                self.vza (or self.vzaDegrees)
                self.sza (or self.szaDegrees)
                self.raa (or self.raaDegrees)
                N.B. Hot spot direction is vza == sza and raa = 0.0
        Kernels integrals are acessible from:
                self.BSAangles (angles in degrees)
                self.BSA_Isotropic (directional-hemispherical integral of self.Isotropic)
                self.BSA_Ross (directional-hemispherical integral of self.Ross)
                self.BSA_Li (directional-hemispherical integral of self.Li)
                self.WSA_Isotropic (bi-hemispherical integral of self.Isotropic)
                self.WSA_Ross (bi-hemispherical integral of self.Ross)
                self.WSA_Li (bi-hemispherical integral of self.Li)
                N.B. You need to set the doIntegrals flag to True on creating an instance of the kernels class if you 
                want access to integrals. The processing takes a bit of time.
        Printing methods are available:
                self.printIntegrals(header=True,reflectance=False)              
                self.printKernels(header=True,reflectance=False)

        Required parameters:

        @param vza: an array containg view zenith angles in degrees
        @param sza: an array containing solar zenith angles in degrees
        @param raa: an array containg relative azimuth angles in degrees

        Options:
        @option critical=1: set to 1 to exit on error, 0 not to
        @option RecipFlag=True: Li reciprocal flag
        @option HB: Li kernel parameter HB 
        @option BR: Li kernel parameter
        @option MODISSPARSE: set to True for default MODIS Li Sparse parameters (overrides BR and HB to 2.0 and 1.0)
        @option MODISDENSE: set to True for default MODIS Li Dense parameters (override BR and HB to 2.0 and 2.5)
        @option RossType: set to 'Thin' for Ross Thin (default) else 'Thick'
        @option LiType: set to 'Sparse' for LiSparse (default). Other options: 'Roujean', 'Dense'
        @option normalise: set to 1 to make kernels 0 at nadir view illumination (default), set to 0 for no normalisation (can also use US spelling, i.e. normalize)
        @option doIntegrals: set to True to calculate integrals of kernels numerically. Set to False not to calculate them. At some point will have Approx flag here as well.
        @option BSAangles: solar zenith angles at which to calculate directional-hemispherical integral of kernels (default 0-89 in steps of 1 degree). Units: degrees.
        @option nbar: the sza at which the isotropic term is set to if normalise=1 is turned on (default 0)

        Notes:
        Requires numpy. If you do integrals, this also requires scipy (or rather scipy.integrate)
        If you want to mimic the results in Wanner et al. 1995, I've set a special function called self.mimic at the end here.
        
        
        '''
        self.__setup(critical=critical,RecipFlag=RecipFlag,RossHS=RossHS,HB=HB,BR=BR,MODISSPARSE=MODISSPARSE,MODISDENSE=MODISDENSE,RossType=RossType,normalise=normalise,normalize=normalize,LiType=LiType,doIntegrals=doIntegrals,BSAangles=BSAangles,nbar=nbar)
        self.shape = vza.shape
        mask = ma.getmask(vza)
        if mask.shape != self.shape:
            mask = (numpy.zeros(vza.shape).astype(bool) + mask).astype(bool)
        VZA = numpy.array(vza[~mask])
        SZA = numpy.array(sza[~mask])
        RAA = numpy.array(raa[~mask])
        self.setAngleInfo(VZA,SZA,RAA)
        self.__doKernels()
        self.__postProcess()
        R = self.Ross.copy()
        L = self.Li.copy()
        self.Ross = numpy.zeros(self.shape)
        self.Ross[~mask] = R
        self.Li = numpy.zeros(self.shape)
        self.Li[~mask] = L
        self.Ross = ma.array(self.Ross,mask=mask,hard_mask=True,copy=True,fill_value=0.)
        self.Li = ma.array(self.Li,mask=mask,hard_mask=True,copy=True,fill_value=0.)

 
    def __setup(self,critical=1,RecipFlag=True,RossHS=True,HB=2.0,BR=1.0,MODISSPARSE=True,MODISDENSE=False,RossType='Thick',normalise=1,normalize=0,LiType='Sparse',doIntegrals=True,BSAangles=[],nbar=0.0):
        self.nbar = nbar
        self.__NEARLYZERO = 1e-20
        self.critical = critical
        self.FILE = -1
        self.outputFile = ''
        # kernel options etc.
        self.LiType = LiType
        self.RossHS = RossHS
        self.doIntegrals = doIntegrals
        if MODISDENSE == True:
            LiType = 'Dense'
            self.HB = 2.0
            self.BR = 2.5
        else:
            if MODISSPARSE == True:
                LiType = 'Sparse'
                self.HB = 2.0
                self.BR = 1.0
            else:
                self.HB = HB
                self.BR = BR
        #self.LiType = LiType
        self.RossType = RossType
        self.normalise = normalise
        self.RecipFlag = RecipFlag
        # some useful numbers
        self.__M_PI = numpy.pi
        self.__M_PI_2 = self.__M_PI * 0.5
        self.__M_PI_4 = self.__M_PI * 0.25
        self.__M_1_PI = 1.0/self.__M_PI

        self.normalise = 0
        self.__integrateKernels(BSAangles=BSAangles)

        if (normalise >= 1 or normalize >= 1):
            self.normalise = max(normalise,normalize)

    def __postProcess(self):
        '''
        Private method for dealing with normalisation
        '''
        self.LiNorm = 0.
        self.RossNorm = 0.
        self.IsotropicNorm = 0.
        # if we are normalising the last element of self.Isotropic, self.Ross and self.Li  contain the nadir-nadir kernel
        if self.normalise >= 1:
            # normalise nbar-nadir (so kernel is 0 at nbar-nadir)
            self.RossNorm = self.Ross[-1]
            self.LiNorm = self.Li[-1]
            self.Ross = self.Ross - self.RossNorm
            self.Li  = self.Li - self.LiNorm
            # depreciate length of arrays (well, teh ones we'll use again in any case)
            self.Ross = self.Ross[0:-1]
            self.Li = self.Li[0:-1]
            self.Isotropic = self.Isotropic[0:-1]
            self.vzaDegrees = self.vzaDegrees[0:-1]
            self.szaDegrees = self.szaDegrees[0:-1]
            self.raaDegrees = self.raaDegrees[0:-1]
            self.N = len(self.vzaDegrees)
            self.vza = self.vza[0:-1]
            self.sza = self.sza[0:-1]
            self.raa = self.raa[0:-1]

    def __doKernels(self):
        '''
        Private method to run the various kernel methods
        '''
        # the kernels
        self.IsotropicKernel()
        self.RossKernel()
        self.LiKernel()

    def setAngleInfo(self,vza,sza,raa):
        '''
        Private method to store and organise the input angle data
        '''
        self.vzaDegrees = numpy.array([vza]).flatten()
        self.szaDegrees = numpy.array([sza]).flatten()
        self.raaDegrees = numpy.array([raa]).flatten()
        self.N = len(self.vzaDegrees)
        
        if(self.N != len(self.szaDegrees) or self.N != len(self.raaDegrees)):
            self.error('kernels: inconsistent number of samples in vza, sza and raa data: ' + str(len(self.vzaDegrees)) + ', ' + str(len(self.szaDegrees)) + ', ' + str(len(self.raaDegrees)),critical=self.critical)
            print(self.vzaDegrees)
            print(self.szaDegrees)
            print(self.raaDegrees)
            return [-1]
        
        if (self.normalise >= 1):
            # calculate nadir term by extending array
            self.vzaDegrees = numpy.array(list(self.vzaDegrees) + [0.0]).flatten()
            self.szaDegrees = numpy.array(list(self.szaDegrees) + [self.nbar]).flatten()
            self.raaDegrees = numpy.array(list(self.raaDegrees) + [0.0]).flatten()
            # not N is one too many now
            self.N = len(self.vzaDegrees)

        self.vza = self.dtor(self.vzaDegrees)
        self.sza = self.dtor(self.szaDegrees) # -1 to make HS direction for raa = 0
        self.raa = self.dtor(self.raaDegrees)
        w = numpy.where(self.vza < 0)[0]
        self.vza[w] = -self.vza[w]
        self.raa[w] = self.raa[w] + numpy.pi
        w = numpy.where(self.sza < 0)[0]
        self.sza[w] = -self.sza[w]
        self.raa[w] = self.raa[w] + numpy.pi


    def __integrateKernels(self,BSAangles=[]):
        '''
        Private method to call integration functions for the kernels


         NB - this overwrites all kernel info ... so be careful how/where you call it
        @option: BSAangles=[] allows the user to set the sza angles at which directional-hemispherical intergal is calculated, else steps of 1 degree from 0 to 89 (though I wouldnt trust it down to 90)
        This function can be rather slow, so using fewer samples or an approximate function may be a god idea
        '''
        if (self.doIntegrals == False):
            return;
        
        import scipy.integrate
        if BSAangles == []:
            BSAangles = numpy.array(range(90))*1.0

        self.BSAangles = numpy.array(BSAangles).flatten()

        # isotropic integral
        self.BSA_Isotropic = numpy.zeros(len(self.BSAangles))+1.0
        self.BSA_Ross = numpy.zeros(len(self.BSAangles))
        self.BSA_Li = numpy.zeros(len(self.BSAangles))
        self.BSA_Isotropic_error = numpy.zeros(len(self.BSAangles))
        self.BSA_Ross_error = numpy.zeros(len(self.BSAangles))
        self.BSA_Li_error = numpy.zeros(len(self.BSAangles))
        
        i = 0
        mu = numpy.cos(self.BSAangles*numpy.pi/180.)
        for sza in self.BSAangles:
            # ross integral
            self.BSA_Ross[i], self.BSA_Ross_error[i] = scipy.integrate.dblquad(RossFunctionForIntegral,0.0, 1.0, __gfun, __hfun, args=(sza,self))
            self.BSA_Li[i], self.BSA_Li_error[i] = scipy.integrate.dblquad(LiFunctionForIntegral,0.0, 1.0, __gfun, __hfun, args=(sza,self))
            i = i + 1
        self.WSA_Ross =  -2.0 * scipy.integrate.simps(self.BSA_Ross * mu,mu)
        self.WSA_Li =  -2.0 * scipy.integrate.simps(self.BSA_Li * mu,mu)
        return
        
    def __GetPhaang(self):
        '''
        Private method to calculate Phase angle component of kernel
        '''
        self.__cosphaang = self.__cos1*self.__cos2 + self.__sin1*self.__sin2*self.__cos3
        # better check the bounds before arccos ... just to be safe
        w = numpy.where(self.__cosphaang < -1)[0]
        self.__cosphaang[w] = -1.0
        w = numpy.where(self.__cosphaang > 1)[0]
        self.__cosphaang[w] = 1.0       
        self.__phaang = numpy.arccos(self.__cosphaang)
        self.__sinphaang = numpy.sin(self.__phaang)
        return

    def __RossKernelPart(self):
        '''
        Private method to calculate main part of Ross kernel
        '''
        self.__cos1 = numpy.cos(self.vza)
        self.__cos2 = numpy.cos(self.sza)
        
        self.__sin1 = numpy.sin(self.vza)
        self.__sin2 = numpy.sin(self.sza)
        self.__cos3 = numpy.cos(self.raa)
        self.__GetPhaang()
        self.rosselement = (self.__M_PI_2 - self.__phaang)*self.__cosphaang+self.__sinphaang
        return

    def GetDistance(self):
        '''
        Private method to get distance component of Li kernels
        '''
        temp = self.__tan1*self.__tan1+self.__tan2*self.__tan2-2.*self.__tan1*self.__tan2*self.__cos3;
        w = numpy.where(temp < 0)[0]
        temp[w] = 0.0
        self.__temp = temp # used by other functions ??
        distance = numpy.sqrt(temp)
        return distance

    def GetpAngles(self, tan1):
        '''
        Private method to do B/R transformation for ellipse shape
        '''
        t = self.BR * tan1
        w = numpy.where( t < 0.)[0]
        t[w] = 0.0
        angp = numpy.arctan(t)
        s = numpy.sin(angp)
        c = numpy.cos(angp)
        # have to make sure c isnt 0
        w = numpy.where(c == 0)[0]
        c[w] = self.__NEARLYZERO
        return c,s,t

    def GetOverlap(self):
        '''
        Private method to do HB ratio transformation
        '''
        self.__temp = 1./self.__cos1 + 1./self.__cos2

        self.__cost =  self.HB * numpy.sqrt(self.__distance * self.__distance + self.__tan1 * self.__tan1 * self.__tan2 * self.__tan2 * self.__sin3 * self.__sin3) / self.__temp;
        w = numpy.where(self.__cost < -1)[0]
        self.__cost[w] = -1.0
        w = numpy.where(self.__cost > 1.0)[0]
        self.__cost[w] = 1.0
        self.__tvar = numpy.arccos(self.__cost)
        self.__sint = numpy.sin(self.__tvar)
        self.__overlap = self.__M_1_PI * (self.__tvar - self.__sint * self.__cost) * self.__temp
        w = numpy.where(self.__overlap < 0)[0]
        self.__overlap[w] = 0.0
        return
         
    def RoujeanKernel(self):
        '''
        Private method - call to calculate Roujean shadowing kernel
        '''
        # first make sure its in range 0 to 2 pi
        self.__phi = numpy.abs((self.raa % (2.*numpy.pi)))
        self.__cos3 = numpy.cos(self.__phi)
        self.__sin3 = numpy.sin(self.__phi)
        self.__tan1 = numpy.tan(self.sza)
        self.__tan2 = numpy.tan(self.vza)

        self.__distance = self.GetDistance()
        self.Li = 0.5 * self.__M_1_PI * ((self.__M_PI - self.__phi) * self.__cos3 + self.__sin3) * self.__tan1 * self.__tan2 - self.__M_1_PI * (self.__tan1 + self.__tan2 + self.__distance);
        return

    def LiKernel(self):
        '''
        Private method - call to calculate Li Kernel
        '''
        # at some point add in LiGround kernel & LiTransit
        if self.LiType == 'Roujean':
            return self.RoujeanKernel()
        # first make sure its in range 0 to 2 pi
        self.__phi = numpy.abs((self.raa % (2.*numpy.pi)))
        self.__cos3 = numpy.cos(self.__phi)
        self.__sin3 = numpy.sin(self.__phi)
        self.__tanti = numpy.tan(self.sza)
        self.__tantv = numpy.tan(self.vza)
        self.__cos1, self.__sin1, self.__tan1 = self.GetpAngles(self.__tantv);
        self.__cos2, self.__sin2, self.__tan2 = self.GetpAngles(self.__tanti);
        self.__GetPhaang(); # sets cos & sin phase angle terms 
        self.__distance = self.GetDistance(); # sets self.temp
        self.GetOverlap(); # also sets self.temp
        if self.LiType == 'Sparse':
            if self.RecipFlag == True:
                self.Li = self.__overlap - self.__temp + 0.5 * (1. + self.__cosphaang) / self.__cos1 / self.__cos2;
            else:
                self.Li = self.__overlap - self.__temp + 0.5 * (1. + self.__cosphaang) / self.__cos1;
        else:
            if self.LiType == 'Dense':
                if self.RecipFlag:
                    self.Li = (1.0 + self.__cosphaang) / (self.__cos1 * self.__cos2 * (self.__temp - self.__overlap)) - 2.0;
                else:
                    self.Li = (1.0 + self.__cosphaang) / (self.__cos1 * (self.__temp - self.__overlap)) - 2.0;
            else:
                B = self.__temp - self.__overlap
                w = numpy.where(B <= 2.0)
                self.Li = B*0.0
                if self.RecipFlag == True:
                    Li = self.__overlap - self.__temp + 0.5 * (1. + self.__cosphaang) / self.__cos1 / self.__cos2;
                else:
                    Li = self.__overlap - self.__temp + 0.5 * (1. + self.__cosphaang) / self.__cos1;
                self.Li[w] = Li[w]

                w = numpy.where(B > 2.0)
                if self.RecipFlag:
                    Li = (1.0 + self.__cosphaang) / (self.__cos1 * self.__cos2 * (self.__temp - self.__overlap)) - 2.0;
                else:
                    Li = (1.0 + self.__cosphaang) / (self.__cos1 * (self.__temp - self.__overlap)) - 2.0;
                self.Li[w] = Li[w]
        return
        
    def IsotropicKernel(self):
        '''
        Public method - call to calculate Isotropic kernel
        '''
        # default behaviour
        self.Isotropic = numpy.zeros(self.N)+1.0
        return

    def RossThin(self):
        '''
        Public method - call to calculate RossThin kernel
        '''
        self.__RossKernelPart()
        self.rosselement = self.rosselement/(self.__cos1*self.__cos2)
        return;

    def RossThick(self):
        '''
        Public method - call to calculate RossThick kernel
        '''
        self.__RossKernelPart()
        self.rosselement = self.rosselement/(self.__cos1+self.__cos2)
        return;

    def RossKernel(self):
        '''
        Public method - call to calculate Ross Kernel
        '''
        if self.RossType == 'Thin':
            self.RossThin()
        else:
            self.RossThick()
        self.Ross = self.rosselement
        if self.RossHS != False:
            if self.RossHS == True:
                self.RossHS = 0.25
            self.Ross = self.Ross * (1 + 1/(1 + self.__phaang/self.RossHS))
    

    def dtor(self,x):
        '''
        Public method to convert degrees to radians
        '''
        return x*numpy.pi/180.0

    def rtod(self,x):
        '''
        Public method to convert radians to degrees
        '''
        return x*180./numpy.pi
    
    def error(self,msg,critical=0,newline=1,code=-1):
        '''
        Public method to do Class error reporting
        @param msg: error message
        @param critical: set to 1 if require exit (default critical=0)
        @param newline: set to 0 if newline not required (default newline=0)
        @param code: error code reported on exit if critical error (default code=-1)
        '''
        if newline == 1:
            nl = '\n'
        else:
            nl = ''
        print(msg + nl)
        if critical == 1:
            exit([code])


    def printIntegrals(self,header=True,reflectance=False):
        '''
        Public method to print kernel integrals (to stdout only at present)
        '''
        if(header == True):
            self.printer('# ' + str(self.N) + ' samples Ross: ' + self.RossType + ' Li: ' + self.LiType + ' Reciprocal: ' + str(self.RecipFlag) + ' normalisation: ' + str(self.normalise) + ' HB ' + str(self.HB) + ' BR ' + str(self.BR) + '\n');
            self.printer('# WSA: Isotropic 1.0 Ross ' + str(self.WSA_Ross) + ' Li ' + str(self.WSA_Li))
            self.printer('# 1: SZA (degrees) 2: BSA Isotropic 3: BSA Ross 4: BSA Li')
            if (reflectance == True):
                self.printer(' ');
            self.printer('\n');
            
        for i in range(len(self.BSAangles)):
            self.printer(str(self.BSAangles[i]) +  ' ' + str(self.BSA_Isotropic[i]) + ' ' + str(self.BSA_Ross[i]) + ' ' + str(self.BSA_Li[i]))
            # print refl data if wanted
            if (reflectance == True):
                self.printer(' ');
            self.printer('\n');
        return

    def printKernels(self,header=True,reflectance=False,file=False):
        '''
        Public method to print kernel values (to stdout only at present)        
        '''
        if(file != False):
            if(file != self.outputFile and self.FILE != -1):
                self.FILE.close()
            self.outputFile = file
            self.FILE = open(self.outputFile,'w')

        if(header == True):
            self.printer('# ' + str(self.N) + ' samples Ross: ' + self.RossType + ' Li: ' + self.LiType + ' Reciprocal: ' + str(self.RecipFlag) + ' normalisation: ' + str(self.normalise) + ' HB ' + str(self.HB) + ' BR ' + str(self.BR) + '\n');
            self.printer('# 1: VZA (degrees) 2: SZA (degrees) 3: RAA (degrees) 4: Isotropic 5: Ross 6: Li')
            if (reflectance == True):
                self.printer(' ');
            self.printer('\n');
            
        for i in range(self.N):
            self.printer(str(self.vzaDegrees[i]) + ' ' + str(self.szaDegrees[i]) + ' ' + str(self.raaDegrees[i]) + ' ' + str(self.Isotropic[i]) + ' ' + str(self.Ross[i]) + ' ' + str(self.Li[i]))
            # print refl data if wanted
            if (reflectance == True):
                self.printer(' ');
            self.printer('\n');
        return

    def printer(self,msg):
            '''
            Public print method ... make more flexible eg for printing to files at some point
            ''' 
            if (self.FILE == -1):
                print(msg),
            else:
                self.FILE.write(msg)


# some things required for the numerical integration

def _Kernels__gfun(x):
    return 0.0

def _Kernels__hfun(x):
    return 2.0*numpy.pi

def RossFunctionForIntegral(phi,mu,sza,self):
    #print phi
    #print mu
    #print sza
    #print '========'
    vza = numpy.arccos(mu)
    raa = self.rtod(phi)
    self.setAngleInfo(vza,sza,raa)
    self.RossKernel()
    return mu * self.Ross[0] / numpy.pi


def LiFunctionForIntegral(phi,mu,sza,self):
    #print phi
    #print mu
    #print sza
    #print '========'
    vza = numpy.arccos(mu)
    raa = self.rtod(phi)
    self.setAngleInfo(vza,sza,raa)
    self.LiKernel()
    return mu * self.Li[0] / numpy.pi


def readASCII(inputFile,dobands=False):
    FILE = open(inputFile,'r')
    header = FILE.readline()
    nBands = int(header.split()[2])
    bands = header.split()[3:3+nBands]
    Bands = numpy.zeros(nBands)
    for i in range(nBands):
        Bands[i] = float(bands[i])
    strdata = FILE.readlines()
    FILE.close()
    N = len(strdata)
    DOY = numpy.zeros(N)
    FLAG = numpy.zeros(N)
    VZA = numpy.zeros(N)
    SZA = numpy.zeros(N)
    RAA = numpy.zeros(N)
    REFL =  numpy.zeros([nBands,N])
    for i in range(N):
        s = strdata[i].split()
        DOY[i] = float(s[0])
        FLAG[i] = int(s[1])
        VZA[i] = float(s[2])
        SZA[i] = float(s[4])
        RAA[i] = float(s[3]) -  float(s[5])
        for j in range(nBands):
                REFL[j,i] = float(s[j+6])
    w = numpy.where(FLAG == 1)
    doy = DOY[w]
    vza = VZA[w]
    sza = SZA[w]
    raa = RAA[w]
    refl = REFL[:,w]
    if dobands == True:
        return vza,sza,raa,refl,doy,Bands
    else:
        return vza,sza,raa,refl,doy

def readPOLDER(inputFile,type=1):
    FILE = open(inputFile,'r')
    strdata = FILE.readlines()
    FILE.close()
    N = len(strdata)
    VZA = numpy.zeros(N)
    SZA = numpy.zeros(N)
    RAA = numpy.zeros(N)
    REFL =  numpy.zeros([5,N])
    for i in range(N):
        s = strdata[i].split()
        if( type == 1):
            VZA[i] = float(s[4])
            SZA[i] = float(s[2])
            RAA[i] = float(s[5])
            for j in range(5):
                REFL[j,i] = float(s[j+6])
        else:
           if (type == 2):
                VZA[i] = float(s[2])
                SZA[i] = float(s[4])
                RAA[i] = float(s[5]) - float(s[3])
                for j in range(5):
                    REFL[j,i] = float(s[j+6])
    return VZA,SZA,RAA,REFL

'''
def legend(*args, **kwargs):
    """
    Overwrites the pylab legend function.

    It adds another location identfier 'outer right'
    which locates the legend on the right side of the plot

    The args and kwargs are forwarded to the pylab legend function

    from http://www.mail-archive.com/matplotlib-users@lists.sourceforge.net/msg04256.html
    """
    #import pylab
    if kwargs.has_key('loc'):
        loc = kwargs['loc']
        loc = loc.split()

        if loc[0] == 'outer':
            # make a legend with out the location
            # remove the location setting from the kwargs
            kwargs.pop('loc')
            leg = pylab.legend(loc=(0,0), *args, **kwargs)
            frame = leg.get_frame()
            currentAxes = pylab.gca()
            currentAxesPos = numpy.array(currentAxes.get_position()).flatten()

            # scale plot by the part which is taken by the legend
            plotScaling = frame.get_width()/currentAxesPos[2]

            if loc[1] == 'right':
                # scale the plot
                currentAxes.set_position((currentAxesPos[0],currentAxesPos[1]-0.05,
                                          currentAxesPos[2] *(1-plotScaling),
                                          currentAxesPos[3]-0.05))
                # set x and y coordinates of legend
                #leg._loc = (1 + leg.axespad, 1 - frame.get_height())
                leg._loc = (1 + leg.axespad, 0)
            # doesn't work
            #if loc[1] == 'left':
            #    # scale the plot
            #    currentAxes.set_position((currentAxesPos[0] +frame.get_width(),
            #                              currentAxesPos[1],
            #                              currentAxesPos[2] *(1-plotScaling),
            #                              currentAxesPos[3]))
            #    # set x and y coordinates of legend
            #    leg._loc = (1 -.05 -  leg.axespad - frame.get_width(), 1 -frame.get_height())

            pylab.draw_if_interactive()
            return leg

    return pylab.legend(*args, **kwargs)
'''


def lutInvertRossHS(VZA,SZA,RAA,REFL,N=1000,fixXi=False,RossType='Thick',LiType='Dense',normalise=1,RecipFlag=True,MODISSPARSE=True):
    if ( fixXi != False ):
        N = 1
        rhs = numpy.array([fixXi])
    else:
        rhs = numpy.array(range(N))*10*(numpy.pi/180.)/N
    rmse = numpy.zeros(N)
    for i in range(N):
      rmse[i],P,FWD,phaseAngle = invertData(VZA,SZA,RAA,REFL,RossType=RossType,LiType=LiType,RossHS=rhs[i],normalise=normalise,RecipFlag=RecipFlag,MODISSPARSE=MODISSPARSE)
    i = numpy.argmin(rmse)
    RMSE,P,FWD,phaseAngle = invertData(VZA,SZA,RAA,REFL,RossType=RossType,LiType=LiType,RossHS=rhs[i],normalise=normalise,RecipFlag=RecipFlag,MODISSPARSE=MODISSPARSE)
    return RMSE,rhs[i],P,numpy.array(FWD),rhs,rmse,phaseAngle 

def testLisa(inputFile,buff=30,LiType='Sparse',RossType='Thick',plot=False,verbose=False,fsza=0.0,forcedoy=False):
    bu = [0.004, 0.015, 0.003, 0.004, 0.013, 0.010, 0.006]
    vza,sza,raa,refl,doy,bands = readASCII(inputFile,dobands=True)

    if type(fsza) == type(True) and  fsza == True:
        msza = numpy.median(sza)
    else:
        msza = fsza
    if verbose == True:
        print('nbar at',msza)
        
    nbands = len(bands)
    if nbands == 4:
        bux = [bu[1], bu[4], bu[0], bu[6]]
    else:
        bux = bu
    mind = min(doy)
    maxd = max(doy)
    w1 = numpy.where(doy >= (mind + buff))
    w2 = numpy.where(doy[w1] <= (maxd - buff))
    sampledays = doy[w1][w2]
    if forcedoy != False:
        sampledays = numpy.array([forcedoy])
    iso = numpy.zeros(len(bux))
    isoPost = numpy.zeros(len(bux))
    sig = numpy.zeros([len(sampledays),len(bux)])
    rel = numpy.zeros([len(sampledays),len(bux)])
    RMSE=1e20 
    count = 0
    mindoy = False 
    minrmse = False
    minP = False
    minFWD = False
    minrefl = False

    # stuff for spectral mixture model
    loff = 400.0
    lmax = 2000.0
    ll = bands-loff
    llmax = lmax-loff
    lk = ll - ll*ll/(2*llmax)
    lk = lk/max(lk)
    K = numpy.matrix(numpy.ones([3,nbands]))
    K[1][:] = lk
    mincount = -1 
    for dos in sampledays:
        rmse,P,FWD,refl,idoy,unc = lisaInvert(vza,sza,raa,refl,doy,dos,LiType=LiType,RossType=RossType,nbar=msza)
        # calculate significance of step change in 1st 2 bands
        # P[i,6] is the magnitude of step change
        dos2 = dos+1
        for i in range(len(bux)):
            iso[i] = P[i,0] + dos*P[i,3] + dos*dos* P[i,4] + dos*dos*dos*P[i,5]
            sig[count][i] = P[i,6] / (unc * bux[i])
            rel[count][i] = P[i,6]/iso[i]
            isoPost[i] = P[i,0] + dos2*P[i,3] + dos2*dos2* P[i,4] + dos2*dos2*dos2*P[i,5] + P[i,6] 
        # do spectral mixture modelling on iso
        if nbands == 7:
            # loff = 400
            # l' = l - loff
            # lmax = 2000 - loff
            # rhoBurn = a0 + a1(l' - l'^2/(2 lmax)) = a0 + a1 * lk
            # lmax = 
            # post = pre * (1-fcc) + fcc * rhoBurn      
            # post = pre * (1-fcc) + fcc * a0 + fcc * a1 * lk
            # post - pre = A + B * lk - fcc * pre
            # where
            # A = fcc * a0
            # B = fcc * a1
            y = numpy.matrix(isoPost - iso)
            K[2] = iso
            M = K * K.transpose()
            MI = M.I
            V = K*y.transpose()   
            # spectral parsamsters
            sP = numpy.array((MI*V).transpose())[0]
            fcc = -sP[2]
            a0 = sP[0]/fcc
            a1 = sP[1]/fcc
            sBurn = a0 + lk*a1 
            sFWD = iso*(1-fcc) + fcc*sBurn
            sPre = iso
            sPost = isoPost
        else:
            fcc = 0
            a0 = 0
            a1 = 0
            sBurn = 0
            sFWD = 0
            sPre = 0
            sPost = 0
        if nbands == 4:
            Test = sig[count][0] <0 and sig[count][1] < 0 and ((sig[count][2]>sig[count][0] and sig[count][2]>sig[count][1]) or (sig[count][3]>sig[count][0] and sig[count][3]>sig[count][1])) 
        else:
            Test = a0 >= 0 and a1 >= 0 and fcc >= 0 and fcc <= 1 and a0 + a1 <= 1.0 and P[1,6] < 0 and P[4,6] < 0

            # put in conditions etc...
        if Test: 
                # valid sample
                rmse1 = numpy.matrix(rmse)
                rmse1 = numpy.array(numpy.sqrt(rmse1*rmse1.transpose()/len(bux)))[0][0]
                thissig = min([sig[count][0],sig[count][1]])
                #print dos,thissig
                #if mindoy == False or thissig < minsig:
                if nbands == 4:
                    Test2 = mindoy == False or rmse1 < minrmsei1
                else:
                    Test2 = mindoy == False or fcc > maxfcc 
                if verbose:
                    print(dos,fcc,a0,a1,thissig,rmse1)
                if Test2:
                        maxpre = sPre
                        maxpost = sPost
                        maxfcc = fcc
                        maxa0 = a0
                        maxa1 = a1
                        maxsBurn = sBurn
                        maxsFWD = sFWD
                        minsig = thissig
                        mindoy = dos
                        minrmse1 = rmse1
                        minrmse = rmse
                        minP = P
                        minFWD = FWD
                        minrefl = refl
                        mincount = count
        count += 1
    if mincount != -1:
        if nbands == 4:
            return doy,minrmse,minP,minFWD,minrefl,mindoy,sig[mincount],rel[mincount]
        else:
            pass
            '''
            if plot:
                #import pylab
                x = [mindoy,mindoy]
                y = [0.0,max(numpy.array([minFWD.flatten(),minrefl.flatten()]).flatten())+0.1]
                pylab.plot(x,y)
                colours = ['k','b','g','r','c','m','y']
                for i in range(nbands):
                    norm = minP[i,0] + doy*minP[i,3] + doy*doy* minP[i,4] + doy*doy*doy*minP[i,5]
                    w = numpy.where(doy > mindoy)
                    norm[w] += minP[i,6]
                    pylab.plot(doy,minFWD[i].flatten(),colours[i] + 's',label='model '+ str(bands[i]))
                    pylab.plot(doy,minrefl[i].flatten(),colours[i] + '^',label='obs '+ str(bands[i]))
                    pylab.plot(doy,norm.flatten(),colours[i] + '-',label='norm '+ str(bands[i]))
                legend(loc='outer right')
                pylab.show()
                preday = mindoy
                postday = mindoy+1
                #print minP[:,0].shape,bands.shape
                prenorm = minP[:,0] + preday*minP[:,3] + preday*preday* minP[:,4] + preday*preday*preday*minP[:,5]
                postnorm =  minP[:,0] + postday*minP[:,3] + postday*postday* minP[:,4] + postday*postday*postday*minP[:,5] + minP[:,6]
                prenorm = numpy.squeeze(numpy.array(prenorm))
                postnorm = numpy.squeeze(numpy.array(postnorm))
                
                pylab.plot(bands,prenorm,'bo',label='pre-burn') 
                pylab.plot(bands,postnorm,'go',label='post-burn')
                pylab.plot(bands,maxsFWD,'g^',label='fwd model')
                pylab.plot(bands,maxfcc*maxsBurn,'rD',label='fcc * burn signal')
                pylab.legend(loc=0)
                pylab.show()
            '''
            return doy,minrmse,minP,minFWD,minrefl,mindoy,sig[mincount],rel[mincount],maxfcc,maxa0,maxa1
    else:
        if nbands == 4:
            return False,False,False,False,False,False,False,False
        else:
            return False,False,False,False,False,False,False,False,False,False,False

def lisaInvert(vza,sza,raa,refl,doy,dos,LiType='Sparse',RossType='Thick',xi=False,nbar=0.0):
    doy2 = doy*doy
    doy3 = doy2*doy
    kk = Kernels(vza ,sza,raa,doIntegrals=False,RossHS=xi,RossType=RossType,LiType=LiType,normalise=1,RecipFlag=True,MODISSPARSE=True,nbar=nbar) 
    K = numpy.ones([7,len(vza)])
    K[1,:] = kk.Ross[:]
    K[2,:] = kk.Li[:]
    K[3,:] = doy
    K[4,:] = doy2
    K[5,:] = doy3
    w = numpy.where(doy <= dos)
    K[6,w] = 0.0
    
    # form matrix
    K = numpy.matrix(K)
    M = K * K.transpose()
    MI = M.I
    nBands = len(refl[:,0])
    P = numpy.matrix(numpy.zeros([nBands,7]))
    for i in range(nBands):
        R = numpy.matrix(refl[i,:])
        V = K*R.transpose()
        P[i,:] = (MI*V).transpose()
    # rmse
    mse = numpy.zeros(nBands)
    FWD = refl.copy() * 0.0
    for i in range(nBands):
        FWD[i,:] = P[i,:] * K
        d = numpy.array((FWD[i,:] - refl[i,:])[0])
        mse[i] =  (d * d).mean()
    rmse =  numpy.sqrt(mse)
    return rmse,P,FWD,refl,doy, numpy.sqrt(MI[6,6])
#pylab.plot(doy,refl[0,:].flatten())
# pylab.plot(doy,FWD[0,:].flatten())
# pylab.show()
 

def testMe(fixXi=.02617993877991494365,LiType='Sparse',RossType='Thick',file='polder.modis.tiles.cover.04.dat.count.top.1.all.h08v05.256.509.dat',ofile=False,type=1,N=1000):
    VZA,SZA,RAA,REFL = readPOLDER(file,type=type)

    rmse,xi,P,FWD,x,y,phi=lutInvertRossHS(VZA,SZA,RAA,REFL,LiType=LiType,RossType=RossType,N=N,fixXi=fixXi)
    if ( ofile == True ):
        aofile = file + '.kernelModelled'
        FILE = open(aofile,'w')
        FILE.write('# xi = ' + str(xi) + ' rmse = ' + str(rmse) + '1:vza 2:sza 3:relphi 4:obs(443) 5:mod(443) 6:obs(565) 7:mod(565) 8:obs(670) 9:mod(670) 10:obs(765) 11:mod(765) 12:obs(865) 13:mod(865)\n')
        for i in range(len(VZA)):
            ostr = str(VZA[i]) + ' ' + str(SZA[i]) + ' ' + str(-RAA[i]) + ' ' 
            for j in range(5):
                ostr = ostr + str(REFL[j,i]) + ' ' + str(FWD[j,i]) + ' '
            ostr = ostr + '\n'
            FILE.write(ostr)
        FILE.close()
    vza = numpy.array(range(141))*1.0 - 70
    raa = vza*0.0 
    sza = raa - int(SZA.mean())
    sza = raa - 40.0
    kk = Kernels(vza ,sza,raa,doIntegrals=False,RossHS=xi,RossType=RossType,LiType=LiType,normalise=1,RecipFlag=True,MODISSPARSE=True)
    K = numpy.ones([3,len(vza)])
    K[1,:] = kk.Ross[:]
    K[2,:] = kk.Li[:]
    fwd = numpy.array(P * K)
    if ( ofile == True ):
        aofile = file + '.kernelPplane'
        FILE = open(aofile,'w')
        FILE.write('# pplane plot at mean sza of observations: sza = ' +  str(sza[0]) + '\n')
        FILE.write('# 1:vza(pplane -ve = hs) 2:r(443) 3:r(565) 4:r(670) 5:r(765) 6:r(865)\n')
        for i in range(len(vza)):
            ostr = str(vza[i]) + ' '      
            for j in range(5):
                ostr = ostr + str(fwd[j,i]) + ' '
            ostr = ostr + '\n'
            FILE.write(ostr)
        FILE.close()
    return P,rmse,xi
 
# w = numpy.where(SZA > 52)
# S = SZA[w]
# w1 = numpy.where(S < 55)
# s = S[w1]
# pylab.plot(phi[w][w1],REFL[0,:][w][w1],'o')
# pylab.show()

def invertData(VZA,SZA,RAA,REFL,RossType='Thick',LiType='Dense',RossHS=False,normalise=1,RecipFlag=True,MODISSPARSE=True):
    # invert    
    kk = Kernels(VZA,SZA,RAA,RossHS=RossHS,MODISSPARSE=MODISSPARSE,RecipFlag=RecipFlag,normalise=normalise,doIntegrals=False,LiType=LiType,RossType=RossType)
    K = numpy.ones([3,len(VZA)])
    K[1,:] = kk.Ross[:]
    K[2,:] = kk.Li[:]
    # form matrix
    K = numpy.matrix(K)
    M = K * K.transpose()
    MI = M.I
    nBands = len(REFL[:,0])
    P = numpy.matrix(numpy.zeros([nBands,3]))
    for i in range(nBands):
        R = numpy.matrix(REFL[i,:])
        V = K*R.transpose()
        P[i,:] = (MI*V).transpose()
    # rmse
    FWD = P * K
    d = FWD - REFL
    e = 0.0
    for i in range(nBands):
      e = e + d[i] * d[i].transpose()
    rmse = numpy.sqrt(e[0,0]/(len(VZA)*nBands))
    phaseAngle = numpy.arctan2(kk._Kernels__sinphaang,kk._Kernels__cosphaang)*180./numpy.pi
    phaseAngle = phaseAngle[0:len(VZA)]
    return rmse,P,FWD,phaseAngle

# test function
def mimic(doPrint=False,doPlot=False,RossHS=False,RecipFlag=False,thisSza=None):
    '''
    A test method to reproduce the results in Wanner et al. 1995.
    There are no parameters and a single option:
            doPrint=True    : print results to stdout (default doPrint=False)

    The method returns:
        VZA,SZA,RAA,RossThick,RossThin,LiSparse,LiDense,Roujean,LiTransit
    where all are numy arrays of dimensions 3 x nSamples 
    so:
        VZA[0,:],RossThick[0,:] are the results for sza = 0.0
        VZA[1,:],RossThick[1,:] are the results for sza = 30.0
        VZA[2,:],RossThick[2,:] are the results for sza = 60.0

    '''
    # set up the angles
    r = 89 # do results for +/- r degrees)
    if thisSza == None:
        SZAS = ma.array([0.0,-30.0,-60.0])  # sza
    else:
        SZAS = ma.array(thisSza)
    vza = numpy.array(range(2*r+1))*1.0 - r
    # set up storage info
    RossThick = ma.zeros([3,len(vza)])
    RossThin = ma.zeros([3,len(vza)])
    LiSparse = ma.zeros([3,len(vza)])
    LiDense = ma.zeros([3,len(vza)])
    Roujean = ma.zeros([3,len(vza)])
    LiTransit = ma.zeros([3,len(vza)])
    SZA = ma.zeros([3,len(vza)])
    VZA = ma.zeros([3,len(vza)])
    RAA = ma.zeros([3,len(vza)])
    # fill the angle info
    RossHS=RossHS
    for i in range(len(SZAS)):
        SZA[i,:] = SZAS[i]
        VZA[i,:] = vza[:] 
        RAA[i,:] = 0.0
        # do the kernels
        kk = Kernels(VZA[i,:] ,SZA[i,:],RAA[i,:],RossHS=RossHS,MODISSPARSE=True,RecipFlag=RecipFlag,normalise=1,doIntegrals=False,LiType='Dense',RossType='Thick')
        RossThick[i,:] = kk.Ross[:]
        LiDense[i,:] = kk.Li[:]
        if doPrint == True:
                kk.printKernels(file='RossThickLiDense.' + str(SZAS[i]) + '.dat')
                kk.printer('')
        kk = Kernels(VZA[i,:] ,SZA[i,:],RAA[i,:],RossHS=RossHS,MODISSPARSE=True,RecipFlag=RecipFlag,normalise=1,doIntegrals=False,LiType='Sparse',RossType='Thin')
        RossThin[i,:] = kk.Ross[:]
        LiSparse[i,:] = kk.Li[:]
        if doPrint == True:
                kk.printKernels(file='RossThinLiSparse.' + str(SZAS[i]) + '.dat')
                kk.printer('')
        kk = Kernels(VZA[i,:] ,SZA[i,:],RAA[i,:],RossHS=RossHS,MODISSPARSE=True,RecipFlag=RecipFlag,normalise=1,doIntegrals=False,LiType='Roujean',RossType='Thin')
        Roujean[i,:] = kk.Li[:]
        if doPrint == True:
                kk.printKernels(file='RossThinRoujean.' + str(SZAS[i]) + '.dat')
                kk.printer('')
        kk = Kernels(VZA[i,:] ,SZA[i,:],RAA[i,:],RossHS=RossHS,MODISSPARSE=True,RecipFlag=RecipFlag,normalise=1,doIntegrals=False,LiType='Transit',RossType='Thin')
        LiTransit[i,:] = kk.Li[:]
        if doPrint == True:
                kk.printKernels(file='RossThinLiTransit.' + str(SZAS[i]) + '.dat')
                kk.printer('')
    '''
    if (doPlot == True):
        #import pylab
        x = [-90.0,90.0]
        y = [0.0,0.0]
        for i in range(len(SZAS)):
            sza = SZAS[i]
            pylab.clf()
            pylab.xlabel('View Zenith Angle')
            pylab.ylabel('Kernel Value')
            pylab.title('Solar Zenith Angle ' + str(sza) + ' Degrees')
            pylab.plot(x,y)
            pylab.plot(kk.vzaDegrees,RossThick[i,:],label='RThick')
            pylab.plot(kk.vzaDegrees,RossThin[i,:],label='RThin')
            pylab.plot(kk.vzaDegrees,LiSparse[i,:],label='LiSp')
            pylab.plot(kk.vzaDegrees,LiDense[i,:],label='LiDen')
            pylab.plot(kk.vzaDegrees,Roujean[i,:],label='Roujean')
            pylab.plot(kk.vzaDegrees,LiTransit[i,:],label='LiTrans')
            pylab.axis([-90.0,90.0,-3.0,3.0])
            pylab.legend(loc=0)
            pylab.show()
    '''
    return VZA,SZA,RAA,RossThick,RossThin,LiSparse,LiDense,Roujean,LiTransit
