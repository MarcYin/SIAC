#/usr/bin/env python 
import numpy as np
import sys
sys.path.insert(0, 'util')
from scipy import ndimage, signal, optimize
from scipy.fftpack import dct, idct
from SIAC.multi_process import parmap
from SIAC.create_training_set import create_training_set
import scipy



def cloud_dilation(cloud_mask, iteration=1):
    '''
    A function for the dilation of cloud mask
   
    '''
    struct = np.ones((3,3)).astype(bool)
    dila_cloud = ndimage.binary_dilation(cloud_mask, structure=struct, iterations=iteration).astype(bool)
    return dila_cloud


class psf_optimize(object):
    def __init__(self, 
		 high_img,
		 high_indexs,
		 low_img,
                 qa,
                 cloud,
                 qa_thresh,
                 xstd = 29.75,
                 ystd = 39,
                 scale=0.95607605898444503,
                 offset=0.0086119174434039214):
       self.high_img    = high_img
       self.Hx, self.Hy = high_indexs
       self.low_img     = low_img
       self.cloud       = cloud
       self.qa_thresh   = qa_thresh
       self.qa          = qa
       self.xstd        = xstd
       self.ystd        = ystd
       self.shape       = self.high_img.shape
       self.parameters  = ['xstd', 'ystd', 'angle', 'xs', 'ys']
       self.slop        = scale
       self.off         = offset
    def _preprocess(self,):
     
        size = 2*int(round(1.96*self.ystd))# set the largest possible PSF size
        #self.high_img[0,:]=self.high_img[-1,:]=self.high_img[:,0]=self.high_img[:,-1]= -9999
        self.bad_pixs = cloud_dilation( (self.high_img < 0.0001) | self.cloud  | (self.high_img >= 1), iteration=int(size/2))
        #xstd, ystd = 29.75, 39
        #ker = self.gaussian(self.xstd, self.ystd, 0)
        gaus_2d = self.dct_gaussian(self.xstd, self.ystd, self.high_img.shape)
        if (self.xstd < 1.) or (self.ystd < 1.):
            self.conved = self.high_img
        else:
            self.conved = idct(idct(dct(dct(self.high_img, axis=0, norm = 'ortho'), axis=1, norm='ortho') * gaus_2d, \
                                                           axis=1, norm = 'ortho'), axis=0, norm='ortho')
        #self.conved = signal.fftconvolve(self.high_img, ker, mode='same')

        l_mask = (~self.low_img.mask) & (self.qa<self.qa_thresh) & (~np.isnan(self.low_img)) & (~np.isnan(self.qa))
        h_mask =  ~self.bad_pixs[self.Hx, self.Hy]
        self.lh_mask = l_mask & h_mask

    def dct_gaussian(self, xstd, ystd, shape):
        win_x, win_y = shape
        xgaus  = np.exp(-2.*(np.pi**2)*(xstd**2)*((0.5 * np.arange(win_x) /win_x)**2))
        ygaus  = np.exp(-2.*(np.pi**2)*(ystd**2)*((0.5 * np.arange(win_y) /win_y)**2))
        gaus_2d = np.outer(xgaus, ygaus)
        return gaus_2d

    def gaussian(self, xstd, ystd, angle, norm = True):
        
        win = 2*int(round(max(1.96*xstd, 1.96*ystd)))
        winx = int(round(win*(2**0.5)))
        winy = int(round(win*(2**0.5)))
        xgaus = signal.gaussian(winx, xstd)
        ygaus = signal.gaussian(winy, ystd)
        gaus  = np.outer(xgaus, ygaus)
        r_gaus = ndimage.interpolation.rotate(gaus, angle, reshape=True)
        center = np.array(r_gaus.shape)/2
        cgaus = r_gaus[int(center[0]-win/2): int(center[0]+win/2), int(center[1]-win/2):int(center[1]+win/2)]
        if norm:
            return cgaus/cgaus.sum()
        else:
            return cgaus 


    def gaus_optimize(self, p0):
        return optimize.fmin_l_bfgs_b(self.gaus_cost, p0, approx_grad=1, iprint=-1,
                                      bounds=self.bounds,maxiter=10, maxfun=10)         


    def shift_optimize(self, p0):
        return optimize.fmin(self.shift_cost, p0, full_output=1, maxiter=100, maxfun=150, disp=0)


    def gaus_cost(self, para):
        # cost for a final psf optimization
        xstd,ystd,angle, xs, ys = para 
        ker = self.gaussian(xstd,ystd,angle,True)                              
        #gaus_2d = self.gaussian(xstd, ystd, self.high_img.shape)
        conved = signal.fftconvolve(self.high_img, ker, mode='same')
        #conved = idct(idct(dct(dct(img, axis=0, norm = 'ortho'), axis=1, norm='ortho') * gaus_2d, \
        #                                axis=1, norm = 'ortho'), axis=0, norm='ortho')

        # mask bad pixels
        cos = self.cost(xs=xs, ys=ys, conved=conved)
        return cos


    def shift_cost(self, shifts):
        # cost with different shits
        xs, ys = shifts
        cos = self.cost(xs=xs, ys=ys, conved=self.conved)
        return cos


    def cost(self, xs=None, ys=None, conved = None):
        # a common cost function can be reused
        shifted_mask = np.logical_and.reduce(((self.Hx+int(xs)>=0),
                                              (self.Hx+int(xs)<self.shape[0]), 
                                              (self.Hy+int(ys)>=0),
                                              (self.Hy+int(ys)<self.shape[1])))
        mask = self.lh_mask & shifted_mask
        x_ind, y_ind = self.Hx + int(xs), self.Hy + int(ys)
        High_resolution_band, Low_resolution_band = conved[x_ind[mask], y_ind[mask]], self.low_img[mask]
        m_fed, s_fed = self.slop*Low_resolution_band+self.off, High_resolution_band
        try:
            r = scipy.stats.linregress(m_fed, s_fed)
            cost = abs(1-r.rvalue)
        except:
            cost = 100000000000.
        return cost


    def fire_shift_optimize(self,):
        #self.S2_PSF_optimization()
        self._preprocess()
        if self.lh_mask.sum() ==0:
            self.costs = np.array([100000000000.,])
            return 0,0
        min_val = [-50,-50]
        max_val = [50,50]
        ps, distributions = create_training_set([ 'xs', 'ys'], min_val, max_val, n_train=50)
        self.shift_solved = parmap(self.shift_optimize, ps)    
        self.paras, self.costs = np.array([i[0] for i in self.shift_solved]), \
                                           np.array([i[1] for i in self.shift_solved])

        if (1 - self.costs.min()) >= 0.6:
            xs, ys = self.paras[self.costs==np.nanmin(self.costs)][0].astype(int)
        else:
            xs, ys = 0, 0
        #print 'Best shift is ', xs, ys, 'with the correlation of', 1-self.costs.min()
        return xs, ys


    def fire_gaus_optimize(self,):
        xs, ys = self.fire_shift_optimize()
        if self.costs.min()<0.1:
            min_val = [4, 4, -15,xs-2,ys-2]
            max_val = [40,40, 15, xs+2,ys+2]
            self.bounds = [4,40],[4,40],[-15,15],[xs-2,xs+2],[ys-2, ys+2]

            ps, distributions = create_training_set(self.parameters, min_val, max_val, n_train=50)
            print('Start solving...')
            self.gaus_solved = parmap(self.gaus_optimize, ps, nprocs=5)
            result = np.array([np.hstack((i[0], i[1])) for i in  self.gaus_solved])
            print('solved psf', dict(zip(self.parameters+['cost',],result[np.argmin(result[:,-1])])))
            return result[np.argmin(result[:,-1]),:]
        else:
            print('Cost is too large, plese check!')
            return []
