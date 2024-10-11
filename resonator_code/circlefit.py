'''
This is a backend module for fitting complex resonator data.
If you're trying to fit a resonator, best to use the class Resonator defined in resonator.py by
from fitTools.resonator import Resonator.
'''
import numpy as np
import scipy.optimize as spopt
import scipy.linalg as lnalg ### added
from scipy import stats


class circlefit(object):
    '''
    contains all the circlefit procedures
    see http://scitation.aip.org/content/aip/journal/rsi/86/2/10.1063/1.4907935
    arxiv version: http://arxiv.org/abs/1410.3365
    '''
    def _remove_cable_delay(self,f_data,z_data, delay):
        return z_data/np.exp(2j*np.pi*f_data*delay)

    def _center(self,z_data,zc):
        return z_data-zc

    def _dist(self,x):
        np.absolute(x,x)
        c = (x > np.pi).astype(int)
        return x+c*(-2.*x+2.*np.pi)

    def _periodic_boundary(self,x,bound):
        return np.fmod(x,bound)-np.trunc(x/bound)*bound

    def _phase_fit_wslope(self,f_data,z_data,theta0, Ql, fr, slope):
        phase = np.angle(z_data)
        def residuals(p,x,y):
            theta0, Ql, fr, slope = p
            err = self._dist(y - (theta0+2.*np.arctan(2.*Ql*(1.-x/fr))-slope*x))
            return err
        p0 = [theta0, Ql, fr, slope]
        p_final = spopt.leastsq(residuals,p0,args=(np.array(f_data),np.array(phase)),maxfev=20000)
        return p_final[0]

    def _phase_fit(self,f_data,z_data,theta0, Ql, fr):
        phase = np.angle(z_data)
        def residuals_1(p,x,y,Ql):
            theta0, fr = p
            err = self._dist(y - (theta0+2.*np.arctan(2.*Ql*(1.-x/fr))))
            return err
        def residuals_2(p,x,y,theta0):
            Ql, fr = p
            err = self._dist(y - (theta0+2.*np.arctan(2.*Ql*(1.-x/fr))))
            return err
        def residuals_3(p,x,y,theta0,Ql):
            fr = p
            err = self._dist(y - (theta0+2.*np.arctan(2.*Ql*(1.-x/fr))))
            return err
        def residuals_4(p,x,y,theta0,fr):
            Ql = p
            err = self._dist(y - (theta0+2.*np.arctan(2.*Ql*(1.-x/fr))))
            return err
        def residuals_5(p,x,y):
            theta0, Ql, fr = p
            err = self._dist(y - (theta0+2.*np.arctan(2.*Ql*(1.-x/fr))))
            return err
        p0 = [theta0, fr]
        p_final = spopt.leastsq(lambda a,b,c: residuals_1(a,b,c,Ql),p0,args=(f_data,phase),maxfev=20000)#,ftol=1e-12,xtol=1e-12)
        theta0, fr = p_final[0]
        p0 = [Ql, fr]
        p_final = spopt.leastsq(lambda a,b,c: residuals_2(a,b,c,theta0),p0,args=(f_data,phase),maxfev=20000)#,ftol=1e-12,xtol=1e-12)
        Ql, fr = p_final[0]
        p0 = fr
        p_final = spopt.leastsq(lambda a,b,c: residuals_3(a,b,c,theta0,Ql),p0,args=(f_data,phase),maxfev=20000)#,ftol=1e-12,xtol=1e-12)
        fr = p_final[0][0]
        p0 = Ql
        p_final = spopt.leastsq(lambda a,b,c: residuals_4(a,b,c,theta0,fr),p0,args=(f_data,phase),maxfev=20000)#,ftol=1e-12,xtol=1e-12)
        Ql = p_final[0][0]
        p0 = np.array([theta0, Ql, fr], dtype='float64')
        p_final = spopt.leastsq(residuals_5,p0,args=(f_data,phase),maxfev=20000)
        return p_final[0]

    def _fit_skewed_lorentzian(self,f_data,z_data):
        amplitude = np.absolute(z_data)
        amplitude_sqr = amplitude**2
        A1a = np.minimum(amplitude_sqr[0],amplitude_sqr[-1])
        A3a = -np.max(amplitude_sqr)
        fra = f_data[np.argmin(amplitude_sqr)]
        def residuals(p,x,y):
            A2, A4, Ql = p
            err = y -(A1a+A2*(x-fra)+(A3a+A4*(x-fra))/(1.+4.*Ql**2*((x-fra)/fra)**2))
            return err
        p0 = [0., 0., 1e3]
        p_final = spopt.leastsq(residuals,p0,args=(np.array(f_data),np.array(amplitude_sqr)),maxfev=20000)
        A2a, A4a, Qla = p_final[0]

        def residuals2(p,x,y):
            A1, A2, A3, A4, fr, Ql = p
            err = y -(A1+A2*(x-fr)+(A3+A4*(x-fr))/(1.+4.*Ql**2*((x-fr)/fr)**2))
            return err
        def fitfunc(x,A1, A2, A3, A4, fr, Ql):
            return A1+A2*(x-fr)+(A3+A4*(x-fr))/(1.+4.*Ql**2*((x-fr)/fr)**2)

        p0 = [A1a, A2a , A3a, A4a, fra, Qla]

        #p_final = spopt.leastsq(residuals2,p0,args=(np.array(f_data),np.array(amplitude_sqr)))
        try:
            # popt, pcov = spopt.curve_fit(fitfunc, np.array(f_data), np.array(amplitude_sqr),method='trf',bounds=[lb,ub],p0=p0,maxfev=20000)
            popt, pcov = spopt.curve_fit(fitfunc, np.array(f_data), np.array(amplitude_sqr),p0=p0,maxfev=20000)
        #A1, A2, A3, A4, fr, Ql = p_final[0]
        #print(p_final[0][5])
            if pcov is not None:
                self.df_error = np.sqrt(pcov[4][4])
                self.dQl_error = np.sqrt(pcov[5][5])
            else:
                self.df_error = np.inf
                self.dQl_error = np.inf
        except:
            popt = p0
            self.df_error = np.inf
            self.dQl_error = np.inf
        #return p_final[0]
        return popt

    def _fit_circle(self,z_data, refine_results=False):
        def calc_moments(z_data):
            xi = z_data.real
            xi_sqr = xi*xi
            yi = z_data.imag
            yi_sqr = yi*yi
            zi = xi_sqr+yi_sqr
            Nd = float(len(xi))
            xi_sum = xi.sum()
            yi_sum = yi.sum()
            zi_sum = zi.sum()
            xiyi_sum = (xi*yi).sum()
            xizi_sum = (xi*zi).sum()
            yizi_sum = (yi*zi).sum()
            return np.array([ [(zi*zi).sum(), xizi_sum, yizi_sum, zi_sum],  \
            [xizi_sum, xi_sqr.sum(), xiyi_sum, xi_sum], \
            [yizi_sum, xiyi_sum, yi_sqr.sum(), yi_sum], \
            [zi_sum, xi_sum, yi_sum, Nd] ])

        M = calc_moments(z_data)
        ###added
        B = np.array([[0,0,0,-2],[0,1,0,0,],[0,0,1,0],[-2,0,0,0]])

        # the code commented out below was messy and has been replaced with scipy linear algebra functions
        '''
        a0 = (
            ((M[2][0]*M[3][2]-M[2][2]*M[3][0])*M[1][1]-M[1][2]*M[2][0]*M[3][1]-M[1][0]*M[2][1]*M[3][2]+M[1][0]*M[2][2]*M[3][1]+M[1][2]*M[2][1]*M[3][0])*M[0][3]
            +(M[0][2]*M[2][3]*M[3][0]-M[0][2]*M[2][0]*M[3][3]+M[0][0]*M[2][2]*M[3][3]-M[0][0]*M[2][3]*M[3][2])*M[1][1]
            +(M[0][1]*M[1][3]*M[3][0]-M[0][1]*M[1][0]*M[3][3]-M[0][0]*M[1][3]*M[3][1])*M[2][2]
            +(-M[0][1]*M[1][2]*M[2][3]-M[0][2]*M[1][3]*M[2][1])*M[3][0]
            +((M[2][3]*M[3][1]-M[2][1]*M[3][3])*M[1][2]+M[2][1]*M[3][2]*M[1][3])*M[0][0]
            +(M[1][0]*M[2][3]*M[3][2]+M[2][0]*(M[1][2]*M[3][3]-M[1][3]*M[3][2]))*M[0][1]
            +((M[2][1]*M[3][3]-M[2][3]*M[3][1])*M[1][0]+M[1][3]*M[2][0]*M[3][1])*M[0][2]
        )
        a1 = (
            ((M[3][0]-2.*M[2][2])*M[1][1]-M[1][0]*M[3][1]+M[2][2]*M[3][0]+2.*M[1][2]*M[2][1]-M[2][0]*M[3][2])*M[0][3]
            +(2.*M[2][0]*M[3][2]-M[0][0]*M[3][3]-2.*M[2][2]*M[3][0]+2.*M[0][2]*M[2][3])*M[1][1]
            +(-M[0][0]*M[3][3]+2.*M[0][1]*M[1][3]+2.*M[1][0]*M[3][1])*M[2][2]
            +(-M[0][1]*M[1][3]+2.*M[1][2]*M[2][1]-M[0][2]*M[2][3])*M[3][0]
            +(M[1][3]*M[3][1]+M[2][3]*M[3][2])*M[0][0]
            +(M[1][0]*M[3][3]-2.*M[1][2]*M[2][3])*M[0][1]
            +(M[2][0]*M[3][3]-2.*M[1][3]*M[2][1])*M[0][2]
            -2.*M[1][2]*M[2][0]*M[3][1]-2.*M[1][0]*M[2][1]*M[3][2]
        )
        a2 = (
            (2.*M[1][1]-M[3][0]+2.*M[2][2])*M[0][3]
            +(2.*M[3][0]-4.*M[2][2])*M[1][1]
            -2.*M[2][0]*M[3][2]+2.*M[2][2]*M[3][0]+M[0][0]*M[3][3]+4.*M[1][2]*M[2][1]-2.*M[0][1]*M[1][3]-2.*M[1][0]*M[3][1]-2.*M[0][2]*M[2][3]
        )
        a3 = (-2.*M[3][0]+4.*M[1][1]+4.*M[2][2]-2.*M[0][3])
        a4 = -4.

        def func(x):
            return a0+a1*x+a2*x*x+a3*x*x*x+a4*x*x*x*x

        def d_func(x):
            return a1+2*a2*x+3*a3*x*x+4*a4*x*x*x

        x0 = spopt.fsolve(func, 0., fprime=d_func)

        def solve_eq_sys(val,M):
            #prepare
            M[3][0] = M[3][0]+2*val
            M[0][3] = M[0][3]+2*val
            M[1][1] = M[1][1]-val
            M[2][2] = M[2][2]-val
            return np.linalg.svd(M)

        U,s,Vt = solve_eq_sys(x0[0],M)
        '''

        ### added start
        # find the eigensystem of the generalized eigenvalue problem Ma = lambda Ba, i.e., find roots of det(M - lambda B)
        spectrum, vectors = lnalg.eig(M,b=B)

        # the smallest positive eigenvalue's eigenvector, A = (A,B,C,D)'
        # minimizes F = A'MA and provides the parameters needed for fit
        # note that the ' signifies transpose above.
        m = min(i for i in spectrum if i >= 0)
        vec_index = list(spectrum).index(m)

        #old code
        #A_vec = Vt[np.argmin(s),:]

        # select the correct eigenvector (which is normalized by default)
        A_vec_normal = vectors[:,vec_index]
        # scalar multiplication of normalized eigenvector to satisfy constraint B^2 + C^2 - 4AC = 1
        A_vec = A_vec_normal/np.sqrt(A_vec_normal[1]**2 + A_vec_normal[2]**2 - 4.*A_vec_normal[0]*A_vec_normal[3])
        ### added stop

        # the center of circle is found as (xc,yc) = (-B/2A,-C/2A)
        xc = -A_vec[1]/(2.*A_vec[0])
        yc = -A_vec[2]/(2.*A_vec[0])

        # find the radius of the circle as r0 = |1/2A|
        r0 = 1./(2.*np.absolute(A_vec[0]))
        #print('norm = {}, sqrt term = {}'.format(np.linalg.norm(A_vec),np.sqrt(
        #   A_vec[1]*A_vec[1]+A_vec[2]*A_vec[2]-4.*A_vec[0]*A_vec[3])))
        if refine_results:
            print("agebraic r0: " + str(r0))
            xc,yc,r0 = self._fit_circle_iter(z_data, xc, yc, r0)
            r0 = self._fit_circle_iter_radialweight(z_data, xc, yc, r0)
            print("iterative r0: " + str(r0))

        # return the center position and radius of the circle
        return xc, yc, r0

    def _guess_delay(self,f_data,z_data):
        phase2 = np.unwrap(np.angle(z_data))
        # fit the first and last 10% of the phase curve.
        window = int(len(f_data)*0.1)
        gradient1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(f_data[:window],phase2[:window])
        gradient2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(f_data[-window:],phase2[-window:])
        gradient = np.mean((gradient1,gradient2))
        return gradient*(-1.)/(np.pi*2.)


    def _fit_delay(self,f_data,z_data,delay=0.,maxiter=0):
        def residuals(p,x,y):
            phasedelay = p
            z_data_temp = y*np.exp(1j*(2.*np.pi*phasedelay*x))
            xc,yc,r0 = self._fit_circle(z_data_temp)
            err = np.sqrt((z_data_temp.real-xc)**2+(z_data_temp.imag-yc)**2)-r0
            return err
        p_final = spopt.leastsq(residuals,delay,args=(f_data,z_data),maxfev=maxiter,ftol=1e-12,xtol=1e-12)
        return p_final[0][0]

    def _fit_delay_alt_bigdata(self,f_data,z_data,delay=0.,maxiter=0):
        def residuals(p,x,y):
            phasedelay = p
            z_data_temp = 1j*2.*np.pi*phasedelay*x
            np.exp(z_data_temp,out=z_data_temp)
            np.multiply(y,z_data_temp,out=z_data_temp)
            #z_data_temp = y*np.exp(1j*(2.*np.pi*phasedelay*x))
            xc,yc,r0 = self._fit_circle(z_data_temp)
            err = np.sqrt((z_data_temp.real-xc)**2+(z_data_temp.imag-yc)**2)-r0
            return err
        p_final = spopt.leastsq(residuals,delay,args=(f_data,z_data),maxfev=maxiter,ftol=1e-12,xtol=1e-12)
        return p_final[0][0]

    def _fit_entire_model(self,f_data,z_data,fr,absQc,Ql,phi0,delay,a=1.,alpha=0.,maxiter=20000):
        '''
        fits the whole model: a*exp(i*alpha)*exp(-2*pi*i*f*delay) * [ 1 - {Ql/Qc*exp(i*phi0)} / {1+2*i*Ql*(f-fr)/fr} ]
        '''
        def funcsqr(p,x):
            fr,absQc,Ql,phi0,delay,a,alpha = p
            return np.array([np.absolute( ( a*np.exp(np.complex(0,alpha))*np.exp(np.complex(0,-2.*np.pi*delay*x[i])) * ( 1 - (Ql/absQc*np.exp(np.complex(0,phi0)))/(np.complex(1,2*Ql*(x[i]-fr)/fr)) )  ) )**2 for i in range(len(x))])
        def residuals(p,x,y):
            fr,absQc,Ql,phi0,delay,a,alpha = p
            err = [np.absolute( y[i] - ( a*np.exp(np.complex(0,alpha))*np.exp(np.complex(0,-2.*np.pi*delay*x[i])) * ( 1 - (Ql/absQc*np.exp(np.complex(0,phi0)))/(np.complex(1,2*Ql*(x[i]-fr)/fr)) )  ) ) for i in range(len(x))]
            return err
        p0 = [fr,absQc,Ql,phi0,delay,a,alpha]
        (popt, params_cov, infodict, errmsg, ier) = spopt.leastsq(residuals,p0,args=(np.array(f_data),np.array(z_data)),full_output=True,maxfev=20000)
        len_ydata = len(np.array(f_data))
        if (len_ydata > len(p0)) and params_cov is not None:  #p_final[1] is cov_x data  #this caculation is from scipy curve_fit routine - no idea if this works correctly...
            s_sq = (funcsqr(popt, np.array(f_data))).sum()/(len_ydata-len(p0))
            params_cov = params_cov * s_sq
        else:
            params_cov = np.inf
        return popt, params_cov, infodict, errmsg, ier

    #

    def _optimizedelay(self,f_data,z_data,Ql,fr,maxiter=20000):
        xc,yc,r0 = self._fit_circle(z_data)
        z_data = self._center(z_data,np.complex(xc,yc))
        theta, Ql, fr, slope = self._phase_fit_wslope(f_data,z_data,0.,Ql,fr,0.)
        delay = 0.
        for i in range(maxiter-1): #interate to get besser phase delay term
            delay = delay - slope/(2.*2.*np.pi)
            z_data_corr = self._remove_cable_delay(f_data,z_data,delay)
            xc, yc, r0 = self._fit_circle(z_data_corr)
            z_data_corr2 = self._center(z_data_corr,np.complex(xc,yc))
            theta0, Ql, fr, slope = self._phase_fit_wslope(f_data,z_data_corr2,0.,Ql,fr,0.)
        delay = delay - slope/(2.*2.*np.pi)  #start final interation
        return delay

    def _fit_circle_iter(self,z_data, xc, yc, rc):
        '''
        this is the radial weighting procedure
        it improves your fitting value for the radius = Ql/Qc
        use this to improve your fit in presence of heavy noise
        after having used the standard algebraic fir_circle() function
        the weight here is: W=1/sqrt((xc-xi)^2+(yc-yi)^2)
        this works, because the center of the circle is usually much less
        corrupted by noise than the radius
        '''
        xdat = z_data.real
        ydat = z_data.imag
        def fitfunc(x,y,xc,yc):
            return np.sqrt((x-xc)**2+(y-yc)**2)
        def residuals(p,x,y):
            xc,yc,r = p
            temp = (r-fitfunc(x,y,xc,yc))
            return temp
        p0 = [xc,yc,rc]
        p_final = spopt.leastsq(residuals,p0,args=(xdat,ydat),maxfev=20000)
        xc,yc,rc = p_final[0]
        return xc,yc,rc

    def _fit_circle_iter_radialweight(self,z_data, xc, yc, rc):
        '''
        this is the radial weighting procedure
        it improves your fitting value for the radius = Ql/Qc
        use this to improve your fit in presence of heavy noise
        after having used the standard algebraic fir_circle() function
        the weight here is: W=1/sqrt((xc-xi)^2+(yc-yi)^2)
        this works, because the center of the circle is usually much less
        corrupted by noise than the radius
        '''
        xdat = z_data.real
        ydat = z_data.imag
        def fitfunc(x,y):
            return np.sqrt((x-xc)**2+(y-yc)**2)
        def weight(x,y):
            try:
                res = 1./np.sqrt((xc-x)**2+(yc-y)**2)
            except:
                res = 1.
            return res
        def residuals(p,x,y):
            r = p[0]
            temp = (r-fitfunc(x,y))*weight(x,y)
            return temp
        p0 = [rc]
        p_final = spopt.leastsq(residuals,p0,args=(xdat,ydat),maxfev=20000)
        return p_final[0][0]

    def _get_errors(self,residual,xdata,ydata,fitparams):
        '''
        wrapper for get_cov, only gives the errors and chisquare
        '''
        chisqr, cov = self._get_cov(residual,xdata,ydata,fitparams)
        if cov is not None:
            errors = np.sqrt(np.diagonal(cov))
        else:
            errors = None
        return chisqr, errors

    def _residuals_notch_full(self,p,x,y):
        fr,absQc,Ql,phi0,delay,a,alpha = p
        err = np.absolute( y - ( a*np.exp(np.complex(0,alpha))*np.exp(np.complex(0,-2.*np.pi*delay*x)) * ( 1 - (Ql/absQc*np.exp(np.complex(0,phi0)))/(np.complex(1,2*Ql*(x-fr)/float(fr))) )  ) )
        return err

    def _residuals_notch_ideal(self,p,x,y):
        fr,absQc,Ql,phi0 = p
        #if fr == 0: print(p)
        err = np.absolute( y - (  ( 1. - (Ql/float(absQc)*np.exp(1j*phi0))/(1+2j*Ql*(x-fr)/float(fr)) )  ) )
        #if np.isinf((np.complex(1,2*Ql*(x-fr)/float(fr))).imag):
         #   print(np.complex(1,2*Ql*(x-fr)/float(fr)))
          #  print("x: " + str(x))
           # print("Ql: " +str(Ql))
            #print("fr: " +str(fr))
        return err

    def _residuals_notch_ideal_complex(self,p,x,y):
        fr,absQc,Ql,phi0 = p
        #if fr == 0: print(p)
        err = y - (  ( 1. - (Ql/float(absQc)*np.exp(1j*phi0))/(1+2j*Ql*(x-fr)/float(fr)) )  )
        #if np.isinf((np.complex(1,2*Ql*(x-fr)/float(fr))).imag):
         #   print(np.complex(1,2*Ql*(x-fr)/float(fr)))
          #  print("x: " + str(x))
           # print("Ql: " +str(Ql))
            #print("fr: " +str(fr))
        return err

    def _residuals_directrefl(self,p,x,y):
        fr,Qc,Ql = p
        #if fr == 0: print(p)
        err = y - ( 2.*Ql/Qc - 1. + 2j*Ql*(fr-x)/fr ) / ( 1. - 2j*Ql*(fr-x)/fr )
        #if np.isinf((np.complex(1,2*Ql*(x-fr)/float(fr))).imag):
         #   print(np.complex(1,2*Ql*(x-fr)/float(fr)))
          #  print("x: " + str(x))
           # print("Ql: " +str(Ql))
            #print("fr: " +str(fr))
        return err

    def _residuals_transm_ideal(self,p,x,y):
        fr,Ql = p
        err = np.absolute( y -   ( 1./(np.complex(1,2*Ql*(x-fr)/float(fr))) )   )
        return err


    def _get_cov_fast_notch(self,xdata,ydata,fitparams): #enhanced by analytical derivatives
        #derivatives of notch_ideal model with respect to parameters
        def dS21_dQl(p,f):
            fr,absQc,Ql,phi0 = p
            return - (np.exp(1j*phi0) * fr**2) / (absQc * (fr+2j*Ql*f-2j*Ql*fr)**2 )

        def dS21_dQc(p,f):
            fr,absQc,Ql,phi0 = p
            return (np.exp(1j*phi0) * Ql*fr) / (2j*(f-fr)*absQc**2*Ql+absQc**2*fr )

        def dS21_dphi0(p,f):
            fr,absQc,Ql,phi0 = p
            return - (1j*Ql*fr*np.exp(1j*phi0) ) / (2j*(f-fr)*absQc*Ql+absQc*fr )

        def dS21_dfr(p,f):
            fr,absQc,Ql,phi0 = p
            return - (2j*Ql**2*f*np.exp(1j*phi0) ) / (absQc * (fr+2j*Ql*f-2j*Ql*fr)**2 )

        u = self._residuals_notch_ideal_complex(fitparams,xdata,ydata)
        chi = np.absolute(u)
        u = u/chi  # unit vector pointing in the correct direction for the derivative

        aa = dS21_dfr(fitparams,xdata)
        bb = dS21_dQc(fitparams,xdata)
        cc = dS21_dQl(fitparams,xdata)
        dd = dS21_dphi0(fitparams,xdata)

        Jt = np.array([aa.real*u.real+aa.imag*u.imag, bb.real*u.real+bb.imag*u.imag\
                , cc.real*u.real+cc.imag*u.imag, dd.real*u.real+dd.imag*u.imag  ])
        A = np.dot(Jt,np.transpose(Jt))
        chisqr = 1./float(len(xdata)-len(fitparams)) * (chi**2).sum()
        try:
            cov = np.linalg.inv(A)*chisqr
        except:
            cov = None
        return chisqr, cov

    def _get_cov_fast_directrefl(self,xdata,ydata,fitparams): #enhanced by analytical derivatives
        #derivatives of notch_ideal model with respect to parameters
        def dS21_dQl(p,f):
            fr,Qc,Ql = p
            return 2.*fr**2/( Qc*(2j*Ql*fr-2j*Ql*f+fr)**2 )

        def dS21_dQc(p,f):
            fr,Qc,Ql = p
            return 2.*Ql*fr/(2j*Qc**2*(f-fr)*Ql-Qc**2*fr)

        def dS21_dfr(p,f):
            fr,Qc,Ql = p
            return - 4j*Ql**2*f/(Qc*(2j*Ql*fr-2j*Ql*f+fr)**2)

        u = self._residuals_directrefl(fitparams,xdata,ydata)
        chi = np.absolute(u)
        u = u/chi  # unit vector pointing in the correct direction for the derivative

        aa = dS21_dfr(fitparams,xdata)
        bb = dS21_dQc(fitparams,xdata)
        cc = dS21_dQl(fitparams,xdata)

        Jt = np.array([aa.real*u.real+aa.imag*u.imag, bb.real*u.real+bb.imag*u.imag\
                , cc.real*u.real+cc.imag*u.imag  ])
        A = np.dot(Jt,np.transpose(Jt))
        chisqr = 1./float(len(xdata)-len(fitparams)) * (chi**2).sum()
        try:
            cov = np.linalg.inv(A)*chisqr
        except:
            cov = None
        return chisqr, cov