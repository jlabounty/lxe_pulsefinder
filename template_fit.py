import numpy as np
import os 
import sys
import numba 
import matplotlib.pyplot as plt
import iminuit 
from iminuit.cost import LeastSquares
from iminuit import Minuit
from scipy.interpolate import CubicSpline, UnivariateSpline, interp1d

class LXeTemplate():
    def __init__(self) -> None:
        pass

    def form_template(self, xs, ys):
        
        #shift such that the peak of the template lies at x = 0
        these_ys = np.array(ys)
        these_xs = np.array(xs) - xs[ np.where(these_ys == np.amax(these_ys))[0] ][0]

        # self.template = CubicSpline(these_xs, these_ys, 
        #             axis=0, bc_type='not-a-knot', 
        #             extrapolate=True    )

        # self.template = UnivariateSpline(these_xs, these_ys)
        self.template = interp1d(these_xs, these_ys, fill_value=0, bounds_error=False)

    def save(self, outfile):
        import pickle 
        with open(outfile, 'wb') as fout:
            pickle.dump(self, fout)

    def load(infile):
        import pickle
        with open(infile, 'rb') as fin:
            return pickle.load(fin)


class TemplateFit():
    def __init__(self, data:np.ndarray, template:LXeTemplate or str, chi2limit=2, scalex=True, scaley=True, verbose=True) -> None:
        self.data = data 
        if(type(template) is str):
            self.template = LXeTemplate.load(template).template
        elif(type(template) is LXeTemplate):
            self.template = template.template
        else:
            self.template = template
            # raise NotImplementedError
        self.current_guess = [0,0,0]
        self.chi2 = -1
        self.verbose=verbose

        self.chi2limit = chi2limit      # limit below which we say a fit has converged
        self.scalex = scalex            # whether we can scale the template amplitude to make the fit converge better
        self.scaley = scaley            # whether we can stretch/compress the template to make the fit better.
        self.convergencelimit = 10      # number of iterations to try before giving up

    def template_function(self, x:np.ndarray, p:np.ndarray):
        '''
            Uses self.template to form a function which can be used by the iMinuit least squares fitter
            Takes in a variable number of parameters, which will be used to determine how many time/amplitude 
            shifted versions of the template to apply.

            Returns: combined template evaluated at the parameters
        '''

        nparams = 3
        number_of_fits = int(len(p) / nparams)
        assert len(p) % nparams == 0, ValueError()

        ys = np.zeros_like(x)
        for i in range(number_of_fits):
            amplitude_factor = p[i*nparams + 0]
            stretch_factor   = p[i*nparams + 1]
            time_offset      = p[i*nparams + 2]

            ys += amplitude_factor * self.template( x + time_offset )

        return ys

    def do_single_fit(self,xs,ys,param_guess):
        yerrs = np.ones_like(xs)
        minimizer = LeastSquares( xs, ys, yerrs, self.template_function  )
        m = Minuit(minimizer, param_guess )  # starting values for minimization
        # m.limits[0]

        nparams = 3
        number_of_fits = int(len(param_guess) / nparams)
        assert len(param_guess) % nparams == 0, ValueError()

        for i in range(number_of_fits):
            m.limits[i*nparams    ] = (0, None)
            m.fixed[i*nparams + 1 ] = True
            # m.limits[i*nparams + 2] = (0.001, None)

        m.migrad(2)  # finds minimum of least_squares function
        m.hesse()   # accurately computes uncertainties

        # print(m)
        return m

    def identify_local_maxima(self, data):
        '''
            Function to pull out the initial conditions for a LXe pulse fit. Identifies local maxima and tags them.
        '''
        maxima = []
        times = []

        from scipy.signal import find_peaks
        # for local maxima
        peaks, _ = find_peaks(data[1], height=100, threshold=10, width=5)
        for x in peaks:
            maxima.append(data[1][x])
            times.append(data[0][x])

        return maxima, times

    def do_fit(self):
        self.chi2 = 1e12
        self.niterations = 0

        self.intermediate_xs = self.data[0]
        self.intermediate_ys = self.data[1]


        maxima, times = self.identify_local_maxima(self.data)

        if(len(maxima) < 2):
            # self.current_guess = [1,0,0]
            new_maximum = np.amax( self.intermediate_ys )
            maximum_time = self.intermediate_xs[ np.where( self.intermediate_ys  == new_maximum )[0] ][0]
            self.current_guess = [new_maximum,0, -1.0*maximum_time]
        else:
            # raise NotImplementedError
            self.current_guess = []
            for time, maxi in zip(times,maxima):
                self.current_guess += [maxi, 0, -time]

        while(True):
            if(self.verbose):
                print("**********************************************************************************")
                print("Fitting with Current guess:", self.current_guess)
                # print(new_maximum, maximum_time)
            self.niterations += 1
            if(self.niterations > self.convergencelimit):
                print(f"Warning: Fit did not converge in {self.convergencelimit} iterations")
                break
                # raise TimeoutError(f"Error: Fit did not converge in {self.convergencelimit} iterations")

            m = self.do_single_fit(self.intermediate_xs, self.intermediate_ys, self.current_guess)
            if(self.verbose):
                print('Fit valid:', m.valid)
                print('Fit params:', m.values)
                print("Fit chi2:", m.fmin.reduced_chi2)
            # input("so?")

            self.current_guess = list(m.values)
            self.chi2 = m.fmin.reduced_chi2

            # if(self.chi2 < self.chi2limit):
            #     break

            # try adding another pulse to the mix, and see if that reduces the chi2
            # print( self.intermediate_ys , self.template_function(self.intermediate_xs, self.current_guess) )
            if(self.verbose):
                plt.plot(self.intermediate_ys , label="Data")
                plt.plot( self.template_function(self.intermediate_xs, self.current_guess), label="Fit (before additional)")
                plt.legend()
                plt.show()
            residuals = self.intermediate_ys - self.template_function(self.intermediate_xs, self.current_guess)
            new_maximum = np.amax( residuals )
            maximum_time = self.intermediate_xs[ np.where( residuals == new_maximum )[0] ][0]
            if(self.verbose):
                print(new_maximum, maximum_time)

            m2 = self.do_single_fit(self.intermediate_xs, residuals, [new_maximum, 0, -1.0*maximum_time])
            if(self.verbose):
                print('   -> Residual Fit valid:', m2.valid)
                print('   -> Residual Fit params:', m2.values)
                print("   -> Residual Fit chi2:", m2.fmin.reduced_chi2)
                # input("so?")

                plt.plot(self.intermediate_xs, self.intermediate_ys, label="Data")
                plt.plot(self.intermediate_xs, residuals, label="Residuals")
                plt.plot(self.intermediate_xs, self.template_function(self.intermediate_xs, [new_maximum,0, -1.0*maximum_time]), label="New Guess")
                plt.plot(self.intermediate_xs, self.template_function(self.intermediate_xs, list(m2.values)), label="New Fit")
                plt.legend()
                plt.show()

            if(m2.fmin.reduced_chi2 < m.fmin.reduced_chi2):
                self.current_guess += list(m2.values)
            else:
                break
                # pass
            
            # break

    def plot(self):
        '''
            Plots the current fit iteration
        '''
        fig, ax = plt.subplots(figsize=(15,5))
        plt.plot(*self.data, label='Data')
        xs = np.linspace(np.amin(self.data[0]), np.amax(self.data[0]), 1000)

        param_string = ''
        ylim = ax.get_ylim()
        for i in range(int(len(self.current_guess)/3)):
            if(self.verbose):
                print(i, self.current_guess[i*3:i*3+3])
            param_string += f'\nTemplate {i} with Parameters: {[round(x,3) for x in self.current_guess[i*3:i*3+3]]}'
            plt.plot([-self.current_guess[i*3+2],-self.current_guess[i*3+2]], ylim, "r:")
            plt.plot(xs, self.template_function(xs, self.current_guess[i*3:i*3+3]), color="xkcd:light grey")

        param_string += f"\nchi2 = {self.chi2}"

        plt.plot(xs, self.template_function(xs, self.current_guess), label="Fit"+param_string)

        ax.set_ylim(*ylim)

        plt.grid()
        plt.legend()
        plt.tight_layout()

        return fig,ax 

            
def make_fake_data_from_template(template, times, amplitudes, noise=True, noisefloor=10):
    '''
        Helper function to generate fake data from a pulse shape template
    '''
    assert len(times) == len(amplitudes)

    limits = (np.amin(times) - 10, np.amax(times) + 400)
    xs = np.linspace(*limits, int(limits[1] - limits[0]))
    ys = np.zeros_like(xs)

    for i, (t, a) in enumerate(zip(times,amplitudes)):
        # print(i, t, a)
        ys += template(xs - t) * a 
        # plt.plot(xs, template(xs - t) * a )

    # plt.show()

    if noise:
        ys += np.random.normal(0, noisefloor, ys.size)

    return xs, ys