import pandas as pd
import matplotlib.pyplot as pl
import numpy as np
from scipy import stats
import scipy.optimize as sciopt
import re
import cv2 

print('For details on the analysis protocol and calculation:\n>>> help(Analysis)\n>>> help(Sample)')


equipment = {
    'filter785TransmitanceAt800': 0.35,
    'filter900TransmitanceAt800':0.83
}


class Sample():
    '''Creates an object to handle its data and all the information 
    related to its experiment.

    To initialise a sample object use the example:
    ```
    >>> ucnp = Sample('path/to/data', expId='20210101-101010', sampleType='ucnp')
    ```
    Sample type can be: ucnp, dye, empty, or diluter. It will define the rules of calculation
    for certain methods.
    '''

    _data = pd.DataFrame()
    _rawdata = pd.DataFrame()
    _dataDeconv = pd.DataFrame({'power_density': [], 'emission': [], 'absorption': []})

    _sampleType = 'empty'
    _refractiveIndex = 1
    _beamWidth = 0*100e-4 # cm
    _cuvetteTransmitance = 0.875
    _fittedEmissionModel = None

    _Diluter = None # sample object
    _Empty = None # sample object
    _BeamProfile = None # Beam profile object
    _dataDetails = pd.Series(dtype='float64')
    _metadata = pd.DataFrame()
    _sampleDetails = pd.Series()
    _chsDictionary = {}
    _dataPath = ''
    _metadataPath = ''
    _samplesCatalogPath = ''
    _dataFileName = ''
    _dataFullPath = ''
    _dataID = ''
    _plotSettings = {}
    _PMCalibrationFactors = [305.70189599790564, 0.5910134151908153] # [mW/V, mW]
    _APDCalibrationFactors = [] #

    _usingSteady = False
    _sampleAbsorbanceMean = None
    _sampleAbsorbanceStd = None
    _absorbanceMean = None
    _absorbanceStd = None
    _musAtAbsorptionWL = 0
    _musAtEmission = 0 # scaterring at emission wavelength to account for emission attenuation through the cuvette.

    def __init__(self, path, expId, sampleType):
 
        assert sampleType in ['ucnp', 'dye', 'diluter', 'empty'], "sampleType should be one of: 'ucnp', 'dye', 'diluter', 'empty'"
        self._settype(sampleType)
        self._plotSettings = {'m': 'o', 'mc': 'k', 'mfc': 'w', 'colour': 'k'}
        self._dataPath = path
        self._dataFileName = f"qy_{expId.replace('-', '_')}.csv"
        self._dataFullPath = path + self._dataFileName
        self._dataID = expId
        self._metadataPath = '/'.join(path.split('/')[:2])
        self._samplesCatalogPath = '/'.join([self._metadataPath, 'samples'])
        self._getmetadata()
        self._loaddata()
        self._getdetails()
        self._getsampledetails()
        self.clean()

    def _settype(self, sampleType):
        self._sampleType = sampleType

    def _loaddata(self):
        self._data = pd.read_csv(self._dataFullPath, usecols=self._metadata['field_name'])
        self._rawdata = self._data.copy()

    def _getmetadata(self):
        self._metadata = pd.read_csv(self._metadataPath + '/metadata.csv')

    def _updatemetadata(self, params):
        """Updates metadata. Sets a dataframe and concatenates with existing one.
        
        Args:
            params (dict): keys are columns names and values are list of items.
            Ex.: {'field_name': ['power]}
        """
        self._metadata = self._metadata.append(pd.DataFrame(params), ignore_index=True)

    def _getdetails(self):

        datalog = pd.read_csv(self._metadataPath + '/datalogs.csv')
        self._dataDetails = datalog[datalog['exp_id'] == self._dataID].iloc[0]
        # expand extra params
        extra_params = self._dataDetails['extra_params'].split('/')
        for param in extra_params:
            key = param.split('=')[0].replace(' ', '_')
            value = param.split('=')[1]
            self._dataDetails[key] = value
        
        # Channels are enumerated in the sorted order
        for item in enumerate(self._dataDetails['in_chs'].strip().lower().split('/')):
            self._chsDictionary[str(item[0])] = item[1].split('=')

    def _getplotlabel(self, col):
        """Gets plotlabel from metadata
        
        Args:
            col (str): column name
        Return:
            str: plot label string from col
        """
        return(
            self._metadata.loc[self._metadata['field_name']==col, 'plot_label'].item()
        )

    def _message(self, msg):
        """Send a message along data id.
        
        Args:
            msg (str): Success or explanation message.
        """
        return(f'{self._dataID}: {self._sampleType} - {msg}')

    def _getsampledetails(self):
        samplesCatalog = pd.read_csv(self._samplesCatalogPath + '/samples_catalog.csv')
        samplesCatalog = samplesCatalog.replace('-', np.nan)
        mask = samplesCatalog['Sample'] == self._dataDetails['sample']
        self._sampleDetails = samplesCatalog[mask].squeeze().copy()
        colsToCheck = ['QuantumYield', 'ExcitationWavelength_nm', 'EmissionWavelength_nm']
        self._sampleDetails[colsToCheck] = self._sampleDetails[colsToCheck].astype(float)
        if self._dataDetails['sample'] == 'empty':
            self._sampleDetails['Sample'] = 'empty'

    def setdiluter(self, diluterObj):
        '''Gives access to all the diluter info through a diluter object.

        Args:
            diluterObj (Sample): If the current sample is a empty or diluter type, 
            setting a diluter won't affect the calculation. 
        '''
        self._Diluter = diluterObj
    
    def setempty(self, emptyObj):
        '''Gives access to all the diluter info through a diluter object.

        Args:
            diluterObj (Sample): If the current sample is a empty or diluter type, 
            setting a diluter won't affect the calculation. 
        '''
        self._Empty = emptyObj

    def setbeamprofile(self, beamprofileObj):
        '''Gives access to all the beam profile info through a beam profile object.

        Args:
            beamprofileObj (BeamProfile): Start the class with the beam profile file name. 
        '''
        self._BeamProfile = beamprofileObj
        return(None, self._message('Success!'))

    def clean(self):
        """Rename cols according to the daq channels connections
        and calculates time based on the sampling rate (daq internal clock).

        Args:
            param1 (int):

        Returns:
            bool: The return value. True for success, False otherwise.
        """
        # get dictionary as ex.: {'0': 'apd', ...}
        colsDict = {k: v[1] for k, v in self._chsDictionary.items()}
        self._data = self._rawdata.rename(columns=colsDict)
        # update metadata
        for col, newCol in colsDict.items():
            self._metadata.loc[self._metadata['field_name'] == col, 'field_name'] = newCol
            plotLabel = f"{newCol.title() if len(newCol) > 3 else newCol.upper()} (V)"
            self._metadata.loc[self._metadata['field_name'] == newCol, 'plot_label'] = plotLabel
        
        # calculate clock time
        iniTime = 0
        # finalTime = len(self._data) / self._dataDetails['sampling_rate']
        # step = 1 / self._dataDetails['sampling_rate']
        step = self._dataDetails['time_per_step'] / self._dataDetails['samples_per_ch']
        finalTime = len(self._data) * step
        self._data['time_daq'] = np.arange(iniTime, finalTime, step)
        # update metadata
        if ~self._metadata['field_name'].isin(['time_daq']).any():
            self._metadata = self._metadata.append(
                pd.DataFrame({
                    'field_name': ['time_daq'],
                    'data_type': ['float'],
                    'data_format': ['^[-]?([0-9]?\.|[1-9][\d]*\.)[0-9]*$'],
                    'example': [0.02938],
                    'standard_units': ['s'],
                    'plot_label': 'Time (s)',
                    'description': 'Time acquired from the Daq internal clock'
                }), ignore_index=True
            )

        # print(f'Data: {self._dataID} update according to ->\n', self._chsDictionary)

    def removebackground(self, channel, baseline='trigger < 0.2'):
        '''Uses first points as background level.
    
        Args:
            channel (str): apd or pm
            baseline (str): condition to select the background value to be removed
        '''
        ch, operator, val = baseline.split()
        mask = eval(f"self._data['{ch}'] {operator} {val}") 
        backgroundlevel = self._data.loc[mask, channel].mean()
        self._data[channel] -= backgroundlevel

        return(backgroundlevel, self._message('Success!'))

    def details(self):
        return(self._dataDetails)

    def data(self):
        return(self._data)

    def setbeamwidth(self, beamWidth):
        """Beam width in cm."""
        self._beamWidth = beamWidth

    def filterpulse(self, precision=4, recalculate=False):
        """Filters data when trigger is set back to 0, as well as eventual outliers.

        Args:
            precision (float): std multiplier to select how far the point is from the fitting curve
            recalculate (bool): used to recalculate this paramenter
        """
        mask = ~((self._data['time_daq'] > 1.5 * self._dataDetails['time_per_step']) \
            & (self._data['trigger'] < self._dataDetails['range_step_size'] * 0.7))
        # eliminates data with trigger at 0 apart from initial 0 setting
        self._data = self._data[mask]

        # eliminating outliers
        coeffs, r_val = self.polyfit('time_daq', 'trigger', degree=1)

        p1 = np.array([0, -coeffs[1] / coeffs[0]])
        p2 = np.array([coeffs[1], 0])

        def distance(point):
            """Calculates distance between a point to the regression line.
            
            Args: 
                point (list of floats): xy, i.e. time, trigger list. Ex.: [time, trigger]
            Returns:
                float: distance
            """
            p = np.array(point)
            return(
                np.abs(np.linalg.norm(np.cross(p2-p1, p1-p))) / np.linalg.norm(p2-p1)
            )
        
        dist = self._data.apply(lambda df: distance([df['time_daq'], df['trigger']]), axis=1)
        self._data = self._data[dist < (dist.mean() + precision * dist.std())].reset_index(drop=True)

        return(None, self._message('Success!'))

    def steadydata(self, fromPt, toPt, precision=2, recalculate=False):
        """Groups the data per trigger point for a given precision. The is calculates the average of the points.
        Use the time plot to vizually identify the best area for the steady data to be evaluated.

        Args:
            fromPt (int): First point in which saturation of the equipment/sample is reached. 
                Use the number of samples per channel to roughly determine this number.
            toPt (int): Last point in which saturation of the equipment/sample is reached.
                Usually this is exactly the number of samples per channel.
            precision (int): The precision of the trigger values used to segregate them. 
                This depend on the step_size of the measurement. 
                Ex.: If step size is 0.001, the recommened precision is 3.
        """
        # newCol = 'apd_std'
        # calculate = recalculate or (newCol not in self._data.columns)
        # if calculate:
        #     if newCol in self._data.columns:
        #         self._data.drop(newCol, axis=1, inplace=True)
        #     elif ~self._metadata['field_name'].isin([newCol]).any():
        #         # update metadata
        #         descripStr = 'Same as transmitted power' if (channel == 'trigger') and (self._sampleType == 'empty') \
        #             else f'Obtained from {channel} calibration.'
        #         self._updatemetadata({
        #             'field_name': [newCol],
        #             'data_type': ['float'],
        #             'data_format': ['^[0-9]*.[0-9]*'],
        #             'example': [5.5],
        #             'standard_units': ['mW'],
        #             'plot_label': 'Power (mW)',
        #             'description': f'Apd std.'
        #         })

        calculate = recalculate or not self._usingSteady
        if calculate:
            self.clean()
            self._usingSteady = True
            self._data['temp_col'] = self._data['trigger'].round(precision)
            group = self._data.groupby('temp_col')
            self._data = group.nth(list(range(fromPt, toPt, 1))).groupby('temp_col').mean()
            apdStd = group.nth(list(range(500, 1000, 1))).groupby('temp_col').std()['apd'].rename('apd_std')
            pmStd = group.nth(list(range(500, 1000, 1))).groupby('temp_col').std()['pm'].rename('pm_std')
            self._data = pd.concat([self._data, apdStd, pmStd], axis=1).reset_index()
            self._data.drop('temp_col', axis=1, inplace=True)

            return(None, self._message('Success!'))
        else: 
            return(None, self._message('Steady data already calculated.'))

    def calibrate(self, calibFunction, channel='', recalculate=False):
        '''Calibrates a channel using a calibration function.
        Any negative resultant value will be replaced by 0.

        Args:
            calibFunction (function): function that converts voltage to power.
            channel (str): 'pm', 'trigger' or 'apd' are the option of channels to be calibrated.
                The power meter calibration (pm) requires a calibration function that correlates the 
                power meter analog output (in V) to power (in mW).
                
                The trigger calibration requires a calibration function that correlates the 
                trigger (daq) analog output (in V) to power (in mW). Fit an empty cuvette holder data
                with the power meter calibrated in order to get the calibration function. 
                Important: the empty data must be acquired on the same conditions (excitation arm) of the 
                intended sample, in which the calibration will be done. 

                The APD calibration is done using a reference dye. Use apd power vs apd signal in order
                to get the proper calibration function. 

            recalculate (bool): If True, the method will drop previous calculation and will recalculate it.         
        '''
        calibDictionary = {
            'pm': 'transmitted_power',
            'trigger': 'laser_power',
            'apd': 'apd_power'
        }
        if channel.lower() not in calibDictionary.keys():
            return(0, self._message('Select one of the channels: pm, trigger or apd'))
        
        newCol = calibDictionary[channel]
        descLabel = newCol.replace('_', ' ').title()
        calculate = recalculate or (newCol not in self._data.columns)
        if calculate:
            if newCol in self._data.columns:
                self._data.drop(newCol, axis=1, inplace=True)
            elif ~self._metadata['field_name'].isin([newCol]).any():
                # update metadata
                descripStr = 'Same as transmitted power' if (channel == 'trigger') and (self._sampleType == 'empty') \
                    else f'Obtained from {channel} calibration.'
                self._updatemetadata({
                    'field_name': [newCol],
                    'data_type': ['float'],
                    'data_format': ['^[0-9]*.[0-9]*'],
                    'example': [5.5],
                    'standard_units': ['mW'],
                    'plot_label': 'Power (mW)',
                    'description': f'{descLabel} - {descripStr}'
                })
        
            if (channel == 'trigger') and (self._sampleType == 'empty'):
                # assuming that empty is used to fit trigger vs transmitted_power beforehand
                self._data[newCol] = self._data['transmitted_power']
            else:
                self._data[newCol] = self._data[channel].apply(calibFunction)
                self._data.loc[self._data[newCol] < 0, newCol] = 0

        return(self._data[newCol], self._message('Success!'))

    def fitXY(self, x, y, condition='', show=False):
        """Fit lenearly cols x and y and plot results if resquested via show.
        (Will be deprecated soon!)
        Args:
            condition (str): condition to filter data to be fitted
                Ex.: 'x>1'
            show (bool): If true, plots results

        Returns:
            slope, intercept, r_val, p_val, std_err (float): 
            return results of fitting
        """
        # loads calibration data 
        if condition:
            mask = eval(condition.replace('x', 'self._data[x]'))
            dataToFit = self._data[mask]
        else: dataToFit = self._data

        slope, intercept, r_val, p_val, std_err = stats.linregress(
            dataToFit[x], dataToFit[y]
        )
        fitting = lambda val: slope * val + intercept

        if show:
            fig, ax = pl.subplots()
            ax.plot(self._data[x], self._data[y],
                label='Data',
                marker='s',
                ms=4,
                color='gray',
                mfc='white',
                alpha=0.5,
            )
            ax.plot(dataToFit[x], dataToFit[y],
                label='Calibration data',
                marker='s',
                ms=4,
                color='black',
                mfc='white',
                alpha=1.0,
            )
            ax.plot(self._data[x], fitting(self._data[x]),
                label='Fitting curve',
                color='firebrick',
                alpha=0.8,
                linestyle='-',
                lw=1.8
            )
            ax.legend()
            ax.set_title(f'y = {slope:.3f}.x + {intercept:.3f} || ' + r'$R^2$ = ' + f"{r_val**2:.6f}", 
                loc='left',
                fontsize= 12
            )
            ax.set_xlabel(self._getplotlabel(x))
            ax.set_ylabel(self._getplotlabel(y))
            fig.tight_layout()
            fig.show()

        # print(f"R^2:{r_val ** 2}")
        return(slope, intercept, r_val, p_val, std_err)

    def polyfit(self, x, y, degree, condition='', show=False):
        """Polynomial fit cols x and y and plot results if resquested via show.

        Args:
            x (str):
            y (str):
            degree (int):
            condition (str): condition to filter data to be fitted
                Ex.: 'x>1'
            show (bool): If true, plots results

        Returns:
            coeffs, r_val, (float, float): Results of fitting
        """
        if condition:
            mask = eval(condition.replace('x', 'self._data[x]'))
            dataToFit = self._data[mask]
        else: dataToFit = self._data

        coeffs = np.polyfit(dataToFit[x], dataToFit[y], degree).tolist()
        r_val = np.corrcoef(dataToFit[x], dataToFit[y])[0,1]
        fitting = np.poly1d(coeffs)

        if show:
            fig, ax = pl.subplots()
            ax.plot(self._data[x], self._data[y],
                label='Data',
                marker='s',
                ms=4,
                color='gray',
                mfc='white',
                alpha=0.5,
            )
            ax.plot(dataToFit[x], dataToFit[y],
                label='Calibration data',
                marker='s',
                ms=4,
                color='black',
                mfc='white',
                alpha=1.0,
            )
            ax.plot(self._data[x], fitting(self._data[x]),
                label='Fitting curve',
                color='firebrick',
                alpha=0.8,
                linestyle='-',
                lw=1.8
            )
            ax.legend()

            xList = [f'$x^{n + 1}$' for n in range(len(coeffs) - 1)[::-1]]
            xList.append('')
            yEq = 'y = ' + ' + '.join([f"{c:.3} " + x for x, c in zip(xList, coeffs)])
            
            ax.set_title(f'{yEq}' + '\n' + r'$R^2$ = ' + f"{r_val**2:.6f}", 
                loc='right',
                fontsize= 12
            )
            ax.set_xlabel(self._getplotlabel(x))
            ax.set_ylabel(self._getplotlabel(y))
            fig.tight_layout()
            fig.show()

        return(coeffs, r_val)

    def transmittedpower(self):
        """Power transmitted thorugh the sample chamber.
        Returns: 
            (transmitted_power (pd.Series) | None, str).
        """
        if 'transmitted_power' in self._data.columns:
            return(self._data['transmitted_power'], self._message('Success!'))
        else:
            return(None, self._message('Calibrate the pm channel before proceding.'))

    def laserpower(self):
        """Laser Power column.
        Returns: 
            (pd.Series | None, str): Laser power column or None and success or failure message.
        """
        if 'laser_power' in self._data.columns:
            return(self._data['laser_power'], self._message('Success!'))
        else:
            return(None, self._message('Calibrate the trigger channel before proceding.'))
        
    def poweratcentre(self, recalculate=False):
        r"""Calculates power at the centre of the cuvette.
        Given by $\sqrt(P_o * P)$
        
        Args:
            recalculate (bool): if True it recalculates the power at the centre of the cuvette.
        """
        newCol = 'power_at_centre'
        calculate = recalculate or (newCol not in self._data.columns)
        if calculate:
            # requires: transmitted_power and laser_power
            transmittedPow, msg1 = self.transmittedpower()
            laserPow, msg2 = self.laserpower()
            if (transmittedPow is None) or (laserPow is None):
                return(None, '\n'.join([msg1, msg2]))

            if newCol in self._data.columns:
                self._data.drop(newCol, axis=1, inplace=True)
            else:
                # update metadata
                self._updatemetadata({
                    'field_name': ['power_at_centre'],
                    'data_type': ['float'],
                    'data_format': ['^[0-9]*.[0-9]*'],
                    'example': [100.565],
                    'standard_units': ['mW'],
                    'plot_label': 'Power (mW)',
                    'description': 'Power at the center of the cuvette.'
                })
            
            self._data['power_at_centre'] = np.sqrt(transmittedPow * laserPow)
        
        return(self._data[newCol], self._message('Success!'))

    def powdensatcentre(self, beamWidth=None, recalculate=False):
        """Calculates power density at the centre of the cuvette.
        Power at the centre/beam area.
        
        Args:
            recalculate (bool): if True it recalculates the power at the centre of the cuvette.
        """
        
        newCol = 'pow_dens_at_centre'
        calculate = recalculate or (newCol not in self._data.columns)
        if calculate:
            # requires: power_at_centre
            powerAtCentre, msg = self.poweratcentre(recalculate=True)
            if powerAtCentre is None:
                return(None, msg)

            beamWidth = beamWidth or self._BeamProfile.beamwidth() or self._beamWidth
            if beamWidth is None:
                return(None, self._message(
                    'Beam width not found. '\
                    + 'Set it with setbeamwidth method or pass it explicitly to powerdensatcentre method.'
                ))
            else: self.setbeamwidth(beamWidth)
            
            if newCol in self._data.columns:
                self._data.drop(newCol, axis=1, inplace=True)
            else:
                # update metadata
                self._updatemetadata({
                    'field_name': [newCol],
                    'data_type': ['float'],
                    'data_format': ['^[0-9]*.[0-9]*'],
                    'example': [10000.565],
                    'standard_units': ['W/cm^2'],
                    'plot_label': r'Power ($W/cm^2$)',
                    'description': 'Power density at the center of the cuvette.'
                })
            self._data[newCol] = powerAtCentre \
                / (1000.0 * np.pi * (beamWidth / 2) ** 2) # W/cm^2

        return(self._data[newCol], self._message('Success!'))

    def sampleabsorbance(self, filterQuery='laser_power > 0', recalculate=False):
        '''Calculates absorbance per cm of the whole sample using power before and after cuvette.
        !! note: this is in ln instead of log and should be reviewed!
        Args:
            mean (bool): If True, it tries to calculate the mean value based on filterQuery.
            filterQuery (str): query to filter data before calculating mean of absorbance.
                It must contain single space between column, operator and value. Ex.: 'column > 1.2' 
        Returns:
            pandas.DataFrame, float, float, str: Data, Absorbance mean, std and success, failure message. 
        '''
        newCol = 'sample_absorbance'
        calculate = recalculate or (newCol not in self._data.columns)
        if calculate:
            # requires: transmitted_power and laser_power
            transmittedPow, msg1 = self.transmittedpower()
            laserPow, msg2 = self.laserpower()
            if (transmittedPow is None) or (laserPow is None):
                return(None, None, None, '\n'.join([msg1, msg2]))

            if newCol in self._data.columns:
                self._data.drop(newCol, axis=1, inplace=True)
            else:
                # update metadata
                self._updatemetadata({
                    'field_name': [newCol],
                    'data_type': ['float'],
                    'data_format': ['^0.[0-9]*'],
                    'example': [0.4],
                    'standard_units': ['$(cm^{-1})$'],
                    'plot_label': 'Absorbance $(cm^{-1})$',
                    'description': 'Absorbance obtained by log(p_o/p), where p_o is power before sample and p is power after. It refers to whatever is inside the cuvette.'
                })

            # if newCol not in self._data.columns:
            mask = laserPow > 0
            self._data.loc[mask, newCol] = np.log(laserPow[mask] / transmittedPow[mask])
            # calculates mean and std
            col, operator, value = filterQuery.split(' ')
            mask = eval(f"self._data['{col}'] {operator} {value}")
            self._sampleAbsorbanceMean = self._data.loc[mask, newCol].mean() 
            self._sampleAbsorbanceStd = self._data.loc[mask, newCol].std()

        return(
            self._data[newCol], 
            self._sampleAbsorbanceMean, 
            self._sampleAbsorbanceStd,
            self._message('Success!')
        )

    def absorbance(self, diluterAbsorbance=None, filterQuery='laser_power > 0', recalculate=False):
        '''Calculates absorbance of the emitter compound present on the aqueous medium.
        
        Args: 
            diluterAbsorbance (float): Diluter absorbance used as a reference to calculate the emitter's absorbance
            recalculate (bool): If True it will recalculate the column

        Returns:
            pandas.DataFrame, float, float, str: Data, Absorbance mean, std and success, failure message.
        '''
        newCol = 'absorbance'
        calculate = recalculate or (newCol not in self._data.columns)
        if calculate:
            sampleAbsorbance, _, _, msg = self.sampleabsorbance(recalculate=False)
            if sampleAbsorbance is None:
                return(None, None, None, msg)
            if (self._sampleType == 'diluter') or (self._sampleType == 'empty'):
                self._data[newCol] = sampleAbsorbance
                self._absorbanceMean = self._sampleAbsorbanceMean
                self._absorbanceStd = self._sampleAbsorbanceStd
            else:
                # check diluter
                if diluterAbsorbance is None:
                    if self._Diluter is None:
                        return(None, None, None, self._message('Diluter not found. Set it using setdiluter method or pass its value explicitaly.'))
                    else: 
                        # check this part. Aparently I'm not using the filtered sample absorbance
                        _, dilAbsorbMean, _, msg = self._Diluter.sampleabsorbance() 
                        if dilAbsorbMean is None:
                            return(None, None, None, 
                                '\n'.join([msg, self._message('Alternatively, you can pass the diluters absorbance explicitaly.')])
                            )
                        else: diluterAbsorbance = dilAbsorbMean

                

                if newCol in self._data.columns:
                    self._data.drop(newCol, axis=1, inplace=True)
                else:
                    # update metadata
                    self._updatemetadata({
                        'field_name': [newCol],
                        'data_type': ['float'],
                        'data_format': ['^0.[0-9]*'],
                        'example': [0.4],
                        'standard_units': ['$cm^{-1}$'],
                        'plot_label': 'Absorbance ($cm^{-1}$)',
                        'description': 'Absorbance obtained by mixture absorbance - diluter absorbance.'
                    })
                # last change 13/05/2021: scattering included
                # self._data[newCol] = sampleAbsorbance - diluterAbsorbance
                self._data[newCol] = sampleAbsorbance - diluterAbsorbance - self._musAtAbsorptionWL
                # calculates average and std
                col, operator, value = filterQuery.split(' ')
                mask = eval(f"self._data['{col}'] {operator} {value}")
                self._absorbanceMean = self._data.loc[mask, newCol].mean()
                self._absorbanceStd = self._data.loc[mask, newCol].std()

        return(
            self._data[newCol], 
            self._absorbanceMean, 
            self._absorbanceStd,
            self._message('Success!') 
        )
        
    def absorbedpower(self, region=[0.45, 0.55], useMean=True, absorptionCoeff=None, filterQuery='laser_power > 0', recalculate=False):
        """Calculates the total power absorbed by the dye on the region: 
        [0.45 to 0.55] cm inside the cuvette, 1 mm slit aperture.
            AbsP(xi, xf) = Pxo (A_u / A_m) * (e^(-A_m * xi) - e^(-A_m * xf))
            With Pxo = Po * e ^ (-A_m * xo) and xo = 0

        Args:
            region ([float, float]): Values in cm between 0 and 1.
            useMean (bool):
            filterQuery (str): used only when useMean is True
            recalculate (bool): If True it will recalculate the column
        """
        newCol = 'absorbed_power'
        calculate = recalculate or (newCol not in self._data.columns)
        if calculate:
            if self._sampleType == 'empty':
                self._data[newCol] = 0
            else:
                absorbDF, absorbMean, _, msg1 = self.absorbance(filterQuery=filterQuery, recalculate=True)
                sampleAbsorbDF, sampAbsorbMean, _, msg2 = self.sampleabsorbance(filterQuery=filterQuery, recalculate=True)
                if absorptionCoeff:
                    absorbance = absorptionCoeff
                else:
                    absorbance = absorbMean if useMean else absorbDF
                sampleAbsorbance = sampAbsorbMean if useMean else sampleAbsorbDF
                laserPow, msg3 = self.laserpower()
                if (absorbance is None) or (sampleAbsorbance is None):
                    return(None, '\n'.join([msg1, msg2, msg3]))

                if newCol in self._data.columns:
                    self._data.drop(newCol, axis=1, inplace=True)
                else:
                    # update metadata
                    self._updatemetadata({
                        'field_name': [newCol],
                        'data_type': ['float'],
                        'data_format': ['^0.[0-9]*'],
                        'example': [0.4],
                        'standard_units': ['mW'],
                        'plot_label': f'Abs. pow. $(mW)$',
                        'description': f'Absorbed power calculated for the region:  xi={region[0]} and xf={region[1]}. AbsP(xi, xf) = Pxo (A_u / A_m) * (e^(-A_m * xi) - e^(-A_m * xf)). With Pxo = Po * e ^ (-A_m * xo) and xo = 0.'
                    })
                
                windowTransmitance = (self._cuvetteTransmitance) ** 0.5
                self._data[newCol] = laserPow * windowTransmitance * (absorbance / sampleAbsorbance) \
                    * (np.exp(- sampleAbsorbance * region[0]) - np.exp(- sampleAbsorbance * region[1]))

        return(self._data[newCol], self._message('Success!'))

    def abspowerperarea(self, recalculate=False):
        '''Divides the total absorbed power by the area of the beam profile.
        Useful to compare the result of two different beams
        Args:
            recalculate (bool): If True it will recalculate the column
        '''
        newCol = 'absorbed_power_per_area'
        calculate = recalculate or (newCol not in self._data.columns)
        if calculate:
            absorbedPow, msg1 = self.absorbedpower()
            beamArea = self._BeamProfile.beamarea()

            if (absorbedPow is None):
                return(None, msg1)
    
            if newCol in self._data.columns:
                self._data.drop(newCol, axis=1, inplace=True)
            else:
                self._updatemetadata({
                    'field_name': [newCol],
                    'data_type': ['float'],
                    'data_format': ['^[0-9]*.[0-9]*'],
                    'example': [10.565],
                    'standard_units': ['mW/cm^2'],
                    'plot_label': 'Abs. Pow. $(mW/cm^2)$',
                    'description': 'Absorbed power per area in cm^2.'
                })

            self._data[newCol] = absorbedPow / beamArea

        return(self._data[newCol], self._message('Success!'))

    def emittedpower(self, scatteringAtEmission=0, recalculate=False):
        """
        for dye: absorb * qy * excitationWL / emissionWL
        fit emitted power x apd
        Args:
            recalculate (bool): If True it will recalculate the column
        """

        newCol = 'emitted_power'
        calculate = recalculate or (newCol not in self._data.columns)
        if calculate:
            if self._sampleType == 'dye':
                absorbedPow, msg1 = self.absorbedpower()
                quantumYield, msg2 = self.quantumyield()
                emissionWL = self._sampleDetails['EmissionWavelength_nm']
                excitationWL = self._sampleDetails['ExcitationWavelength_nm']
                if (absorbedPow is None) or (quantumYield is None):
                    return(None, '\n'.join([msg1, msg2]))
            else:
                apdPow, msg = self.apdpower()
                filterTransmitance = equipment['filter900TransmitanceAt800']
                self._musAtEmission = scatteringAtEmission
                if apdPow is None:
                    return(None, msg)
                
            if newCol in self._data.columns:
                self._data.drop(newCol, axis=1, inplace=True)
            else:
                self._updatemetadata({
                    'field_name': ['emitted_power'],
                    'data_type': ['float'],
                    'data_format': ['^[0-9]*.[0-9]*'],
                    'example': [10.565],
                    'standard_units': ['mW'],
                    'plot_label': 'Emission $(mW)$',
                    'description': 'Total emitted power in all directions.'
                })

            self._data[newCol] = absorbedPow * quantumYield * excitationWL / emissionWL if self._sampleType == 'dye' else apdPow * np.exp(scatteringAtEmission * 1 / 2) / filterTransmitance
            # water absorbance at 800nm not taken into account here. It should change emission in around 4%
            # scattering is taken into account from another measurement.

        return(self._data[newCol], self._message('Success!'))

    def _emissionModel(self, power, phib, rhob):
        """
        """
        # review this equation
        emissionWL = self._sampleDetails['EmissionWavelength_nm']
        excitationWL = self._sampleDetails['ExcitationWavelength_nm']
        _, profile = self._BeamProfile.histprofile(normalised=True)
        powProf = np.outer(power, profile[:-1])
        return(2 * phib * self.aCoeff * (powProf ** 2) * (excitationWL / emissionWL) / (self._BeamProfile.camPxSize ** 2 * rhob + powProf))

    def _emissionModel1(self, power, phib, rhob, rba):
        """Considers the std rhob = rhob * rba
        """
        _, profile = self._BeamProfile.histprofile(normalised=True)
        powProf = np.outer(power, profile[:-1])
        return(phib * self.aCoeff * (powProf ** 2) / ((self._BeamProfile.camPxSize ** 2 * rba * rhob + powProf) * (self._BeamProfile.camPxSize ** 2 * rba + powProf)))

    def _emissionModel2(self, power, phib, rhob, a, b):
        """Considers the std rhob = rhob * rba
        """
        _, profile = self._BeamProfile.histprofile(normalised=True)
        powProf = np.outer(power, profile[:-1])
        phi_b * x /(rho_b + a/(x+b) + x)
        return(phib * self.aCoeff * (powProf ** 2) / (self._BeamProfile.camPxSize ** 2 * (rhob + a / (powProf + b)) + powProf))

    def fitemission(self, bounds=(0, [1, np.inf]), model='std', show=False):
        '''
        '''

        self.aCoeff, pcov = sciopt.curve_fit(
            lambda aCoeff, power: aCoeff * power,
            self._data['power_at_centre'],
            self._data['absorbed_power']
        )
 
        counts, _ = self._BeamProfile.histprofile(normalised=True)
        if model == 'std':
            emission = lambda power, phib, rhob: np.sum(counts * self._emissionModel(power, phib, rhob), axis=1)
        elif model == 'mod':
            emission = lambda power, phib, rhob, rba: np.sum(counts * self._emissionModel1(power, phib, rhob, rba), axis=1)

        popt, pcov = sciopt.curve_fit(emission, self._data['power_at_centre'], self._data['emitted_power'], bounds=bounds)
        self._fittedEmissionModel = lambda power: self._emissionModel(power, *popt)
        emissionFunc = lambda power: emission(power, *popt)

        if show:
            fig, ax = pl.subplots()
            ax.plot(self._data['power_at_centre'], self._data['emitted_power'],
                label='Data',
                marker='s',
                ms=4,
                color='black',
                mfc='white',
                alpha=1.0,
            )
            ax.plot(self._data['power_at_centre'], emission(self._data['power_at_centre'], *popt),
                label='Fitting curve',
                color='firebrick',
                alpha=0.8,
                linestyle='-',
                lw=1.8
            )
            ax.legend()
            # ax.set_title(f'y = {slope:.3f}.x + {intercept:.3f} || ' + r'$R^2$ = ' + f"{r_val**2:.6f}", 
            #     loc='left',
            #     fontsize= 12
            # )
            ax.set_xlabel(self._getplotlabel('power_at_centre'))
            ax.set_ylabel(self._getplotlabel('emitted_power'))
            fig.tight_layout()
            fig.show()

        return(popt, pcov, emissionFunc, self._message('Success!'))

    def deconvemission(self, recalculate=False):
        """
        """
        calculate = recalculate or (self._dataDeconv.empty)
        if calculate:
            if self._dataDeconv.empty:
                self._updatemetadata({
                    'field_name': ['power_density'],
                    'data_type': ['float'],
                    'data_format': ['^[0-9]*.[0-9]*'],
                    'example': [5.5],
                    'standard_units': ['W/cm^2'],
                    'plot_label': 'Power $(W/cm^{2})$',
                    'description': 'Power density - Obtained from beam profile deconvolution using the pixels area of the camera.'
                })
                self._updatemetadata({
                    'field_name': ['absorption'],
                    'data_type': ['float'],
                    'data_format': ['^[0-9]*.[0-9]*'],
                    'example': [5.5],
                    'standard_units': ['mW/cm^2'],
                    'plot_label': 'Absorption $(mW/cm^{3})$',
                    'description': 'Absorption - Absorbed power per volume unit in cm^3.'
                })
                self._updatemetadata({
                    'field_name': ['emission'],
                    'data_type': ['float'],
                    'data_format': ['^[0-9]*.[0-9]*'],
                    'example': [5.5],
                    'standard_units': ['mW/cm^2'],
                    'plot_label': 'Emission $(mW/cm^{3})$',
                    'description': 'Emission - Emitted power per volume unit in cm^3. The area used is the pixel area.'
                })
            # power, emission
            self._dataDeconv = {'power_density': [], 'emission': [], 'absorption': []}

            for power in self._data['power_at_centre']:
                counts, profile = self._BeamProfile.histprofile(normalised=True)
                powDensProf = power * profile[:-1] / (self._BeamProfile.camPxSize ** 2)
                self._dataDeconv['absorption'] += list(powDensProf * self.aCoeff / 0.1) # 0.1 cm length or 1 mm
                self._dataDeconv['power_density'] += list(powDensProf / 1000) # converts from mW/cm^2 to W/cm^2
                self._dataDeconv['emission'] += list(self._fittedEmissionModel(power)[0] / (0.1 * self._BeamProfile.camPxSize ** 2))

            self._dataDeconv = pd.DataFrame(self._dataDeconv)
            emissionWL = self._sampleDetails['EmissionWavelength_nm']
            excitationWL = self._sampleDetails['ExcitationWavelength_nm']
            self._dataDeconv['quantum_yield'] = self._dataDeconv['emission'] * emissionWL / (excitationWL * self._dataDeconv['absorption'])
        
        return(self._dataDeconv, self._message('Success!'))

    def apdpower(self, recalculate=False):
        """Calculates power after filter attenuation in 360 degrees solid angle.
        Indeed the power reaching the apd is only a tiny fraction of total emitted power, 
        however this is fine for the purpose of the quantum yield calculation. 

        Args:
            filterTransmitance (float): Lower than 1, indicates the filter attenuation of the emitted light by the sample.
            recalculate (bool): If True it will recalculate the column
        """
        newCol = 'apd_power'
        calculate = recalculate or (newCol not in self._data.columns)
        if calculate:
            if self._sampleType != 'dye':
                return(None, self._message('Use the calibrate method to calculate the Apd power.'))
            else:
                reabsorption = 0.103 # 10% for the dye-781. This shouldn't be declared here in future to make it easier to read the code and update it. 
                # calculated in notebook: dye-reabsorption
                emittedPow, msg = self.emittedpower()
                print('ok')
                filterTransmitance = equipment['filter785TransmitanceAt800']
                if emittedPow is None:
                    return(None, msg)
                if newCol in self._data.columns:
                    self._data.drop(newCol, axis=1, inplace=True)
                else:
                    self._updatemetadata({
                        'field_name': [newCol],
                        'data_type': ['float'],
                        'data_format': ['^[0-9]*.[0-9]*'],
                        'example': [5.5],
                        'standard_units': ['mW'],
                        'plot_label': 'APD Power (mW)',
                        'description': 'APD power - Power reaching the apd expanded to all space.'
                    })
                self._data[newCol] = emittedPow * filterTransmitance * (1 - reabsorption)

        return(self._data[newCol], self._message('Success!'))
            
    def quantumyield(self, recalculate=False):
        """Calculates the quantum yield of the sample by:
            ```QY = E lambda_e / (A lambda_a)```
            With E the total emission at the wavelength nu_e and A the total absorbed power at the wavelength nu_a. 

        Args:
            recalculate (bool): If True it will recalculate the column
        """
        newCol = 'quantum_yield'
        calculate = recalculate or (newCol not in self._data.columns)
        if calculate:
            if self._sampleType == 'dye':
                return(self._sampleDetails['QuantumYield'], self._message('Success!'))
            elif self._sampleType != 'ucnp':
                return(0, self._message('Success!'))
            else:
                emissionWL = self._sampleDetails['EmissionWavelength_nm']
                excitationWL = self._sampleDetails['ExcitationWavelength_nm']
                emittedPow, msg1 = self.emittedpower()
                absorbedPow, msg2 = self.absorbedpower()
                if (emittedPow is None) or (absorbedPow is None):
                    return(None, '\n'.join([msg1, msg2]))
            if newCol in self._data.columns:
                self._data.drop(newCol, axis=1, inplace=True)
            else:
                self._updatemetadata({
                    'field_name': [newCol],
                    'data_type': ['float'],
                    'data_format': ['^[0-9]*.[0-9]*'],
                    'example': [10000.565],
                    'standard_units': ['-'],
                    'plot_label': r'$\phi$',
                    'description': 'Quantum yield: (E * lambda_e) / (A * lambda_a) '
                })
            self._data[newCol] = emittedPow * emissionWL / (absorbedPow * excitationWL)

        return(self._data[newCol], self._message('Success!'))

    def kappa(self, period=5, recalculate=False):
        """Emission slope in log scale.
        
        Args:
            period (int): Number of points for linear approximation
            recalculate (bool): If true it recalculates the property
        """
        newCol = 'kappa'
        calculate = recalculate or (newCol not in self._data.columns)
        if calculate:
            emittedPow, msg1 = self.emittedpower()
            powDensAtCentre, msg2 = self.powdensatcentre()
            if (emittedPow is None) or (powDensAtCentre is None):
                return(None, '\n'.join([msg1, msg2]))
            
            if newCol in self._data.columns:
                self._data.drop(newCol, axis=1, inplace=True)
            else:
                # update metadata
                self._updatemetadata({
                    'field_name': [newCol],
                    'data_type': ['float'],
                    'data_format': ['^[0-9]*.[0-9]*'],
                    'example': [10000.565],
                    'standard_units': ['-'],
                    'plot_label': r'$\kappa$',
                    'description': 'Emission slope in the logarithmic scale.'
                })

            # calculate kappa
            df = self._data.copy()
            df['index_'] = df.index

            df = df[df['pow_dens_at_centre'] > 0].set_index('pow_dens_at_centre')
            # log10 or log?
            getslope = lambda df: np.polyfit(np.log10(df.index.values), np.log10(df.values), 1)[0]
            kappa = df['emitted_power'].rolling(period).apply(getslope)

            kappa = pd.concat([kappa, df['index_']], axis=1)

            kappa.set_index(kappa['index_'], inplace=True)
            kappa.rename(columns={'emitted_power': newCol}, inplace=True)
            kappa.drop('index_', axis=1, inplace=True)

            self._data = pd.concat([self._data, kappa], axis=1)
        
        return(self._data[newCol], self._message('Success!'))

    def view(self, x='time', yList=['trigger', 'pm', 'apd'], label='', deconv=False):
        '''
        
        Args:
            label (str): label+label2
            deconv (bool): if true it uses _dataDeconv dataframe to search the data to be plotted
        '''
        fig01, axs = pl.subplots(len(yList), sharex=True, figsize=[7,6])
        dataToPlot = self._dataDeconv if deconv else self._data
        
        axs = [axs] if len(yList) == 1 else axs
    
        for n in range(len(yList)):
            axs[n].plot(dataToPlot[x], dataToPlot[yList[n]],
                label=f"{self._details()[label]}" if label else f"{self._dataID}",
                marker=self._plotSettings['m'],
                ms=4,
                color=self._plotSettings['colour'],
                mfc='w',
                alpha=0.7,
                linestyle='' if deconv else '-',
                lw=1.0
            )

        for n in range(len(yList)):
            ylabel = self._metadata.loc[self._metadata['field_name']==yList[n], 'plot_label'].item()
            axs[n].set_ylabel(ylabel)
            axs[n].grid()
        handles, labels = axs[0].get_legend_handles_labels()
        fig01.legend(handles, labels, title=label if label else 'Exp Id', loc='center', bbox_to_anchor=(0.88, 0.6),
            ncol=1, fancybox=True, shadow=True, fontsize=12)

        xlabel = self._metadata.loc[self._metadata['field_name']==x, 'plot_label'].item()
        axs[-1].set_xlabel(xlabel)

        # axs[0].set_title(f'Measurement ID: {self._dataID}')
        fig01.tight_layout()
        fig01.subplots_adjust(right=0.76)
        fig01.show()
        print(self._dataDetails)
        return(fig01, axs)
    

class Analysis():
    """Creates and handle different Samples.
    Ex.:
    >>> exp = Analysis()
    >>> exp.loaddata(
        sampleType='dye',
        dataIds=['20201111-101010'],
        diluterId='20201111-101110',
        emptyId='20201111-101200'
    )

    **Analysis protocol**:
        In order to get the quantum yield values for a ucnp sample, 
        the steps bellow should be followed for a reference dye data
        and subsequently for ucnps data. 

        1. Get steady data - it averages data points after the equipmment and samples stabilised. 
        2. Calibrate the power meter output (pm) - Volts to mW using standard curves for that.
            This is the transmitted power through the cuvette.
        3. Calibrate the DAQ output (trigger) - Volts to mW using empty measuremnts.
            This is the laser power, i.e. the power before reaching the cuvette.
        4. Get the power density in the center - It uses the power before and after the cuvette 
            along the beam width.
        5. Get absorbance - Sample-Absorbance is the abosorbance of the whole sample.
            On the other hand Absorbance is the absorbance of the compound without its diluter.
            It's necessary to have the diluter abosrobance in order to calculate this property.
        6. Get absorbed power - Using absorbance it calculates the total power absorbed by the 
            compound in a certain range within the cuvette. 
        **Dye**:
        7. Get emitted power - Knowing its quantum yield the emitted power is calculated.
        8. Get apd power - It calculates the power reaching the apd after the light passes
            through a long pass filter.
        9. Fit apd power vs apd voltage (apd) - It calibrates the apd output.
        **UCNP**:
        7. Calibrate the apd - Using a dye calibration curve it calculates the power reaching the 
            apd.
        8. Get emitted power - Calculates the power before the short pass filter, i.e. emitted power.
        9. Get kappa - Calculates the slope of the emission curve vs power density at centre in log
            scale.
        10. Get quantum yield - Calculates the quantum yield of a sample from its emitted power and
            aborbed power.
    
    """

    _details = pd.DataFrame(columns=['type'])
    _sampleList = np.array([])
    datalog = pd.read_csv('../data/datalogs.csv')
    _dataPath = ''
    _acceptedTypes = ['ucnp', 'dye', 'diluter', 'empty']
    _markersList = ['o', 'v', '>', 's', '<', '^', '*', 'p', '+']
    _coloursList = [
        'rebeccapurple',
        'firebrick', 
        'black', 
        'magenta',
        'crimson', 
        'orangered', 
        'slateblue',
        'seagreen', 
        'teal',
        'gray',
        'saddlebrown'
    ]

    def __init__(self, dataPath='../data/raw-data/'):
        """Builds an object for analysis

        Args: 
            dataPath (str): Path to the data files.
        """
        self._dataPath = dataPath
        pl.rcParams.update({'font.size': 15})
        # self.datalog = pd.read_csv('../data/datalogs.csv')

    def _load(self, sampleType, dataId):
        """Loads a sample and stores in a numpy array.

        Args: 
            sampleType (str): 'ucnp', 'dye', 'diluter' or 'empty'.
            dataId (str): id of the data.
        """
        if dataId is not None:
            cnt = len(self._sampleList)
            if self.search(dataId) is None:
                self._sampleList = np.append(self._sampleList, Sample(path=self._dataPath, expId=dataId, sampleType=sampleType))
        
            dataIdx = self.search(dataId)
            self._sampleList[dataIdx]._plotSettings['m'] = \
                self._markersList[cnt % len(self._markersList)]
            self._sampleList[dataIdx]._plotSettings['colour'] = \
                self._coloursList[cnt % len(self._coloursList)]
            # update details
            self.details()
        
    def _isid(self, sampleId):
        if type(sampleId) is str:
            return(re.match('^202[0-9][0-9]{4}-[0-9]{6}$', sampleId))
        else: return(None)

    def search(self, sampleId=None):
        """Search for sample with given id.

        Args:
            sampleId (str): id of experiment
        Returns:
            (int | None): Index of sample position or None in case nothing is found
        """
        # if sampleId in not None:
        assert self._isid(sampleId), 'The id should have the format: ########-######'
        for sample in self._sampleList:
            if sample._dataID == sampleId:
                result = np.where(self._sampleList == sample)
                assert len(result[0]) == 1, f'Found multiple samples for {sampleId}'
                idx = result[0][0]
                return(idx)
            
        return(None)
        # elif sampleType is not None:
        #     assert sampleType in self._acceptedTypes, f'Select one of: {self._acceptedTypes}.'
        #     for sample in self._sampleList:
        #         if sample._sampleType == sampleType:
        #             result = np.where(self._sampleList == sample)
        #             idxs = result[0]
        #             return(idxs)
        #     return(None)

    def loaddata(self, sampleType, dataIds, diluterId=None, emptyId=None):
        """
        
        """
        dataIds = dataIds if type(dataIds) is list else [dataIds]
        for sampleId in dataIds + [diluterId, emptyId]:
            if sampleId is not None:
                assert self.datalog['exp_id'].str.contains(sampleId).any(), f'{sampleId} not found on datalog'
                
        self._load('empty', emptyId)
        self._load('diluter', diluterId)
        if diluterId is not None:
            self.sample(diluterId).setempty(self.sample(emptyId))
        for dataId in dataIds:
            self._load(sampleType, dataId)
            if diluterId: self.sample(dataId).setdiluter(self.sample(diluterId))
            if emptyId: self.sample(dataId).setempty(self.sample(emptyId))

    def details(self, cols=None, full=False):
        """Shows the details of the experiments. 

        Args:
            cols (list): List of columns to show the details
            full (bool): If None, it shows predefined info
        """
        stdCols = ['exp_id', 'type', 'sample']
        colsToShow = ['apd_gain', 'comments', 'density_filter',
            'laser_wavelength', 'power_meter_range', 'range_start',
            'range_end', 'range_step_size',
            'samples_per_ch', 'sampling_rate', 'time_per_step'
        ] if cols is None else cols
        colsToShow = stdCols + colsToShow
        self._details = pd.DataFrame(columns=['type'])
        for sample in self._sampleList:
            self._details = self._details.append(
                sample.details(),
                ignore_index=True
            )
            self._details.iloc[-1, self._details.columns.get_loc('type')] = sample._sampleType
        if full:
            return(self._details)
        else:
            return(self._details[colsToShow])

    def sampleidx(self, id_x=None):
        '''Look up for sample index given its index or id.

        Args:
            id_x (str | int): id or data index
        Returns:
            (int): The found sample index.
        '''
        if (type(id_x) is int) and (0 <= id_x < len(self._sampleList)):
            idx = id_x
        elif self._isid(id_x):
            idx = self.search(id_x)
        else: idx = None
        
        assert idx is not None, f'Sample {id_x} not found. Try one of those:\n{self.details()}'
        return(idx)
        
    def sample(self, id_x):
        '''Look up for sample given its index or id.

        Args:
            id_x (str | int): id or data index
        Returns:
            (pandas.DataFrame): The found sample index and sample.
        '''
        return(self._sampleList[self.sampleidx(id_x)])
    
    def drop(self, id_x):
        """Drops the data selected by ID or idx, then updates details.

        Args:
            id_x (str | int): data id or index
        """
        self._sampleList = np.delete(self._sampleList, self.sampleidx(id_x))

    # improvements: args as dictionary
    def get(self, attr, kwargs={}, which=[], includeEmpty=False, includeDiluter=False):
        """Calls method from data objects.

        Args:
            attr:
            args:
            which (list): list of indexes for which the calculation will be done
            includeEmpty (bool): 
            includeDiluter (bool): 
        """
        self.details()
        include = []
        if includeDiluter: include.append('diluter')
        if includeEmpty: include.append('empty')
        if which:
            mask = self._details.index.isin(which) | self._details['type'].isin(include)
        else:
            include += ['ucnp', 'dye']
            mask = self._details['type'].isin(include)

        which = self._details[mask].index
        msgs = []
        for idx in which:
            result = getattr(self.sample(idx), attr)(**kwargs)
            if result[-1] not in msgs:
                msgs += [result[-1]]
        print('\n'.join(msgs))

    def view(self, x='time', yList=['trigger','pm','apd'], which=[], label='', deconv=False, includeEmpty=False, includeDiluter=False):
        '''
        Args:
            which (list): list of integers containing data idx
        '''
        nChannels = len(yList)
        fig01, axs = pl.subplots(nChannels, sharex=True) #, figsize=[10,8])
        
        axs = [axs] if nChannels == 1 else axs
        include = []
        if includeDiluter: include.append('diluter')
        if includeEmpty: include.append('empty')
        if which:
            mask = self._details.index.isin(which) | self._details['type'].isin(include)
        else: 
            include += ['ucnp', 'dye']
            mask = self._details['type'].isin(include)

        which = self._details[mask].index.values
        dataArr = self._sampleList[which]
        
        for dataObj in dataArr:
            dataToPlot = dataObj._dataDeconv if deconv else dataObj._data
            for n in range(nChannels):
                axs[n].plot(dataToPlot[x], dataToPlot[yList[n]],
                    label=f"{dataObj.details()[label]}" if label else f"{dataObj._dataID}",
                    marker=dataObj._plotSettings['m'],
                    ms=4,
                    color=dataObj._plotSettings['colour'],
                    mfc='w',
                    alpha=0.7,
                    linestyle= '' if deconv else '-',
                    lw=1.0
                )

        for n in range(nChannels):
            ylabel = dataObj._metadata.loc[dataObj._metadata['field_name']==yList[n], 'plot_label'].item()
            axs[n].set_ylabel(ylabel)
            axs[n].grid()
        handles, labels = axs[0].get_legend_handles_labels()
        fig01.legend(handles, labels, title=label.replace('_', ' ').strip().title() if label else 'Exp Id',
            loc='center', bbox_to_anchor=(0.88, 0.6),
            ncol=1, fancybox=True, shadow=True, fontsize=12
        )

        xlabel = dataObj._metadata.loc[dataObj._metadata['field_name']==x, 'plot_label'].item()
        axs[-1].set_xlabel(xlabel)

        # axs[0].set_title(f'Measurement ID: {self._dataID}')
        fig01.tight_layout()
        fig01.subplots_adjust(right=0.76)
        fig01.show()

        return(fig01, axs)


class BeamProfile():
    camPxSize = 5.3e-4 # cm (Px size updated on 17/05/2021: previous was 4.3e-4)
    beamWidth = 0 # cm
    _rawImg = {}
    _img = {}

    def __init__(self, imagename, path='../data/beam-profile/'):
        self._rawImg = cv2.imread(path+imagename, cv2.IMREAD_GRAYSCALE)
        self._img = self._rawImg.copy()

    def reset(self):
        self._img = self._rawImg

    def data(self):
        return(self._img)

    def removebackground(self, baseline):
        self._img[np.where(self._img < baseline)] = baseline
        self._img = self._img - baseline

    def removebackground_(self):

        _, row, col = self.centrepeak()
        imgSize = self._img.shape

        samplesize = round(2 * self.beamwidth(pxSize=1))
        cornerRow = 0 if imgSize[0] / 2 - row < 0 else round(imgSize[0] - samplesize)
        cornerCol = 0 if imgSize[1] / 2 - col < 0 else round(imgSize[1] - samplesize)

        sample = self._img[cornerRow : cornerRow + samplesize, cornerCol : cornerCol + samplesize]
        self._img = self._img - sample.mean()
        # self._img[np.where(self._img < 0)] = 0

    def normalise(self):
        self._img = self._img / self._img.sum()

    def show(self):
        cv2.imshow('Beam stop',self._img)
        cv2.waitKey(0)
    
    def centrepeak(self):
        pos = np.where(self._img == self._img.max())
        row, col = pos[0][0], pos[1][0]

        return(self._img.max(), row, col)

    def beamwidth(self, method=None, pxSize=None):

        pxSize = pxSize if pxSize else self.camPxSize
        if (method == 'fwhm') or (method is None):
            coeff = 0.5
        elif method == '1/e2':
            coeff = 0.135
        maxPeak, row, col = self.centrepeak()
        threshold = coeff * maxPeak
        rowWidth = pxSize * len(np.where(self._img[row, :] > threshold)[0])
        colWidth = pxSize * len(np.where(self._img[:, col] > threshold)[0])
        # print(rowWidth, colWidth)
        self.beamWidth = (colWidth + rowWidth) / 2
        return(self.beamWidth)

    def beamarea(self):
        '''Calculates the beam profile area: pi*r^2
        '''
        return(np.pi * (self.beamwidth(method='fwhm') / 2) ** 2)

    def trim(self, ratio=1.3):

        peak, row, col = self.centrepeak()
        beamwidth = round(ratio * self.beamwidth(method='fwhm', pxSize=1))

        self._img = self._img[row - beamwidth : row + beamwidth, col - beamwidth : col + beamwidth]

    def histprofile(self, normalised=False, show=False):

        counts, binEdges = np.histogram(self._img, bins=range(1, round(self._img.max() + 2), 1))
        if show:
            fig, ax = pl.subplots()
            ax.barh(binEdges[:-1], counts)
            ax.set_ylabel('Intensity (counts)')
            ax.set_xlabel('Number of pixels')
            pl.show()

        if normalised:
            binEdges = binEdges / np.sum([c * b for c, b in zip(counts, binEdges[:-1])])
        return(counts, binEdges)

    def plotsurf(self):
        
        (x, y) = np.meshgrid(np.arange(self._img.shape[0]), np.arange(self._img.shape[1]))
        fig = pl.figure()
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(x, y, self._img, cmap=pl.cm.coolwarm)
        fig.colorbar(surf)

        ax.set_xlabel('X (px)')
        ax.set_ylabel('Y (px)')
        ax.set_zlabel('Intensity (counts)')

        pl.show()
    