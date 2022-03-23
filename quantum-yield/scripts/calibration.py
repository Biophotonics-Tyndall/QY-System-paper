import pandas as pd
import matplotlib.pyplot as pl
from scipy import stats
import numpy as np

   
def polyfit(x, y, degree):
    """
    """
    # Polynomial Coefficients
    coeffs = np.polyfit(x, y, degree).tolist()
    
    # correlation
    r_val = np.corrcoef(x, y)[0,1]

    return(coeffs, r_val)


def calibratePM(laser, pmRange, densityFilter=0, degree=1):
    """
    """
    PMCalibrationFile = '../data/equipment-calibration/calibration_data.csv'
    metadataFile = '../data/equipment-calibration/calibration_data_metadata.csv'
    # loads metadata
    metadata = pd.read_csv(metadataFile)
     # loads calibration data 
    calibrationData = pd.read_csv(PMCalibrationFile)
    filterData = (calibrationData['pm_range'] == pmRange) \
        & (calibrationData['laser'] == laser) \
        & (calibrationData['density_filter'] == densityFilter)
    calibrationData = calibrationData[filterData]

    def getlabel(col, metadata=metadata):
        return(
            metadata.loc[metadata['column'] == col, 'plot_label'].item()
        )


    coeffs, r_val = polyfit(calibrationData['daq_input'], calibrationData['pm_reading'], degree=degree)
    fitting = np.poly1d(coeffs)

    fig, ax = pl.subplots()
    dateColours = ['steelblue', 'cadetblue', 'dodgerblue', 'MediumTurquoise', 'Navy']
    k = 0
    for date in calibrationData['acquisition_date'].unique():
        mask = calibrationData['acquisition_date'] == date
        ax.plot(calibrationData.loc[mask, 'daq_input'], calibrationData.loc[mask, 'pm_reading'],
            label=f'Calibration data: {date}',
            marker='s',
            ms=4,
            color='black',
            mfc=dateColours[k % len(dateColours)],
            alpha=1.0,
            linestyle='',
            lw=1.5
        )
        k+=1
    
    ax.plot(calibrationData['daq_input'], fitting(calibrationData['daq_input']),
        label='Fitting curve',
        marker='',
        ms=4,
        color='firebrick',
        alpha=0.8,
        linestyle='-',
        lw=1.8
    )
    ax.legend()
    ax.set_title("Callibration curve")
    ax.set_ylabel(getlabel('pm_reading'))
    ax.set_xlabel(getlabel('daq_input'))
    
    xList = [f'$x^{n + 1}$' for n in range(len(coeffs) - 1)[::-1]]
    xList.append('')
    yEq = 'y = ' + ' + '.join([f"{c:.3} " + x for x, c in zip(xList, coeffs)])
    ax.text(
        0.1 * calibrationData['daq_input'].max(),
        0.5 * calibrationData['pm_reading'].max(),
        f'{yEq}\n' + r'$R^2$ = ' + f"{r_val**2:.6f}"
    )
    
    fig.tight_layout()
    fig.show()

    # print(f"R^2:{r_val ** 2}")
    return(coeffs, r_val)
    