# Data folder

All acquired data is save in the [raw-data](./raw-data) folder with name in the following format: **_qy_date_timeID.csv_**. All the details of the experiments are logged into the [datalogs](datalogs.csv) file and one can relate the specific data file to its log details by file name and _saved_name_ column. The [metadata](metadata.csv) file contains all the details regarding the columns of the saved data. 

|field_name|data_type|data_format                        |example      |standard_units|description                                                           |
|----------|---------|-----------------------------------|-------------|--------------|----------------------------------------------------------------------|
|time      |float    |^([0-9]?\.&#124;[1-9][\d]*\.)[0-9]*$    |123423.982734|seconds       |Positive float number for time acquired with python using the OS clock|
|0         |float    |^[-]?([0-9]?\.&#124;[1-9][\d]*\.)[0-9]*$|-1.99834     |volts         |Daq input channel                                                     |
|1         |float    |^[-]?([0-9]?\.&#124;[1-9][\d]*\.)[0-9]*$|-1.99834     |volts         |Daq input channel                                                     |
|2         |float    |^[-]?([0-9]?\.&#124;[1-9][\d]*\.)[0-9]*$|-1.99834     |volts         |Daq input channel                                                     |

Beam profile folder and equipment calibration folders contain images of the beam profile and csv data file used to calibrate the analog output of the power meter. 

