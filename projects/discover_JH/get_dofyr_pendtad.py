import numpy as np

def get_dofyr_pentad(dates):
    """
    This function converts dates in YYYYMMDD format to decimal year and pentad (5-day) of year.

    Parameters:
        dates (ndarray): an array of dates in YYYYMMDD format

    Returns:
        dofyrs (ndarray): an array of decimal year and pentad (5-day) of year
    """

    # Convert dates to decimal year
    year = np.floor(dates / 10000)
    month = np.floor((dates - year * 10000) / 100)
    day = dates - year * 10000 - month * 100
    doy = np.floor((275 * month) / 9) - ((month + 9) / 12) * 275 + day - 30
    leap_years = np.floor((year - 2000) / 4) - np.floor((year - 2000) / 100) + np.floor((year - 1600) / 400)
    doy = doy + leap_years
    days_in_year = 365 + leap_years
    dofyrs = year + doy / days_in_year

    # Convert decimal year to pentad of year
    pentad = np.floor((doy - 1) / 5) + 1
    dofyrs = np.concatenate((dofyrs[:, np.newaxis], pentad[:, np.newaxis]), axis=1)

    return dofyrs
