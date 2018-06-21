
Welcome to meteo_si's documentation!
====================================

`Meteo SI` is a collection of routines for atmospheric sciences. Unless noted
 otherwise, all units are SI. 


Download
--------

The code is available at https://github.com/maahn/meteo_si

Installation
------------

Change to the folder containing the project and do ::

  python setup.py install

in the terminal. If you do not have root privileges, you can also do ::

  python setup.py install --user

which will install meteo_si in userbase/lib/pythonX.Y/site-packages or ::

  python setup.py install --home=~

which will install meteo_si in ~/lib/python.



Density
-------
Routines to estimate the density of air:

.. automodule:: meteo_si.density

Humidity
--------
Routines to work with different units of humidity:

.. automodule:: meteo_si.humidity

Temperature
-----------
Routines to convert from and to virtual temperature and Celsius.

.. automodule:: meteo_si.temperature

Wind
----
Functions to work with wind observations. 

.. automodule:: meteo_si.wind

Indices and tables
==================
 
* :ref:`genindex`
* :ref:`search`


