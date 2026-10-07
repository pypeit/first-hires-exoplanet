.. _installation:

Installation
============

To install **pyodine**, the only way right now is to simply download or clone it from its Gitlab repository:

.. code-block:: console

   $ git clone https://gitlab.com/Heeren/pyodine.git

Now you can use it by importing it into your own Python code where needed (check out the *Tutorial* section on how to import and use it: :doc:`tutorial/preparation`).

Of course you should make sure that you have all dependencies installed - these will only be Python packages:

.. IMPORTANT::
   List of required **packages**

    * astropy >= 4.2.1
    * h5py >= 2.10.0
    * numpy >= 1.20.1
    * lmfit >= 1.0.2
    * barycorrpy >= 0.4.4
    * pathos >= 0.2.8
    * argparse >= 1.1
    * matplotlib >= 3.3.4
    * progressbar2 >= 3.37.1
    * dill >= 0.3.4

Also we want to warn you that we have noticed small differences in the results between using different versions of Python 3 - specifically, between Python 3.8 and 3.9. At the moment, we do not know the exact cause(s), and have not been able to assess yet whether the results from one of these versions are better (i.e. closer to the truth). In any case, the differences that we have observed so far are small enough (on the cm/s level in the final RV estimates) that all physical interpretations of the results should not be influenced.
