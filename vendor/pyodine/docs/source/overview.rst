.. _overview:

Overview of main routines
=========================

The main routines for template creation, observation modelling and combination of chunk velocities all reside in the top level of the **pyodine** repository. We have built them with high flexibility between instruments in mind, meaning that no matter for which instrument you want to run the software, you always use the same main routines - all instrument-specific code and parameters are defined separately in the ``utilities`` modules (see :ref:`overview_utilities`).

Create a deconvolved stellar template
-------------------------------------

Before you can model observations of a star, you need to create a deconvolved stellar template, i.e. a high-resolution spectrum of the same star **WITHOUT** any I2 features and cleaned of the instrumental line-spread function (LSF). This is done using the :func:`pyodine_create_templates.create_template` function:

.. autofunction:: pyodine_create_templates.create_template


Model a single observation
--------------------------

After you have created a deconvolved stellar template, you can use it to model an observation of the same star, obtained **WITH** the I2 cell in the light path. This is done with the function :func:`pyodine_model_observations.model_single_observation`:

.. autofunction:: pyodine_model_observations.model_single_observation


Model multiple observations
---------------------------

In many cases you will have more than one single observation of the same star that you want to model (or even observations of several different stars, each of which you already have a deconvolved stellar template for). Then we can take advantage of Python's parallelizing capabilities and use the function :func:`pyodine_model_observations.model_multi_observations` to model the observations on several cores at the same time:

.. autofunction:: pyodine_model_observations.model_multi_observations


Combine chunk velocities to RV timeseries
-----------------------------------------

When you've modelled a number of observations and saved the fit results, you can use the function :func:`pyodine_combine_vels.combine_velocity_results` to combine the chunk velocities from the observations of any one star to a RV timeseries:

.. autofunction:: pyodine_combine_vels.combine_velocity_results
