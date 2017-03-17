==================
Configuration File
==================

All commands may be presented in a configuration file fed to bipype with @ prefix.
A file with configuration should contain all desired commands and their options (if applicable) one per line.

Example 1
---------

``bipype @fancy_filename``

Content of file with fancy filename:

.. code-block:: none

   -e
   --mode run
   --threads 4
   --output_type 'ITS'

So, this is equivalent to

``bipype -e --mode run --threads 4 --output_type 'ITS'``

Example 2
---------

``bipype @another_fancy_filename``

Content of file with another fancy filename

.. code-block:: none

   --mode run
   --threads 4

So, this is equivalent to

``bipype -e --mode run --threads 4 --output_type 'ITS'``
