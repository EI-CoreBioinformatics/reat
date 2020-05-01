# Cromwell configuration files

The files in this directory are example configuration files for running cromwell under different conditions.

* cromwell_noserver_slurm.conf - SLURM configuration using cromwell's **run** mode, there is no cromwell server and the database is backed on disk.
* cromwell_noserver.conf - Localhost configuration using cromwell's **run** mode, no cromwell server configuration and disk backed database.
* cromwell_server_options.conf - Localhost configuration with mysql database for call caching configured also on the localhost, this configuration can be used with a cromwell **server** instead of just the run mode.