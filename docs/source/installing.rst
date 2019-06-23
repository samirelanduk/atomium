Installing
----------

pip
~~~

atomium can be installed using pip:

``$ pip3 install atomium``

atomium is written for Python 3, and does not support Python 2.

If you get permission errors, try using ``sudo``:

``$ sudo pip3 install atomium``


Development
~~~~~~~~~~~

The repository for atomium, containing the most recent iteration, can be
found `here <http://github.com/samirelanduk/atomium/>`_. To clone the
atomium repository directly from there, use:

``$ git clone git://github.com/samirelanduk/atomium.git``


Requirements
~~~~~~~~~~~~

atomium requires `requests <http://docs.python-requests.org/>`_ for fetching
structures from the RCSB, `paramiko <http://www.paramiko.org//>`_ for
fetching structures over SSH,
`msgpack <https://github.com/msgpack/msgpack-python>`_ for parsing .mmtf files,
and `valerius <https://valerius.samireland.com>`_ for dealing with sequences.


Testing
~~~~~~~

To test a local version of atomium, cd to the atomium directory and run:

``$ python -m unittest discover tests``

You can opt to only run unit tests or integration tests:

``$ python -m unittest discover tests.unit``
``$ python -m unittest discover tests.integration``

You can run the 'big test' to get a random 1000 structures, parse them all, and
report any problems:

``$ python tests/big.py``

Finally, to perform speed profiles you can run:

``$ python tests/time/time.py``

...which creates various profiles that SnakeViz can visualise.

