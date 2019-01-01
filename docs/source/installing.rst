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
