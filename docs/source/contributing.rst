Contributing to atomium
=======================

Thank you for considering spending your time on improving this project!
You can help by raising an issue to report a bug or suggest a new
feature, or by creating a pull request to add a new feature yourself
(subject to approval).

Raising an Issue
----------------

The `GitHub issue
tracker <https://github.com/samirelanduk/atomium/issues>`__ is the place
to either report a bug, or propose a new feature that you'd like to see
added. There will be a style guide for the text of the issue when you
open the tracker - just delete the feature request section if you are
reporting a bug, and *vice versa*.

Pull Requests
-------------

Pull requests are GitHub's mechanism for allowing willing contributors
to add features or make fixes themselves, without waiting for me to get
around to it.

How do Pull Requests work?
~~~~~~~~~~~~~~~~~~~~~~~~~~

The first step is to fork the atomium repository to your own GitHub
account. This will give you your own copy of the repository, where you
can make changes that don't affect the original.

The second step is to make your desired changes. To do this you need to
clone your new copy of atomium to your local machine, and to do this
you simply enter ``clone https://github.com/YOUR_USERNAME/atomium`` at
your terminal, and you will have a local copy to work with.

Make commits at sensible points while making your changes, with
informative commit messages in the present tense.

Once you are finished, push your changes to your forked repository on
GitHub. You are now ready to make the pull request.

Find the 'New Pull Request' button and click it, and you will be given
an overview of the pull request. The default branch to merge into is
master, which is fine - I will manually change this to whatever the
current branch being worked on is at the time. That way your feature (if
accepted) will appear in the next release (you will be credited).

Look over the changes one last time, and if it's all good, click 'create
pull request'. It then gets sent to me to look over, and either accept
and merge, close, or request changes.


What does and doesn't make a good Pull Request?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An ideal pull request for atomium is one that either adds a function (or
set of functions which carry out a single piece of functionality), or
modifies the minimal amount of code to fix one bug. One pull request,
one feature.

If your pull request modifies large parts of atomium, changes the way the
core API works for the user, or adds an entirely new class, there's a
good chance it won't be merged.

If the pull request resolves a current, open issue, there's a good
chance it will be accepted. If in doubt, asking beforehand is always a
good idea!

Requirements
~~~~~~~~~~~~

Style
^^^^^

-  Two lines between functions, three lines between classes.

-  Lines strictly no more than 80 characters long. This is less strict
   in test files but still try to keep them below 80 characters.

-  underscore\_naming\_convention

-  Docstrings do not have a line break between the triple quote marks
   and the start of the string.

Documentation
^^^^^^^^^^^^^

If a function or class doesn't have documentation, it doesn't exist. Not
to the user anyway. All functions should begin with a docstring which
summarises what the function does, lists the parameters it takes,
outlines any exceptions it can raise, and specifies what the function
returns.

The text should be in RST format. For example:

.. code::

    def calculate_hypotenuse(side1, side2):
        """Takes the lengths of two sides of a right-angled triangle, and uses the
        Pythagorean theorem to calculate the hypotenuse length.

        :param float side1: The length of the first side.
        :param float side2: The length of the second side.
        :raises ValueError: if either side is negative.
        :rtype: ``float``"""

        if side1 < 0 or side2 < 0:
            raise ValueError("Sides can't be negative")
        return math.sqrt((side1 ** 2) + (side2 ** 2))

If in doubt, just look at the other functions.

Tests
^^^^^

If a function doesn't have tests, the function doesn't work. atomium has
integration tests and unit tests.

Suppose the hypotenuse function above already exists, and you want to
add a function for getting the distance between two points, using this
function. You might add the following:

.. code::

    def distance_between(point1, point2):
    	"""Takes two (x, y) points and returns the distance between them.

    	:param list point1: The first point in the form ``[x, y]``.
    	:param list point2: The first point in the form ``[x, y]``.
    	:rtype: ``float``"""

    	delta_x = point2[0] - point1[0]
    	delta_y = point2[1] - point1[1]
    	return calculate_hypotenuse(delta_x, delta_y)

How should this be tested?

Unit Tests
''''''''''

Unit tests test a function *in isolation*. In this case, the unit test
would check that the function works but it should not execute
``calculate_hypoteneuse``! The test might look like this:

.. code::

    from unittest import TestCase
    from unittest.mock import patch

    class DistanceTests(TestCase):

        @patch("calculate_hypotenuse")
        def test_can_get_distance_between_points(self, mock_hyp):
            point1 = [0, 0]
            point2 = [4, 3]
            mock_hyp.return_value = 5
            distance = distance_between(point1, point2)
            mock_hyp.assert_called_with(point1, point2)
            self.assertEqual(distance, 5)

The ``calculate_hypoteneuse`` function is patched with a mock object
here. We set its return value and just ensure that it was called, and
that what it returns is what our function returns.

Note that this way if ``calculate_hypoteneuse`` is broken, the tests for
``distance_between`` will still pass - they are isolated.

Unit tests live in ``tests/unit``. Each class/collection of functions
gets its own test file, each function gets its own test class, with
different test functions for each possible use of the function.

Again, see existing tests for numerous examples.

Integration Tests
'''''''''''''''''

Integration tests check that the code works when called as the user
would call it. Nothing is mocked or patched - this is a test that all
the functions work together to do what the user wants.

If your pull request is to add a function that works 'under the hood'
and which the user never uses, you don't need to add an integration test
(the existing tests will cover it). If you've added user-facing code, it
does need a few lines. Just fine somewhere suitable in one of the
``tests/integration`` files and add it in - don't worry too much about
putting it in the right place as I move things around pretty often
anyway.

So in this case, you might just add the line:

.. code::

    self.assertEqual(distance_between([0, 0], [3, 4]), 5)


Final Checks
^^^^^^^^^^^^

All tests should be run before submitting the pull request.

Unit tests are run with:

.. code::

    $ python -m unittest discover tests/unit


Integration tests are run with:

.. code::

    $ python -m unittest discover tests/integration
