=====================
Writing Pythonic Code
=====================

Code is often described as "Pythonic" or "not Pythonic" (with the
implication that "Pythonic" is better.) The description is often
applied to refer to code that reflects best practices which have
emerged from the Python community and have become almost second nature
to experienced programmers.

Reading lots of Python code (particularly from well-respected long
maintained community projects) is the best way to develop this sense,
but some principles are described here.

Good coding practice
====================

Familiarity with modern coding practices that apply across most
languages is a good start:

    * Compact but descriptive names for variables, functions, etc.
    * Succinct comments where necessary
    * Encapsulation of data and abstraction through functions and classes
    * Use of existing libraries rather than reimplementing

Using language features
=======================

Where Python or its standard library provides a means of accomplishing
a task, it is generally preferred to use that means rather than
reimplementing the wheel. The canonical example is using list
comprehensions rather than for loops to transform a list::

  >>> newlist = [i + 1 for i in oldlist]

not::

  >>> newlist = []
  >>> for i in range(len(oldlist)):
  ...     newlist.append(oldlist[i] + 1)


Note there are several non-Pythonic ways to perform this task.

For those not familiar with the features of the lanaguage and the
standard library, this does represent a barrier to entry. However once
that knowledge is built, using the features of the language makes
one's intention much more clear. It often also results in shorter code
that is easier to comprehend.

See several examples in :doc:`tips`.


Idiom and communication
=======================

Because "code is more often read than written," anything that improves
clarity is beneficial. A list comprehension and a for loop may have
the same result, but the use of a list comprehension immediately makes
it apparent to the reader that the code is intended to create a new
list based on some element-by-element translation of the input
list. It is a pattern with a common solution, and sticking to the
common solution helps make the pattern apparent so the reader of the
code understands the underlying problem.

Generally this choice of the common way is referred to as "idiomatic
Python." This can be expanded to conventions such as the use of "self"
as the first argument in instance methods, even though such choice is
generally free.


Further Reading
===============

A web search for "pythonic" will give a wealth of opinions. These references are a good starting point.

    * `Python Style Guide (PEP 8) <https://www.python.org/dev/peps/pep-0008/>`_
    * `Zen of Python (PEP 20) <https://www.python.org/dev/peps/pep-0020/>`_
    * `What is Pythonic? <https://blog.startifact.com/posts/older/what-is-pythonic.html>`_
    * `Examples of idiomatic and nonidiomatic Python <https://medium.com/the-andela-way/idiomatic-python-coding-the-smart-way-cc560fa5f1d6>`_
    * `Idomatic Python from Wikibooks <https://en.wikibooks.org/wiki/Python_Programming/Idioms>`_


--------------------------

:Release: |version|
:Doc generation date: |today|
