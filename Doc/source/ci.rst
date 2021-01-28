======================
Continuous Integration
======================

SpacePy uses `GitHub Actions <https://docs.github.com/en/actions>`_
for continuous integration. Most of the relevant information is
checked into the repository: the configuration file `CI.yml
<https://github.com/spacepy/spacepy/blob/master/.github/workflows/ci.yml>`_
manages the CI process, which ultimately runs `the unit tests
<https://github.com/spacepy/spacepy/blob/master/tests/test_all.py>`_. However
a few elements of the setup are not in the repository and are
documented here. This may be useful if this ever has to be set up in
the future, or if you want to run SpacePy CI tests on your fork before
opening a pull request.

.. contents::
   :local:

Initial run
===========

A workflow cannot be run `until it has been run once against the
default branch <https://github.community/t/
workflow-dispatch-event-not-working/128856/2>`_ (``master``). This makes
it somewhat hard to test the workflow before merging; in SpacePy this was
handled by merging a tiny workflow first (PR `496 <https://github.com/
spacepy/spacepy/pull/496>`_).
      
Once this first run has been made, updated versions of the workflow
can be run from topic branches. It will show up under the ``Actions``
tab of a repository and the branch to use can be selected. It is also
possible to specify a branch `by using the REST API <https://
github.community/t/workflow-dispatch-workflow-not-showing-in-actions-tab/
130088/15>`_; you will need `an access token <https://docs.github.com/en
github/authenticating-to-github/creating-a-personal-access-token>`_ with
``workflow`` scope.

Merging rules
=============

PRs require CI to pass before merging; this is managed with a `branch
protection rule <https://docs.github.com/en/github/
administering-a-repository/managing-a-branch-protection-rule>`_
(``Settings`` from the tab at the top of a repository, ``Branches`` from
the left menu.) The relevant choices is "Require status checks to pass
before merging." Every variant of the unit testing job (``test (2.7,
ubuntu-18.04...`` etc.) will be in the list of checks; leave these alone and
select only the ``All tests`` job. The name of this won't change and it
will always depend on *all* the jobs in the workflow.
"Require branches to be up to date" should *not* be selected;
this encourages merging rather than our preferred rebase, and the tests
will run against a (temporary) merge regardless.

Rerunning CI on a PR
====================

There is no way to manually trigger a workflow run on a pull request.
SpacePy's CI workflow is set up to `trigger the workflow <https://
docs.github.com/en/actions/reference/events-that-trigger-workflows
#pull_request>`_ when a PR is marked ready for review, so one way to
force a run is to mark the PR as a draft, and then as ready again.

Note that a PR will not trigger the CI `if there is a merge conflict
<https://github.community/t/run-actions-on-pull-requests-with-merge-conflicts/
17104>`_.

Usage
=====

SpacePy administrators can view the usage minutes, storage (for caches),
etc. on our `billing page <https://github.com/organizations/spacepy/settings/
billing>`_.

--------------------------

:Release: |version|
:Doc generation date: |today|
