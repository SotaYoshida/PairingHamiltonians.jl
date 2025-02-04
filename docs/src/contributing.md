# Contributing to PairingHamiltonians.jl

Thank you for considering contributing to this package.

Feedbacks and contributions to this package are very welcome.
These can be:
- bug report
- submitting a new function or a patch to the bug
- documentation issue
- feature request
- etc.

For these contributions, it would be nice to let you know a basic policy (the workflow of code development and LICENSE, shown below) in this package.
Note that the package has been developed by a single author (@SotaYoshida) so far, and thereby the followings are just the author's policy. Comments on the development policy itself are also welcome.


## Workfklow of code development

We use the GitHub to host the package, to track issues/pull requests.

The document is built using [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/),
which is a package for building documentation from docstrings and markdown files.
It is automized to build and deploy the document using GitHub Actions.

A typical workfklow of code development is the following:

1. clone or folk the repository and modify the code and/or document locally

2. propose modifications through a pull request (PR) to **`develop` branch**

3. if the PR passes the automated tests, it will be merged 

4. At some point, we tag and release, which is done automatically by instructing the JuliaRegistrator bot to do so.


### Automated tests on GitHub Actions

We use GitHub Actions to run the test codes and to build/deploy the document.
When some changes are submitted through a pull request, the test codes are run to check that the changes are not destructive.

The test jobs are specified in yml files like `.github/workflows/CI.yml` and one can find the test code in `test/` of the repository.

!!! note
    If you submit a major change to the code, please consider to update the test codes matching to your modifications.

## Reporting bugs by opening a new issue

Thank you so much for considering to report bugs!
When you report a bug in the code, please open an new issue at the repository.
Please make sure to include any information necessary for us to reproduce the errors. Thanks!

## Propose modifications through Pull Requests (PRs)

You can propose modifications you made on the package through pull requests.

* As stated above, please consider to make test codes matching to your modifications.
 
* Please make sure to submit your PR to `develop` branch. The 'main' branch will be protected by 'github branch protection'.

As this package is currently being developed by a single author, branching rules such as git-flow and GitHub Flow have not been adopted. When we got contributors, the branch-rule will be set upon discussions.

## LICENSE

Any contribution from you will be under the [MIT License](https://opensource.org/licenses/MIT), as well as the package itself.
Feel free to contact to @SotaYoshida if that's a concern.