<p align="center">
  <a href="http://mathlab.github.io/ITHACA-SEM/" target="_blank" >
    <img alt="ITHACA-SEM" src="./readme/logo_ithaca-sem.png" width="200" />
  </a>
</p>


# ITHACA-SEM - In real Time Highly Advanced Computational Applications with Spectral Element Methods - Reduced Order Models for Nektar++

## Table of contents
* [Description](#description)
* [Dependencies and installation](#dependencies-and-installation)
* [Examples](#examples)
* [How to cite](#how-to-cite)
* [Authors and contributors](#authors-and-contributors)
* [How to contribute](#how-to-contribute)
	* [Submitting a patch](#submitting-a-patch) 
* [License](#license)


## Description
**ITHACA-SEM** is a C++ implementation  of several reduced order modelling techniques. **ITHACA-SEM** is designed to work with [**Nektar++**](www.nektar.info) simulations.
The current Nektar master branch is periodically merged to this code (last time August 2021). A previous version based on [**Nektar++ 4.4.0**] is found in the deprecated folder.

## Dependencies and installation
Please see the installHowTo.txt in the root directory.

Nektar++ and Eigen is provided within **ITHACA-SEM**.

## Examples
Example session files for Nektar++ are provided in the ITHACA_Test_cases subfolder.

## Authors and contributors
**ITHACA-SEM** is currently developed and mantained at [SISSA mathLab](http://mathlab.sissa.it/) by [Dr. Martin Hess](mailto:mhess@sissa.it) under the supervision of [Prof. Gianluigi Rozza](mailto:gianluigi.rozza@sissa.it).

Contact us by email for further information or questions about **ITHACA-SEM**.
 **ITHACA-SEM** is at an early development stage, so contributions improving either the code or the documentation are welcome.

## How to contribute
We'd love to accept your patches and contributions to this project. There are just a few small guidelines you need to follow.

### Submitting a patch

  1. It's generally best to start by opening a new issue describing the bug or
     feature you're intending to fix.  Even if you think it's relatively minor,
     it's helpful to know what people are working on.  Mention in the initial
     issue that you are planning to work on that bug or feature so that it can
     be assigned to you.

  2. Follow the normal process of [forking][] the project, and setup a new
     branch to work in.  It's important that each group of changes be done in
     separate branches in order to ensure that a pull request only includes the
     commits related to that bug or feature.

  3. Do your best to have [well-formed commit messages][] for each change.
     This provides consistency throughout the project, and ensures that commit
     messages are able to be formatted properly by various git tools.

  4. Finally, push the commits to your fork and submit a [pull request][]. Please,
     remember to rebase properly in order to maintain a clean, linear git history.

[forking]: https://help.github.com/articles/fork-a-repo
[well-formed commit messages]: http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html
[pull request]: https://help.github.com/articles/creating-a-pull-request

## License
**ITHACA-SEM** is freely available under the MIT license.
