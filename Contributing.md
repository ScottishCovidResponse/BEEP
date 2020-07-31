# Contributing to BEEPmbp

At present, contributions are limited to members of Scottish COVID-19
Response Consortium (SCRC).  If you are a member, please read the
below on how best to make progress on this project.

## Pull Request Process

1. Please branch from the dev branch within
   ScottishCovidResponse/BEEPmbp

   - Please **don't** fork this repository to your own
     github.com/username or other organisation.

2. Read the issues and choose/be assigned an appropriate issue. This
   is likely to be either a `Starter issue` or one assigned by
   chatting in `zulip/BEEPmbp`. You could also create an issue if
   there is something you want to work on that has not been mentioned.

3. Name the branch `username/featurename`

4. If referencing the issue from within a commit message then it
   should have the correct path to the issue on the
   `ScottishCovidResponse/SCRCIssueTracking` repository. An autolink
   reference has been created to make a shortcut;
   `ScottishCovidResponse/SCRCIssueTracking#123` -> `SCRC-123`. Please
   put the reference in the commit description (not the first line of
   the commit) i.e.:

   - `git commit -m "pithy description" -m "SCRC-123 and more details"` to refer to issue number 123. 

    - writing the commit message and description in an editor is also
      possible with `export EDITOR=vi` for example, and then just
      doing a `git commit`.

5. Please run the tests (see [README.md](README.md)) before pushing.

6. When your feature is ready to merge into the `dev` branch
   create a pull request against `ScottishCovidResponse/BEEPmbp#dev`
   and assign `ScottishCovidResponse/beepmbp` as reviewers.

7. Semver will be handled in PRs from `dev` into `master`.

8. Thanks for your time and effort!

## Coding style

- The code uses hard tabs of width 2. Please configure your editor to
  use this.

- In general, follow the style of the code surrounding the code you
  are writing.

- Code documentation is generated using Doxygen.  Please add doxygen
  comments to your code if possible.

## Licensing

By contributing to this project (e.g. by merging, submitting a pull
request or providing advice on code), you agree - unless
simultaneously and expressly stated otherwise - that your contribution
may be included in the source code of the project and published under
the 2-Clause BSD License and that the contribution was created in
whole or in part by you and you have the right to submit it under the
open source license indicated above.
