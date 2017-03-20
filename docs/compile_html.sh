#!/usr/bin/env sh

pandoc --from=markdown --to=rst --output=annotation_resources.rst annotation_resources.md
pandoc --from=markdown --to=rst --output=getting_started.rst getting_started.md
pandoc --from=markdown --to=rst --output=about.rst about.md
pandoc --from=markdown --to=rst --output=output.rst output.md
#pandoc --from=markdown --to=rst --output=index.rst index.md


make html

