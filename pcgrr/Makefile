.PHONY: all pkgdown

check:
	@R -e "devtools::check()" --quiet --no-restore --no-save

pkgdown:
	@R -e "pkgdown::build_site()" --quiet --no-restore --no-save

roxydoc:
	@R -e "devtools::document()" --quiet --no-restore --no-save

build:
	@R -e "pak::local_install(upgrade = FALSE, dependencies = FALSE)" --quiet --no-restore --no-save

