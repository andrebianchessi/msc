pandoc -o dist/thesis.pdf --from markdown --template eisvogel --listings -F pandoc-crossref --citeproc \
config.md \
introduction.md \
objectives.md \
exampleSection/exampleSection.md \
references.md