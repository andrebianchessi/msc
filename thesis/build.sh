pandoc -F pandoc-crossref --citeproc \
config.md \
introduction.md \
objectives.md \
exampleSection/exampleSection.md \
references.md \
-o dist/main.pdf