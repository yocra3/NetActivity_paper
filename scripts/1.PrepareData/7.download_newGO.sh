## Download GO ontology
wget http://current.geneontology.org/ontology/go-basic.obo -P data/GO_terms/
wget https://raw.githubusercontent.com/sgrote/OboToTerm/master/obo_to_term_tables.py  -P data/GO_terms/
wget https://raw.githubusercontent.com/sgrote/OboToTerm/master/obo_to_term_functions.py  -P data/GO_terms/


python3 data/GO_terms/obo_to_term_tables.py data/GO_terms/go-basic.obo data/GO_terms/


## Remove links differents to is_a
grep -v ^relationship data/GO_terms/go-basic.obo > data/GO_terms/go-basic-mod.obo

## Remove multicellular organismal process and rhythmic process -> reduces distance for too many GOs
grep -v 'is_a:.*GO:0032501' data/GO_terms/go-basic-mod.obo | grep -v 'is_a:.*GO:0048511' > data/GO_terms/go-basic-mod2.obo
mkdir data/GO_terms/mod_graph/
python3 data/GO_terms/obo_to_term_tables.py data/GO_terms/go-basic-mod2.obo data/GO_terms/mod_graph
