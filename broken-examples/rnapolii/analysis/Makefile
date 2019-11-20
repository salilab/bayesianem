step=2
N=1000

DIR := $(shell dirname ${MAKEFILE_LIST})

all: table_A-Z.csv

.SECONDARY:

table_%.csv:${DIR}/make_table_new.py ../data/native.rmf3 data_%.pkl clusters_%.pkl
	${IMP} $^ ${N} 3 | tee $@

clusters_%.pkl: data_%.pkl
	${IMP} ${DIR}/clustering_new.py $< $@ ${N}

data_%.pkl: [%]/_output*_${step}/rmfs/*.rmf3
	${IMP} ${DIR}/get_data.py $* $@ ${N} $(shell ${DIR}/get_stat_or_rmf3.sh $^)

clean:
	${RM} table_*.csv clusters_*.pkl data_*.pkl

IMP=~/imp-fast/setup_environment.sh python

.PHONY: all clean
