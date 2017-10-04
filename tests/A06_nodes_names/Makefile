.PHONY: all clean index readcount nodenames

include ../conf.mk

K=10

all: readcount nodenames

readcount: index
	La=$$(cat $(FQ) | grep "^@" | wc -l); \
	Lb=$$(cat _prophyle_index_output.txt | awk '{print $$2}' | sort | uniq | wc -l); \
	if ! [[ "$${La}" -eq "$${Lb}" ]]; then \
		echo "Files have a different read count"; \
		exit 1; \
	fi;

nodenames: index
	cat _prophyle_index_output.txt | ./extract_reported_contigs.py | sort | uniq | grep -E -v "^[0A]$$" > _nodenames.reported.txt
	cat index.fa | grep -e "^>" | perl -pe "s/>(.*)@.*/\1/g" | sort | uniq > _nodenames.present.txt
	cat _nodenames.present.txt _nodenames.reported.txt | sort | uniq > _nodenames.present_reported.txt
	diff -c _nodenames.present.txt _nodenames.present_reported.txt

index:
	$(1STEP) -k $(K) $(FA) $(FQ) > _prophyle_index_output.txt

clean:
	rm -f _* index.fa.*
