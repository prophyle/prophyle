#! /usr/bin/env python3

import re
#from datetime import datetime
import datetime
import sys
from io import StringIO

#r=re.compile(r"(\d+\.\d+)user (\d+\.\d+)system (\d+):(\d+\.\d+)elapsed (\d+)%CPU \(\d+avgtext\+\d+avgdata (\d+)maxresident\)")
r=re.compile(r"([0-9.]+)user ([0-9.]+)system ([0-9.:]+)elapsed (\d+)%CPU \(\d+avgtext\+\d+avgdata (\d+)maxresident\)")


def rm_ms(delta):
    return delta - datetime.timedelta(microseconds=delta.microseconds)

last_section = None


old_stdout = sys.stdout
sys.stdout = mystdout = StringIO()

toc=[]

for line in sys.stdin:

	#s="1048.49user 6.51system 18:00.09elapsed 97%CPU (0avgtext+0avgdata 6269888maxresident)"
	
	a=r.match(line)
	#print(a)
	if a is not None:
		g=a.groups()
		#print(a)
		#print(g)

		user_time=datetime.timedelta(seconds=float(g[0]))
		system_time=datetime.timedelta(seconds=float(g[1]))

		el_time_s=0

		for p in g[2].split(":"):
			el_time_s=el_time_s*60 + float(p)

		elapsed_time=datetime.timedelta(seconds=el_time_s)

		#user_time=datetime.strptime(g[0],"%S")
		#elapsed_time=datetime.strptime(g[2],"%M:%S.")

		memory_kb=float(g[-1])/4
		memory_gb=memory_kb/(1024*1024)

		print()
		print("```")
		print("User time:     {}".format(rm_ms(user_time)))
		print("System time:   {}".format(rm_ms(system_time)))
		print("Elapsed time:  {}".format(rm_ms(elapsed_time)))
		print("CPU usage:     {:7}%".format(round(100*(user_time.total_seconds() + system_time.total_seconds()) / elapsed_time.total_seconds())))
		print()
		print("Memory peak:   {:7.2f} GB".format(round(memory_gb,2)))
		print("```")
		print()
	else:
		line=line.replace(" <==","").replace("==>","##")
		if line[0]=="#":

			parts=line.split("/")
			name=parts[0].replace("## ","").strip()
			full_name=line.replace("## ","").strip()

			if last_section!=parts[0] and len(parts)>1:
				last_section=parts[0]
				print("***")
				print(last_section)
				toc.append("* [{}](#{})".format(name,name))
			toc.append("  * [{}](#{})".format(full_name,full_name.replace("/","").replace(".","").lower()))
			print()
			print("#"+line,end="")
		elif line[0].strip()=="":
			print()
		else:
			print("* "+line,end="")

sys.stdout = old_stdout

print("Table of Contents")
print("=================")
print()
print("\n".join(toc))
print()
print(mystdout.getvalue())
