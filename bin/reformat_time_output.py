#! /usr/bin/env python3

import re
#from datetime import datetime
import datetime

r=re.compile(r"(\d+\.\d+)user (\d+\.\d+)system (\d+:\d+\.\d+)elapsed (\d+)%CPU \(\d+avgtext\+\d+avgdata (\d+)maxresident\)")

s="1048.49user 6.51system 18:00.09elapsed 97%CPU (0avgtext+0avgdata 6269888maxresident)"

a=r.match(s)
g=a.groups()
print(a)
print(g)

user_time=datetime.timedelta(seconds=float(g[0]))
system_time=datetime.timedelta(seconds=float(g[1]))

#user_time=datetime.strptime(g[0],"%S")
#elapsed_time=datetime.strptime(g[2],"%M:%S.")

memory_kb=float(g[-1])/4
memory_gb=memory_kb/(1024*1024)

print("User time:   {}".format(user_time))
print("System time: {}".format(system_time))
print("Memory peak: {} GB".format(round(memory_gb,2)))
