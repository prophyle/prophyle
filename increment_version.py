#! /usr/bin/env python3

vfn="prophyle/version.py"

exec(open(vfn).read())

numbers=VERSION.split(".")
numbers[-1]=str(int(numbers[-1])+1)

version=".".join(numbers)

with open(vfn,"w") as f:
	f.write('VERSION="{}"'.format(version))
