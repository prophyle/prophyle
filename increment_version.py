#! /usr/bin/env python3

exec(open("prophyle/version.py").read())

numbers=VERSION.split(".")
numbers[-1]=str(int(numbers[-1])+1)

version=".".join(numbers)

with open("prophyle/version.py","w") as f:
	f.write('VERSION="{}"'.format(version))
