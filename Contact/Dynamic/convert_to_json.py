import argparse
import io
import yaml
import json

parser = argparse.ArgumentParser(description="convert config files to json format")
parser.add_argument("--input")
parser.add_argument("--output")

args = parser.parse_args()

filename = args.input
filename_out = args.output

print("input = ", args.input)
print("output = ", args.output)

libconf_found = False
try:
  import libconf
  libconf_found = True
except:
  print("libconf not found!")  

fh = io.open(filename, "r")

if filename.endswith(".yaml"):
  settings = yaml.load(fh)
elif filename.endswith(".cfg") and libconf_found:
  settings = libconf.load(fh)

fh_out = io.open(filename_out, "w+")
json.dump(settings, fh_out)
fh_out.close()

fh.close()

