#!/usr/bin/env python3

# import necessary modules
import sys
from collections import defaultdict
import pickle

ReferenceSegmentsPickle,LinksPickle,UpstreamOutputFile,DownstreamOutputFile = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]

def rec_dd():
    return defaultdict(rec_dd)

ref_dict = {}
link_dict = {}

print("Adding links segments...")
with open(LinksPickle, 'rb') as f:
    link_dict = pickle.load(f)
f.close()

print("Adding reference segments...")
with open(ReferenceSegmentsPickle, 'rb') as f:
    ref_dict = pickle.load(f)
f.close()


# Build arrays of upstream and downstream keys from link_dict
print("Building arrays of upstream and downstream keys from link_dict...")
upstreamkeyarray = []
downstreamkeyarray = []

for key, value in link_dict.items():
    # Filter out links that aren't reference segments to reference segments
    # print(ref_dict[value["RefSegmentID"]])
    # print(value["RefSegmentID"])
    # print(type(value["RefSegmentID"]))
    # print(key)
    # print(type(key))

    if ref_dict[value["RefSegmentID"]] and ref_dict[key]:
        downstreamkeyarray.append(key)
        upstreamkeyarray.append(value["RefSegmentID"])

    #31611814
    #31611815
    # if value["RefSegmentID"] == "31611815" and key == "31611814":
    #     #ref_dict[value["RefSegmentID"]] and ref_dict[key]:
    #     print("It matches")
    #     print(ref_dict[value["RefSegmentID"]])
    #     print(value["RefSegmentID"])
    #     print(type(value["RefSegmentID"]))
    #     print(repr(value["RefSegmentID"]))
    #     print(key)
    #     print(type(key))
    #     print(ref_dict[key])
        #print(upstreamkeyarray.index('31611815'))


# Save nested dictionary as pickle file
print("Saving to: ")

print(UpstreamOutputFile)
f = open(UpstreamOutputFile, "wb")
pickle.dump(upstreamkeyarray, f)
f.close()

print(DownstreamOutputFile)
f = open(DownstreamOutputFile, "wb")
pickle.dump(downstreamkeyarray, f)
f.close()
