#!/usr/bin/env python3

## Author: Steffen Heyne, MPI-IE Freiburg, Germany
## script that parses a local directory structure and maps this to a UCSC trackDB.txt config
## https://github.com/steffenheyne/UCSC_trackHub_generator

import argparse
import os.path
import glob
import pprint
import functools
import re
import yaml
import sys

args = None
trackCounter = 1

pp = pprint.PrettyPrinter()

## bigwig colors
## keys are used as regex pattern against filename/ 'track' key in dict
bigwig_colors = {
    "input":    "150,150,150",
    "H3K4me1":  "65,171,93",
    "H3K4me2":  "161,217,155",
    "H3K27ac":  "252,78,42",
    "H3K4me3":  "203,24,29",
    "H3K36me3": "254,196,79",
    "H3K27me3": "140,107,177",
    "H3K9me3":  "29,145,192",
    "H3K9ac":   "164,0,0", # 252,146,114",
    "CTCF":     "106,81,163",
    "WGBS":     "0,102,255",
    "methyl":   "0,102,255",
    "RNA.*fwd": "0,102,0",
    "RNA.*rev": "153,51,0",
    "RNA.*RPKM":"71,107,107",
    "DNase":    "0,204,102",
    "Hp1a":     "0,128,255",
    "H1":       "255,102,255",
    "H3K9me2":  "51,51,255"
    }


multiwig_default = {"track": None,
                    "type": "bigWig",
                    "container": "multiWig",
                    "parent": None,
                    "shortLabel":None,
                    "longLabel": None,
                    "aggregate": "transparentOverlay",
                    "showSubtrackColorOnUi": "on",
                    "priority": 1,
                    "html": "examplePage"
}


bigwig_default = {"track": None,
                  "type": "bigWig",
                  "parent": None,
                  "bigDataUrl": None,
                  "shortLabel":None,
                  "longLabel": None,
                  "color": "255,0,0"
                  }

## bigwig configuration that can be part of multiwig container or 
## individual bigwig track (if track is not in multiwig container) 
bigwig_combined = {"visibility": "hide",
                  "maxHeightPixels": "500:20:8",
                  "viewLimits": "0:20",
                  "alwaysZero": "on",
                  "autoScale": "off",
                  "windowingFunction": "mean+whiskers",
                  "priority": 1
                  }

bigbed_default = {
    "track": None,
    "parent": None,
    "bigDataUrl": None,
    "shortLabel":None,
    "longLabel": None,
    "type": "bigBed 3 +",
    "color": "255,0,0",
    "visibility": "squish",
    #"colorByStrand": "255,0,0 0,0,255"
    }

## specific bigwig configurations
## either for multiwig container or individual bigwig tracks
## keys are used as regex pattern against filename/ 'track' key in dict
## more specific patterns should come first in dict as we break loop after first match
bigwig_specific = {
    "methyl|WGBS": {
        "viewLimits": "0:100",
        "maxHeightPixels": "500:30:8"},
    "male.*H3K9me3": {
        "viewLimits": "0:15"},
    "H3K9me3": {
        "viewLimits": "0:8"},
    "H3K27me3": {
        "viewLimits": "0:14"},
    "snRNA": {
        "viewLimits": "0:15",
        "maxHeightPixels": "500:30:8",
        "transformFunc": "LOG"},
    "RNA": {
        "viewLimits": "0:30"},
    }


composite_default = {"track": None,
                     "parent": None,
                     "type": None,
                     "compositeTrack": "on",
                     "shortLabel": None,
                     "longLabel": None,
                     "visibility": "hide",
                     "priority": 1,
                     "centerLabelsDense": "on",
                     "html": "examplePage"
}    


super_default = {"track": None,
                 "superTrack": "on show",
                 "parent": None,
                 "shortLabel": None,
                 "longLabel": None,
                 "priority": 1,
                 "html": "examplePage"
}

## adapted from http://code.activestate.com/recipes/577879-create-a-nested-dictionary-from-oswalk/
def get_directory_structure(rootdir):
    """
    Creates a nested dictionary that represents the folder structure of rootdir
    """
    rootdir = rootdir.rstrip(os.sep)
    start = rootdir.rfind(os.sep) + 1
    dir= {"containers": [rootdir]} 
    for path, dirs, files in os.walk(rootdir):
        folders = path[start:].split(os.sep)

        subdir = dict.fromkeys(files)
        parent = functools.reduce(dict.get, folders[:-1], dir)
        
        config = get_container_config(path, folders, subdir)
        
        parent[folders[-1]] = {'containers': dirs}
        parent[folders[-1]].update(config)
        
    return dir


def get_container_config(path, parents, files):
    """
    Creates a trackhub container and tracks config based on current 
    directory name and file content
    
    Allowed dir names: *.multiwig
                       *.composite
                       *.super
    
    Only one nesting level is supported by UCSC! So only super containers can hold one 
    level of multiwig or composite containers, not more! 
    Composite/multiwig containers can only contain tracks.
    
    The toplevel can contain tracks (not related to any container)
    
    Current code only supports *.bw|*.bigwig or *.bb|*.bigbed tracks!
    """    
    config = { 'tracks': {} }
    container_config = {}    
    generatorType = None
    
    if re.match(".*\.multiwig$",parents[-1],re.IGNORECASE):
        container_config.update(multiwig_default)
        container_config.update(bigwig_combined)
        generatorType = "multiwig"
    elif re.match(".*\.composite$",parents[-1],re.IGNORECASE):
        container_config.update(composite_default)
        generatorType = "composite"
    elif re.match(".*\.super$",parents[-1],re.IGNORECASE):
        container_config.update(super_default)
        generatorType = "super"
    elif len(parents)>1: 
        sys.exit("Every subdir needs to be a multiwig, composite or super container!")
        
    container_config["track"] = parents[-1]
    container_config["shortLabel"] = parents[-1]
    container_config["longLabel"] = parents[-1]
    
    container_config["parent"] = parents[len(parents)-2]
    
    if generatorType == "multiwig":
        for pat in bigwig_specific:
            if re.match(".*("+pat+")",container_config["track"],re.IGNORECASE):
                print(" ".join(["match ",container_config["track"]," ",pat]))
                container_config.update(bigwig_specific[pat])
                break
            
    ## toplevel must not have a parent entry
    if  len(parents)-2 <= 0:
        container_config.pop('parent',None)
    
    config['tracks'][parents[-1]] = container_config
    
    ## get per track config
    tracks = get_tracks_config(files, generatorType, parents)
    config['tracks'].update(tracks)
    
    ## set type for specific containers 
    if generatorType == "composite" or generatorType == "multiwig":
        tmp = set([v for k in tracks.keys() for kk,v in tracks[k].items() if kk == 'type'])
        multi_track_type = 'bigWig'
        if len(tmp)>0:
            multi_track_type = tmp.pop()
        if len(tmp) > 0:
            sys.exit("Only one tracktype allowed in composite or multiwig containers!")
        if generatorType == "multiwig" and multi_track_type != "bigWig":
            sys.exit("Only bigWig tracks are allowed in multiwig containers!")
        
        config['tracks'][parents[-1]]['type'] = multi_track_type
    
    ## update configs from config files if found (first *.yaml) in current path
    config = update_config_from_file(path, config)

    ## just dump config of current container into its directory
    ## file can be used as starting point to modify/add specific options
    with open(os.path.join(path,"container_config.used"), 'w') as f:
       yaml.dump(config['tracks'], f, default_flow_style=False)
    
    return config


## configure tracks for current container
def get_tracks_config(files, type, parents):
    """
    Creates a config per track from 'files'
    'type' is current container type
    'parents' is used to get the path and right parent name etc 
    
    Current code only supports *.bw|*.bigwig or *.bb|*.bigbed tracks!
    
    """
    tracks_config = {}
    global trackCounter
    
    for track_file in files:
        track_config = {}
        ## we have a bigwig file
        if re.match(".*\.(bw|bigwig)$",track_file,re.IGNORECASE):
            track_config.update(bigwig_default)
            if type != "multiwig":
                track_config.update(bigwig_combined)

            ## toplevel tracks have no parent entry
            if len(parents)-2 > -1:
                track_config["parent"] = parents[-1]
            else:
                track_config.pop('parent',None)
            
            track_config["track"] = "_".join(["track",str(trackCounter)])
            track_config["bigDataUrl"] = os.path.join(*parents[1:]+[track_file])
            track_config["shortLabel"] = track_file
            track_config["longLabel"] = track_file
            track_config["color"] = get_bigwig_color(track_file,parents[-1])
            trackCounter += 1
            
            if type != "multiwig":
                for pat in bigwig_specific:
                    if re.match(".*("+pat+")",track_file,re.IGNORECASE):
                        print(" ".join(["match ",track_file," ",pat]))
                        track_config.update(bigwig_specific[pat])
                        break
            
            tracks_config[track_file] = track_config
        ## we have a bigbed file
        elif re.match(".*\.(bb|bigbed)$",track_file,re.IGNORECASE):
            track_config.update(bigbed_default)
            if len(parents)-2 > -1:
                track_config["parent"] = parents[-1]
            else:
                track_config.pop('parent',None)
            
            track_config["track"] = "_".join(["track",str(trackCounter)])
            track_config["bigDataUrl"] = os.path.join(*parents[1:]+[track_file])
            track_config["shortLabel"] = track_file
            track_config["longLabel"] = track_file
            track_config["color"] = get_bigwig_color(track_file,parents[-1],bigbed_default['color'])
            trackCounter += 1
            
            tracks_config[track_file] = track_config
        
    return tracks_config


def get_bigwig_color(filename, parent, default="255,0,0"):
    #print([filename,parent])
    for pattern,color in bigwig_colors.items():
        if (re.search(pattern, filename, re.IGNORECASE) or re.search(pattern, parent, re.IGNORECASE)):
            #print(["match",pattern])
            return color
    return default


def update_config_from_file(path, config):
    
    ## take first yaml file that is found in path
    config_files = glob.glob(os.path.join(path,"*.yaml"))
    if config_files:
        config_file = config_files[0]
    else:
        return config 
    
    configFromFile = {}
    if os.path.isfile(config_file):
        with open(config_file, "r") as f:
            configFromFile = yaml.load(f)
     
    for tr in config['tracks']:
        if tr in configFromFile:
            config['tracks'][tr].update(configFromFile[tr])
    
    return config


def write_hub(file, hub, depth, in_root, outdir):
    
    ## write out container config section
    for container in hub['containers']:
        
        if depth>0:
            for k,v in hub[container]['tracks'][container].items():
                file.write("{m: <{de}}".format(m='',de=str((depth-1)*5))) 
                file.write("{} {}\n".format(k,v))
            file.write("\n")
            
        write_hub(file, hub[container], depth+1, in_root, outdir)
        
        ## write out all 'child' tracks of container
        for track in hub[container]['tracks']:
            if container != track:
                for k,v in hub[container]['tracks'][track].items():
                    if k=='bigDataUrl':
                        v = hub[container]['tracks'][track][k].split(os.sep)[-1]
                    ## indentation
                    file.write("{m: <{de}}".format(m='',de=str(depth*5)))
                    file.write("{} {}\n".format(k,v))
                ## remove link if we have on old link with same name
                if os.path.islink(os.path.join(os.path.abspath(outdir),track)):
                    os.remove(os.path.join(outdir,track))
                ## link track into output dir
                os.symlink(os.path.realpath(os.path.join(in_root,hub[container]['tracks'][track]['bigDataUrl'])), os.path.join(outdir,track))
                file.write("\n")

def main():

    global trackCounter

    parser = argparse.ArgumentParser() 
    
    parser.add_argument("indir",
                        help="input directory")

    parser.add_argument("-o", "--outputDir",
                        dest="outdir",
                        required=True,
                        help="output directory")
    
    parser.add_argument("-t", "--trackDbFilename",
                        dest="trackDbFilename",
                        default="trackDb.txt",
                        help="filename of trackhub config, useful if you use multiple trackDb files "
                        "for one organism, see also --postContent! (default: '%(default)s')")
    
    
    parser.add_argument("-i", "--startIndex",
                        dest="startIndex",
                        default=1,
                        type=int,
                        help="numerical index for first track, important if multiple trackDb files are used (default: '%(default)s')")
    
    parser.add_argument("-p", "--postContent",
                        dest="postContent",
                        default='',
                        help="string/text that is inserted at the end of generated trackDb, use eg. 'include trackDb.test.txt' to include an additional track config file; Note: you likely need to specify -t -i when you generate 'trackDb.test.txt' ! (default: '%(default)s')")
        
    args = parser.parse_args()
        
    if args.trackDbFilename == '':
        sys.exit("Please provide a filename for -t") 
        
    print(args.outdir)
    print(os.path.abspath(args.indir))
    
    trackCounter = args.startIndex
    
    ## get the hub by parsing directory structure
    hub = get_directory_structure(args.indir)

    os.makedirs(args.outdir,exist_ok = True)
    
    ## write hub config to output dir and link all files for upload into it
    with open(os.path.join(args.outdir,args.trackDbFilename), 'w') as f:
        write_hub(f,hub,0, args.indir,args.outdir)
        f.write(args.postContent)
        f.close()
    
    ## just dump hub as yaml file for inspection/debugging
    with open(os.path.join(args.outdir,args.trackDbFilename+".hub_dict.yaml"), 'w') as f:
       yaml.dump(hub, f, default_flow_style=False)
        

if __name__ == "__main__":
    main()
