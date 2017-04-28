import argparse
import collections
import json
import sys

def load_remove_comments( f):
    s= ""
    for line in f:
        idx= line.find('//')   # find C++ comment
        if idx == -1:  # no comment
            s = s + line + '\n'
        else:
            s = s + line[:idx] + '\n' # strip comment
    return s
        
def get_dict( dic, keylist):
    for key in keylist:
        dic= dic[key];
    return dic;

def set_dict( dic, keylist, val):
    for key in keylist[:-1]: 
        dic= dic.setdefault(key, {}) # dive into recursively, except for the rightmost key
    dic[ keylist[-1] ]= val 


def get_converted_key( conv, keystring): # keystring, possibly starting with '.'
    keylist= keystring.strip('.').split('.')
    
    try:    
        convkey= get_dict( conv, keylist)
    except:
        return keylist
    else:
        print "converted", '"'+keystring.strip('.')+'"', "to", '"'+convkey+'"'
        return convkey.strip('.').split('.')


def convert( din, dconv, dout, keystring):
    for (key,val) in din.items():
        if isinstance(val, dict): # recursion step
            convert( val, dconv, dout, keystring+'.'+key)
        else:
            convKey= get_converted_key( dconv, keystring+'.'+key) 
            set_dict( dout, convKey, val)
                    

def main():
    prog = 'python convert_json.py'
    description = 'A simple tool to convert JSON files with a prescribed key conversion.'
    parser = argparse.ArgumentParser(prog=prog, description=description)
    parser.add_argument('infile', nargs='?', type=argparse.FileType(),
                        help='a JSON file to be converted')
    parser.add_argument('convfile', nargs='?', type=argparse.FileType(),
                        help='a JSON file containing the conversion of JSON keys')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                        help='write the output to outfile')
    parser.add_argument('--sort-keys', action='store_true', default=False,
                        help='sort the output of dictionaries alphabetically by key')
    options = parser.parse_args()
    infile = options.infile or sys.stdin
    convfile = options.convfile
    outfile = options.outfile or sys.stdout
    sort_keys = options.sort_keys

    # read files
    with infile:
        try:
            print "Reading input file..."
            sinfile= load_remove_comments( infile)
            din = json.loads(sinfile,
                            object_pairs_hook=collections.OrderedDict)

        except ValueError as e:
            raise SystemExit(e)
            
    with convfile:
        try:
            print "Reading conversion file..."
            sconvfile= load_remove_comments( convfile)
            dconv= json.loads(sconvfile)
            
        except ValueError as e:
            raise SystemExit(e)
    
    # do conversion
    dout= collections.OrderedDict()
    print "Converting..."
    convert( din, dconv, dout, "")

    # write output
    print "Writing output file..."
    with outfile:
        json.dump(dout, outfile, sort_keys=sort_keys, indent=4)
        outfile.write('\n')

if __name__ == '__main__':
    main()

