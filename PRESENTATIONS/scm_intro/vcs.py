#!/usr/bin/env python
import os
import errno
import sys

layers=11

slide_template = r"""\begin{frame}[plain]{%(title)s}

  \centering
  
  \includegraphics[height=.95\textheight]{output/%(file)s}


  
\end{frame}

"""

def make_slides(basename, title, slides):
    outfile = file(basename+".tex", "w")

    try:
        os.mkdir("output")
    except os.error as e:
        if e.errno != errno.EEXIST:
            raise

    for i in range(slides):
        filename = basename+`i`
        outfile.write(slide_template%{
                "title" : title,
                "file"  : filename+".pdf"
                })
    
        os.system("dia --show-layers=%(layers)s --filter=eps --export=output/%(file)s.eps %(basename)s.dia"%
                  {"file" : filename,
                   "basename" : basename,
                   "layers" : ",".join(["Layer%d"%j for j in range(i+1)])})
        os.system("convert output/%(file)s.eps output/%(file)s.pdf"%{"file":filename})

sets = {
    "vcs" : [11, "Centralised version control"],
    "branching" : [7, "Making a feature branch in a centralised system"],
    "branch_merging" : [4, "Merging back feature branch in a centralised system"],
    "dvcs" : [8, "Distributed version control systems"],
    "pull_request"  : [4, "Pull request and code review"],
    "pull_request_2"  : [3, "Pull request and code review"]
    }

try:
    set = sys.argv[1]
    if set == "--print_targets":
        
        print " ".join([name+".tex" for name in sets.keys()])

    else:
        slides, title = sets[set]

        make_slides(set, title, slides)
except IndexError:
    for set, (slides, title) in sets.iter_items():
        make_slides(set, title, slides)
    
