import gc
import depend as dp
import codecs

import pdb
def flatten(l, ltypes=(list)):
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)
    

def node_id(arg, useport=False):
   if isinstance(arg,int):     return ('o%d' % arg)
   if isinstance(arg,dp.depend_array):
      if useport:  return ('o%d:o%d' % (id(arg._bval),id(arg))).replace('-', '_')
      else: return ('o%d' % (id(arg._bval))).replace('-', '_')
   else:                       return ('o%d' % id(arg)).replace('-', '_')

def depsame(d1, d2):
   if d1 is d2: return True
   if isinstance(d1,dp.depend_array) and isinstance(d2,dp.depend_array) and d1._bval is d2._bval: return True
   return False

def plot_deps(filename="dep_dump.dot"):
   gc.collect()
   objects = gc.get_objects()   
   
   # find objects which are depend objects   
   deps = [ o for o in objects if isinstance(o,dp.depend_base)]
   depids = [ id(o) for o in deps ]  
               
   # find edges
   dd_edges = []
   for o in deps:
      for do in o._dependants:
         dd_edges.append((o,do))      
   print "Found ", len(deps), " dependency objects, and ", len(dd_edges), " dep-dep edges"
      
   # remove multiple edges / edges between deparrays sharing the same storage
   tmp=[]
   for e in dd_edges:
      rmedge=False   
      for t in tmp:
         if depsame(e[0],t[0]) and depsame(e[1], t[1]): 
            rmedge=True; break           
      if not rmedge: tmp.append(e)
   dd_edges=tmp
   print len(dd_edges), " dep-dep edges left after cleanup"
   
   # find and refine sync edges
   sy_edges=[]
   for o1 in deps:
      if not o1._synchro is None:
         for o2 in o1._synchro.synced.values():
            sy_edges.append((o1,o2))
   tmp=[]
   for e in sy_edges:
      rmedge=False   
      for t in tmp:
         if depsame(e[0],t[0]) and depsame(e[1], t[1]): 
            rmedge=True; break
      if depsame(e[0],e[1]): rmedge=True             
      if not rmedge: tmp.append(e)
   sy_edges=tmp


   
   # find objects which contain depend objects
   depnames = {}
   conts = []
   contids = []
   cd_edges = []
   for o in objects:
      hasdep=False
      if hasattr(o,"__dict__"):
         for k,d in o.__dict__.items():
            if id(d) in depids:
               hasdep=True
               depnames[id(d)]=k
               cd_edges.append((o,d))
      if hasdep: 
         conts.append(o); contids.append(id(o))
               
   print "Found ", len(conts), " objects containing deps."
   
   cc_edges = []
   for o in conts:   
      for o2 in conts: 
         for k,d in o.__dict__.items():
            if id(o2) == id(d) :
               cc_edges.append((o,o2,k))
            if isinstance(d, list):   #also search lists for connections
               for k2 in d:
                  if id(o2) == id(k2) :
                     cc_edges.append((o,o2,k))
               

   
   dfile=codecs.open(filename, 'w', encoding='utf-8')
   dfile.write('digraph DepGraph {\n'
            '  node[shape=box, style=filled, fillcolor=white];\n')

   stores={}
   for o in deps: 
      if o._name is None or o._name == "":  name=node_id(o)
      else:  name=o._name         
      autofunc=(not o._func is None)
      autosync=(not o._synchro is None)                     
      if autosync: style='fillcolor="#40A080",fontcolor=black'
      elif autofunc: style='fillcolor="#C05010",fontcolor=black'
      else: style='fillcolor="#B0B0B0",fontcolor=black'         


      if isinstance(o,dp.depend_array):
         if not id(o._bval) in stores: stores[id(o._bval)]=([],(name,style))
         stores[id(o._bval)][0].append(o)
      else:
         dfile.write('  %s[shape=oval,label=%s,%s];\n' % (node_id(o),name,style))         
         
   for s,r in stores.items():
      dfile.write('  %s[shape=record,style="rounded,filled",%s,label="{' % (node_id(s),r[1][1]))
      for l in r[0]:
         dfile.write('<%s>%s' % (node_id(id(l)), node_id(id(l)) if (l._name is None or l._name == "" ) else l._name ) )
         if not l is r[0][-1]: dfile.write(' | ')
         else: dfile.write(' }"')
      dfile.write('];\n')
         
   for o in dd_edges:   
      dfile.write('  %s -> %s [color="#400000"];\n' % (node_id(o[0]),node_id(o[1])))      

   for o in sy_edges:   
      dfile.write('  %s -> %s [color="#006020"];\n' % (node_id(o[0]),node_id(o[1])))      

   for o in conts:   
      dfile.write('  %s[fontcolor=gray,label="%s (%d)"];\n' % (node_id(o), type(o).__name__, id(o)))

   for o in cd_edges:   
      dfile.write('  %s -> %s [color=gray];\n' % (node_id(o[0],useport=True), node_id(o[1],useport=True)))

   for o in cc_edges:
      dfile.write('  %s -> %s [color=black, label="%s"];\n' % (node_id(o[0]), node_id(o[1]), o[2]))             
   

   dfile.write('}'); dfile.close()
