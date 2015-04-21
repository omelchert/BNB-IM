""" main_bnb_sg_f.py

all-in-one file containing classes and function definitions for a simple
branch and bound algorithm to compute groundstate spin configurations
for Ising models with local fields for arbitrary bond configurations,
neighbor-relations and field strengths. 

Contains class defintions for Edge and UndirectedMultiGraph. Notation
consistent with 

"Some Basic Definitions in Graph Theory",
John W. Essam and Michael E. Fisher,
Review of Modern Physics, 42 (1970) 271

Branch-and-bound algorith as detailed in 

"Exact ground state of Ising systems",
S. Kobe, A. Hartwig,
Computer Physics Communications 16 (1978) 1-4

Author : Oliver Melchert
Date   : 04/20/2015 
"""
import copy


class Edge(object):
        """unordered pair of vertices

        EF1970 1.5: unordered pair of distinct vertices from  
        vertex set V. [i,j] is said to be incident with vertices
        i,j and connects them.
        """
        def __init__(self,i=None,j=None):
                """default constructor for a new instance of class Edge

                \prop i vertex from vertex set V
                \prop j vertex from vertex set V
                """
                self.i = i              # unordered pair of vertices
                self.j = j
                self._id   = None       # integer id for edge bookkeeping

        @property
        def id(self):
                return self._id

        @property
        def isLoop(self):
                return True if self.i==self.j else False

        def otherEnd(self,this):
                if this in [self.i,self.j]:
                        return self.i if this==self.j else self.j
                else: 
                        return None

        def __repr__(self):
                """string representation of edge instance"""
                return "[%d,%d].%s"%(self.i,self.j,str(self.id) if self.id!=None else 'X')


class UndirectedMultiGraph(object):
        '''undirected multi graph G = (V,E) without loops

        Contains class defintions for Edge and UndirectedMultiGraph. Notation
        consistent with 

        "Some Basic Definitions in Graph Theory",
        John W. Essam and Michael E. Fisher,
        Review of Modern Physics, 42 (1970) 271

        EF1970 1.12: undirected graph G = (V,E) is an vertex set V, having at
        least one member, together with an associated edge set E. The term 
        multigraph may be used to emphasize that multiedges are allowed.

        Author : Oliver Melchert
        Date   : Feb. 05, 2015 
        '''
	def __init__(self):
		"""default constructor for a new instance of class myGraph"""
		self._nVertices = 0	# number of vertices 
		self._nEdges    = 0  	# number of edges 
                self._eCtr      = 0     # counter for integer edgeIds
		self._incEdgId  = {}    # unordered set of node:(edgeId,otherEnd) pairs 
                self._E         = {}    # unordered set of edgeId:edge pairs 

        @property	
	def V(self): 
		"""returns vertex set of the graph, see EF1970 1.1"""
		return self._incEdgId.keys()

        @property
        def E(self):
                """returns set of individual edges of the graph, see EF1970 1.7"""
                return self._E.values()

        @property
	def v(self): 
                """returns number of vertices in graph, see EF1970 1.1"""
		return self._nVertices

        @property
	def e(self): 
                """returns number of edges in graph, see EF1970 1.7"""
		return self._nEdges

	def adjVertices(self,i):
                """returns set of vertices j adjacent to vertex j, see EF1970 1.31"""
		return set([otherEnd for (eId,otherEnd) in self._incEdgId[i].iteritems()]) 

	def incEdges(self,i):
                """returns set of edges [i,j] incident to vertex i, see EF1970 1.5"""
		return set([self._E[eId] for (eId,otherEnd) in self._incEdgId[i].iteritems()]) 

	def deg(self,i):
                """degree of a vertex, see EF1970 1.19"""
		return len(self._incEdgId[i])

	def addVertex(self,i):
		if i not in self.V: 
			self._incEdgId[i]={}
			self._nVertices += 1

	def delVertex(self,i):
                for eij in self.incEdges(i):
                        self.delEdge(eij)
                del self._incEdgId[i]
                self._nVertices -= 1

	def addEdge(self,e):
                if not e.isLoop:
                   self.addVertex(e.i)
                   self.addVertex(e.j)
                   e._id = self._eCtr
                   self._eCtr+=1
                   self._E[e._id] = e
                   self._incEdgId[e.j][e._id]=e.i
                   self._incEdgId[e.i][e._id]=e.j
                   self._nEdges += 1

	def delEdge(self,e):
                del self._incEdgId[e.i][e._id]
                del self._incEdgId[e.j][e._id]
                del self._E[e._id]
                self._nEdges -= 1

        def contains(self,e):
                """True if graph contains edge.
                
                Comparison is made based on the terminal vertices, not the edge id
                """
                return True if (e.j in self._incEdgId[e.i].values()) else False

        def __str__(self):
                """return string with graph in lemon graph format"""
                string = '@nodes\nlabel deg\n'
                for i in self.V:
                        string += '%3d %3d\n'%(i,self.deg(i))

                string += '@edges\n     label\n'
                for e in self.E:
                        string += '%3d  %3d  %3d\n'%(e.i,e.j,e.id)

                return string


def readGraph_bondList_field(fName):
        """read bond configuration and field-strengths from file 
        in DIMACS like format
        """
        G = UndirectedMultiGraph()
        h = dict()
        for line in open(fName):
               c = line.split()
               if c[0]=='h' and len(c)==3:
                       si      = int(c[1])
                       hi      = int(c[2])
                       h[si]   = hi
               if c[0]=='e' and len(c)==4:
                       eij     = Edge(int(c[1]), int(c[2]))
                       eij.wgt = int(c[3])
                       G.addEdge(eij)
        return G,h 


def energy_f(G,h,spin):
        """energy for given bond, spin, and field setup
        
        NOTE:
        -# edges need to have a ".wgt" attribute

        INPUT:
        \param[in] G Undirected weighted graph containing bond setup
        \param[in] h Local fields at spin locations 
        \param[in] spin Spin configuration 
        
        RETURNS: (E)
        \param[ret] E Energy of spin configuration
        """
        locErg  = [spin[e.i]*spin[e.j]*e.wgt for e in G.E]
        locPiv  = [h[i]*spin[i] for i in range(G.v)]
        return -sum(locErg)-sum(locPiv)


def bnb_sg_field_bounds(G,h):
        """compute upper and lower bound on configurational energy for
        spinglass with local fields for use in branch and bound algorithm

        Implemented according to KH1978:

        "Exact ground state of Ising systems",
        S. Kobe, A. Hartwig,
        Computer Physics Communications 16 (1978) 1-4

        see section 3 "Improved branch and bound technique" of the article.

        NOTE:
        -# edges need to have a ".wgt" attribute
        -# Code written at London Heathrow Airport while waiting for a connecting
           flight to Hannover. Hence, it might be a bit flakey and would benefit 
           from further refactoring

        INPUT:
        \param[in] G Undirected weighted graph containing bond setup
        \param[in] h Local fields at spin locations 
        
        RETURNS: (EMin,(sMax,EMax))
        \param[ret] EMin Lower bound for having all bonds satisfied and pivots maximised
        \param[ret] EMax Upper bound resulting from heuristic approach
        \param[ret] sMax Spin configuration with EMax, resulting from heuristic approach
        """
        # Lower bound for having all bonds satisfied and 
        # pivots maximised, see eq. (5) of KH1978
        EMin=-sum([abs(e.wgt) for e in G.E])-sum([abs(hi) for hi in h.values()])
        # Upper bound resulting from heuristic solution,
        # see algorithm description in section 3
        (sMax,EMax) = bnb_sg_field_heuristic(G,h)
        return EMin,(sMax,EMax)


def bnb_sg_field_heuristic(G,h):
        """heuristic approach to find spin configuration with low energy
        for spinglass with arbitraty bond configuration, neighborhood relations 
        and local fields.

        Implemented according to KH1978:

        "Exact ground state of Ising systems",
        S. Kobe, A. Hartwig,
        Computer Physics Communications 16 (1978) 1-4

        see section 3 "Improved branch and bound technique" of the article.

        NOTE:
        -# Code written at London Heathrow Airport while waiting for a connecting
           flight to Hannover. Hence, it might be a bit flakey and would benefit 
           from further refactoring

        INPUT:
        \param[in] G Undirected weighted graph containing bond setup
        \param[in] h Local fields at spin locations 
        
        RETURNS: (sMax,EMax)
        \param[ret] EMax Upper bound resulting from heuristic approach
        \param[ret] sMax Spin configuration resulting from heuristic approach

        AUTHOR : O. Melchert
        DATE   : 04/19/2015
        """
        # initialize empty spin array
        s    = [0 for i in range(G.v)]
        # determine id of spin with maximal local field
        i    = max([(abs(h[i]),i) for i in h.keys()])[1]
        # set spin so as to maximize the "pivot" si*hi
        s[i] = 1 if h[i]>0 else -1

        
        spinSet = set(h.keys())-set([i])
        for ii in range(G.v-1):
                sId     = None
                hLocMax = 0
                for k in spinSet:
                        hLoc = -sum([ e.wgt*s[e.otherEnd(k)] for e in G.incEdges(k)]) - h[k]
                        if abs(hLoc)>=abs(hLocMax):
                                hLocMax = hLoc
                                sId     = k
                s[sId]=1 if hLocMax<0 else -1
                spinSet -= set([sId])

        return s,energy_f(G,h,s)


def bnb_sg_field(i,G,h,sIni,EMin,EMax,res={}):
        """branch and bound recursion algorithm for Ising model with arbitrary
        bond configuration, neighbor-relations and local fields. Implemented 
        according to KH1978:

        "Exact ground state of Ising systems",
        S. Kobe, A. Hartwig,
        Computer Physics Communications 16 (1978) 1-4

        NOTE:
        -# For a simple O(N^2) algorithm to yield heuristic spin cfg with low
           energy, see KH1978, section 3 "Improved branch and bound technique" 
        -# To obtain all cfgs with energy <= EMax, disable (2.1) 
        -# Code written at London Heathrow Airport while waiting for a connecting
           flight to Hannover. Hence, it might be a bit flakey and would benefit 
           from further refactoring

        INPUT:
        \param[in] i Current spin id in recursion process
        \param[in] G Undirected weighted graph containing bond setup
        \param[in] h Local fields at spin locations 
        \param[in] sIni Spin configuration at current recursion depth 
        \param[in] EMin Lower bound that might at best be reached in current branch 
        \param[in,out] EMax Upper bound that should not be exeeded
        \param[in,out] res Dictionary of (E,[cfgs satisfying given bounds])-pairs 
        """
        # start brancing procedure at recursion depth i
        for si in [-1,1]:

                sIni[i]=si
                
                # change in current lower best bound due to new unsatisfied bonds
                dE_b = sum([ 2*abs(e.wgt) for e in G.incEdges(i) if e.wgt*si*sIni[e.otherEnd(i)] < 0])
                # change in current lower best bound due to negative pivots
                dE_f = 2*abs(h[i]) if h[i]*si < 0 else 0
                # update best bound for current branch at given current recurstion depth
                ECurr = EMin + dE_b + dE_f
                
                if ECurr <= EMax:
                   
                   # (1) further recursion step if current bounds are satisfied  
                   if i < G.v-1:
                        res,EMax = bnb_sg_field(i+1,G,h,sIni,ECurr,EMax,res)
                   
                   # (2) bookkeeping at maximal depth, i.e. for "full" cfgs
                   else:
                        # (2.1) update heuristic upper bound, if valid cfg
                        #       with lower energy is encountered (see last 
                        #       paragraph in section 2 of KH1978)
                        EMax=ECurr
                        # (2.2) update dictionary of resulting ECurr,sCurr pairs
                        ECurr = EMin+dE_b+dE_f
                        if ECurr not in res: res[ECurr]=[]
                        res[ECurr].append(copy.deepcopy(sIni))

                # clear spin before returning to previous recurstion depth
                sIni[i]=0

        return res,EMax


def main_bnb_sg_f():

        # read graph from specified bond-file
        bondFile = "./example_KH1978.cfg"
        #bondFile = "./example_HR6.5.cfg"
        G,h = readGraph_bondList_field(bondFile)

        # obtain bounds for branch-and-bound recursion
        EMin,(sMax,EMax) = bnb_sg_field_bounds(G,h)

        # set initial "empty" spin configuration and 
        # kick off branch-and-bound recursion
        sIni = [0 for i in range(G.v)]
        res,EMax  = bnb_sg_field(0,G,h,sIni,EMin,EMax)

        print "# BRANCH-AND-BOUND ENUMERATION OF SPIN CONFIGURATIONS TO"
        print "# YIELD GROUNDSTATE OF ISING MODELS WITH BOND CONFIGURATION"
        print "# AND LOCAL FIELDS AS SPECIFIED IN BOND-FILE:", bondFile
        print
        print "# LOWER BOUND: EMin = ", EMin
        print "# UPPER BOUND: EMax = ", EMax
        print 
        print "# ALL CFGS WITH E<=EMax (FIRST LISTED = GROUNDSTATE):"
        for E,cfgList in sorted(res.iteritems()):
                print
                print "# E =",E,":"
                print "   # (idx): (cfg)"
                for idx,cfg in enumerate(cfgList):
                        print "   %3d:  %s"%(idx+1,''.join(['+' if si > 0 else '-' for si in cfg]))


main_bnb_sg_f()
# EOF: main_bnb_sg_f.py
